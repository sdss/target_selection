#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2023-05-09
# @Filename: sdss_id.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import re

from typing import TYPE_CHECKING

import networkx
import peewee

from sdssdb.peewee.sdss5db import catalogdb as cdb
from sdssdb.peewee.sdss5db import targetdb as tdb

from target_selection import log
from target_selection.exceptions import TargetSelectionError


if TYPE_CHECKING:
    from sdssdb.connection import PeeweeDatabaseConnection


class SDSS_ID_Input(peewee.Model):
    """Model for the sdss_id_input table."""

    catalogid = peewee.BigIntegerField(primary_key=True)

    class Meta:
        table_name = 'sdss_id_input'
        schema = 'sandbox'


class SDSS_ID_To_Catalogid(peewee.Model):
    """Model for the sdss_id_to_catalogid table."""

    catalogid = peewee.BigIntegerField()
    sdssid = peewee.BigIntegerField()

    class Meta:
        table_name = 'sdss_id_to_catalogid'
        schema = 'sandbox'
        primary_key = peewee.CompositeKey('sdss_id', 'catalogid')


def update_sdss_ids(database: PeeweeDatabaseConnection,
                    target_selection_plan: str | None = None,
                    cross_match_plan: str | None = None):
    """Updates sdss_ids."""

    log.info('Generating list of input targets.')
    # generate_inputs(database,
    #                 target_selection_plan=target_selection_plan,
    #                 cross_match_plan=cross_match_plan)

    graph = SchemaGraph(database, 'catalogdb')


def generate_inputs(database: PeeweeDatabaseConnection,
                    target_selection_plan: str | None = None,
                    cross_match_plan: str | None = None):
    """Generate the table with the list of input targets for which to generate sdss_ids."""

    if not target_selection_plan and not cross_match_plan:
        raise TargetSelectionError('target_selection_plan or cross_match_plan are needed.')

    if cross_match_plan and '%' in cross_match_plan:
        raise TargetSelectionError('Templates cannot be used with cross_match_plan.')

    SDSS_ID_Input._meta.set_database(database)

    if SDSS_ID_Input.table_exists():
        log.warning('sdss_id_input exists. Dropping it.')
        SDSS_ID_Input.drop_table()

    log.info('Creating sdss_id_input.')
    SDSS_ID_Input.create_table()

    log.info('Selecting targets to process into sdss_id_input.')
    with database.atomic():
        database.execute_sql('SET LOCAL enable_bitmapscan = false')
        database.execute_sql("SET LOCAL work_mem = '20GB'")

        if target_selection_plan is not None:
            query = (tdb.Target
                     .select(tdb.Target.catalogid)
                     .join(tdb.CartonToTarget)
                     .join(tdb.Carton)
                     .join(tdb.Version)
                     .join_from(tdb.Target, SDSS_ID_To_Catalogid,
                                on=(tdb.Target.catalogid == SDSS_ID_To_Catalogid.catalogid),
                                join_type=peewee.JOIN.LEFT_OUTER)
                     .where(SDSS_ID_To_Catalogid.catalogid.is_null())
                     .distinct(tdb.Target.catalogid))

            if '%' in target_selection_plan:
                query = query.where(tdb.Version.plan % target_selection_plan)
            else:
                query = query.where(tdb.Version.plan == target_selection_plan)

            SDSS_ID_Input.insert_from(query, fields=[SDSS_ID_Input.catalogid]).execute()

        if SDSS_ID_Input.select().count() == 0:
            raise TargetSelectionError('All targets already have sdss_ids assigned to them.')

        # Run a sanity check. All targets to process must correspond to the same cross-match.
        log.debug('Checking that all targets correspond to the same cross-match.')
        n_cross_match = (SDSS_ID_Input
                         .select(cdb.Catalog.version_id)
                         .join(cdb.Catalog, on=(SDSS_ID_Input.catalogid == cdb.Catalog.catalogid))
                         .group_by(cdb.Catalog.version_id)
                         .count())

        if n_cross_match > 1:
            raise TargetSelectionError('Multiple cross-match versions detected.')


class SchemaGraph:

    def __init__(self, database: PeeweeDatabaseConnection, schema: str, join_weight='table_size'):

        self.database = database

        self.schema = schema
        self.graph = networkx.Graph()
        self.join_weight = join_weight

        self.update_graph()

    def update_graph(self):
        """Updates the graph adding nodes and edges from the schema."""

        tables = self.database.get_tables(self.schema)
        relsize = self.get_relation_size() if self.join_weight == 'table_size' else None

        for table in tables:
            if table == 'catalog':
                continue

            node_attrs = {}
            if self.join_weight == 'table_size':
                node_attrs['nrows'] = relsize.get(table, 1)
            self.graph.add_node(table, **node_attrs)

        for otable, fks in self.get_foreign_keys().items():
            for (ocols, dtable, dcols) in fks:
                self.graph.add_edge(otable, dtable, columns=ocols, dest_columns=dcols)

        # catalog_to_X does not have FKs, so add those manually.
        for table in tables:
            if table.startswith('catalog_to'):
                dtable = table[11:]
                dtable_pk = self.database.get_primary_keys(dtable, self.schema)
                self.graph.add_edge(table,
                                    dtable,
                                    columns=['target_id'],
                                    dest_columns=dtable_pk)

        print(
            list(
                networkx.shortest_simple_paths(
                    self.graph,
                    'catalog_to_gaia_dr3_source',
                    'catalog_to_legacy_survey_dr10',
                    weight='nrows')))

    def get_foreign_keys(self):
        """Returns a dictionary of foreign keys in ``schema``.

        Returns
        -------
        fks
            A dictionary of all the foreign keys in ``schema`` keyed by table. The value
            is a list of tuples in which the first element is the column(s) in the table,
            the second is the referenced table, and the third is the referenced column(s).

        """

        fks_query = self.database.execute_sql(
            f"""
            SELECT conrelid::regclass AS table_from, conname, pg_get_constraintdef(oid)
            FROM   pg_constraint
            WHERE  contype IN ('f', 'p ')
            AND    connamespace = '{self.schema}'::regnamespace
            AND    pg_get_constraintdef(oid) LIKE 'FOREIGN KEY %%'
            ORDER  BY conrelid::regclass::text, contype DESC;
            """
        )

        fks = {}
        for fk in fks_query:
            table = fk[0]
            match = re.match(
                r'FOREIGN KEY \(([\w_,\ ]+)\) REFERENCES ([\w_]+)\(([\w_,\ ]+)\)',
                fk[2]
            )
            if match:
                if table not in fks:
                    fks[table] = []
                fks[table].append(
                    (
                        list(map(lambda x: x.strip(), match.group(1).split(','))),
                        match.group(2),
                        list(map(lambda x: x.strip(), match.group(3).split(','))),
                    )
                )

        return fks

    def get_relation_size(self, relation: str | None = None):
        """Returns the number of tuples of a relation in the database.

        Parameters
        ----------
        relation
            The name of the relation. Otherwise returns a dictionary of all relations.

        Returns
        -------
        relsize
            If ``relation=None`` returns a dictionary of relation name to number of tuples
            (which can be considered an estimate of the number of rows for tables).
            Otherwise returns the number of tuples for that relation.

        """

        relsize = self.database.execute_sql(
            f"""
            SELECT pg_class.relname, pg_class.reltuples
            FROM pg_class
            JOIN pg_catalog.pg_namespace n ON n.oid = pg_class.relnamespace
            WHERE n.nspname='{self.schema}';
            """
        )

        relsize_dict = dict(relsize)

        if relation is None:
            return relsize_dict

        return relsize_dict[relation]
