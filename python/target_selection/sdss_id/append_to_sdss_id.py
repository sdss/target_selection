#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Zach Way (zway1@gsu.edu)
# @Date: 2024-04-23
# @Filename: append_to_sdss_id.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os
from datetime import date

import peewee
from peewee import JOIN, fn

from sdssdb.peewee.sdss5db.catalogdb import Catalog, database
from sdssdb.peewee.sdss5db.targetdb import Target

from .create_catalogidx_to_catalogidy import (
    MetaXMatch,
    TempMatch,
    UniqueMatch,
    create_unique_from_region,
)


class TempCatalogidV21(peewee.Model):
    """Model for a temporary table of all catalogids within a version."""

    pk = peewee.PrimaryKeyField()
    catalogid21 = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "temp_catalogid_v21"


class TempCatalogidV25(peewee.Model):
    """Model for a temporary table of all catalogids within a version."""

    pk = peewee.PrimaryKeyField()
    catalogid25 = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "temp_catalogid_v25"


class TempCatalogidV31(peewee.Model):
    """Model for a temporary table of all catalogids within a version."""

    pk = peewee.PrimaryKeyField()
    catalogid31 = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "temp_catalogid_v31"


class SdssIdStacked(peewee.Model):
    """Model for catalogdb sdss_id_stacked"""

    sdss_id = peewee.BigAutoField(primary_key=True)
    catalogid21 = peewee.BigIntegerField()
    catalogid25 = peewee.BigIntegerField()
    catalogid31 = peewee.BigIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()
    last_updated = peewee.DateField()

    class Meta:
        database = database
        schema = "sandbox"  # "catalogdb"
        table_name = "sdss_id_stacked"


class SdssIdFlat(peewee.Model):
    """Model for catalogdb sdss_id_flat"""

    sdss_id = peewee.BigIntegerField()
    catalogid = peewee.BigIntegerField()
    version_id = peewee.SmallIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()
    n_associated = peewee.SmallIntegerField()
    ra_catalogid = peewee.DoubleField()
    dec_catalogid = peewee.DoubleField()
    rank = peewee.SmallIntegerField(null=False)
    pk = peewee.BigAutoField(primary_key=True)

    class Meta:
        database = database
        schema = "sandbox"  # "catalogdb"
        table_name = "sdss_id_flat"


class SdssIdStackedAddendum(peewee.Model):
    """Model for addendum to sdss_id_stacked"""

    catalogid21 = peewee.BigIntegerField(index=True)
    catalogid25 = peewee.BigIntegerField(index=True)
    catalogid31 = peewee.BigIntegerField(index=True)
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()
    last_updated = peewee.DateField()

    class Meta:
        database = database
        schema = "sandbox"


class SdssIdFlatAddendum(peewee.Model):
    """Model for addendum to sdss_id_flat"""

    sdss_id = peewee.BigIntegerField()
    catalogid = peewee.BigIntegerField()
    version_id = peewee.SmallIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()
    n_associated = peewee.SmallIntegerField(null=True)
    ra_catalogid = peewee.DoubleField(null=True)
    dec_catalogid = peewee.DoubleField(null=True)
    rank = peewee.SmallIntegerField(null=True)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "sdss_id_flat_addendum"


class AppendToTables:
    """Adds new objects to the sdss_id tables.

    This class automatically uses the correlation method MetaXMatch and
    then has methods to put these objects through the sdss_id process. This
    can be done for 1) a list of catalogids, 2) a catalog_to_? table in
    catalogdb, or 3) all of targetdb.target

    Examples of how to run this class for each case are included at the end of
    the file append_to_sdss_id.py

    Parameters
    ----------
    database : peewee:PostgresqlDatabase
        A `peewee:PostgresqlDatabase` to the database the tables to
        cross-match.
    individual_table : str
        The table to add to the crossmatch in the form of catalogdb.catalog_to_?
        For example individual_table="catalogdb.catalog_to_too_target"
    catalogid_list : list
        A list of catalogids from the database
    """

    def __init__(self, database, individual_table=None, catalogid_list=None):
        self.database = database
        self.individual_table = individual_table
        self.ind_table_clean = self.individual_table
        self.catalogid_list = catalogid_list
        self.dir_path = os.path.dirname(__file__)

        if self.catalogid_list is not None:
            config = {
                "version_ids_to_match": [21, 25, 31],
                "individual_xmatch_config": self.dir_path + "/config/individual_crossmatches.yml",
                "log_file": "catalogidx_to_catalogidy_from_list.log",
                "show_first": 20,
                "split_insert_nunmber": 100000,
                "database_options": {"enable_hashjoin": "false"},
                # "split_query": [['panstarrs1',522000000000000,5000000000000]],
                "catalogid_list": catalogid_list,
            }

        elif self.individual_table is not None:
            if "." in self.individual_table:
                self.ind_table_clean = self.individual_table.split(".")[-1]

            config = {
                "version_ids_to_match": [21, 25, 31],
                "individual_xmatch_config": self.dir_path + "/config/individual_crossmatches.yml",
                "log_file": f"catalogidx_to_catalogidy_{self.ind_table_clean}.log",
                "show_first": 20,
                "split_insert_nunmber": 100000,
                "database_options": {"enable_hashjoin": "false"},
                # "split_query": [['panstarrs1', 522000000000000, 5000000000000]],
                "individual_table": individual_table,
            }
        else:
            config = {
                "version_ids_to_match": [21, 25, 31],
                "individual_xmatch_config": self.dir_path + "/config/individual_crossmatches.yml",
                "log_file": "catalogidx_to_catalogidy_all_target.log",
                "show_first": 20,
                "split_insert_nunmber": 100000,
                "database_options": {"enable_hashjoin": "false"},
            }
        self.config = config

    def run_MetaXMatch(self, database):
        """This takes in an individual catalog_to_? table (labeled 'individual_table') or a
        list of catalogids and creates the catalogidx_to_catalogidy_?
        and catalogidx_to_catalogidy_?_unique tables.
        """
        metax = MetaXMatch(
            config_filename=None,
            database=database,
            from_yaml=False,
            from_dict=True,
            config_dict=self.config,
            outer_join_sdss_id=True,
        )

        metax.run()
        create_unique_from_region(metax.output_name)
        return metax.output_name

    def create_temp_catalogid_lists(self, database, output_name):
        """This method takes the newly created catalogidx_to_catalogidy table and
        creates three aggregate tables in the sandbox containing all of the possible
        catalogids to be matched. This is also where the process restricts the matching of
        catalogids to 1) those which are labeled a 'best' match and 2) catalogids which have an
        inter-crossmatch match on that catalogid's lead catalog"""

        TempMatch._meta.table_name = output_name
        UniqueMatch._meta.table_name = output_name + "_unique"

        v21_cid_query_x = (
            TempMatch.select(TempMatch.catalogidx.alias("catalogid21"))
            .join(
                Catalog,
                on=(
                    (Catalog.catalogid == TempMatch.catalogidx) & (Catalog.lead == TempMatch.lead)
                ),
            )
            .switch(TempMatch)
            .join(
                SdssIdStacked,
                join_type=JOIN.LEFT_OUTER,
                on=(TempMatch.catalogidx == SdssIdStacked.catalogid21),
            )
            .where((TempMatch.version_idx == 21) & (SdssIdStacked.catalogid21.is_null()))
        )

        v25_cid_query_x = (
            TempMatch.select(TempMatch.catalogidx.alias("catalogid25"))
            .join(
                Catalog,
                on=(
                    (Catalog.catalogid == TempMatch.catalogidx) & (Catalog.lead == TempMatch.lead)
                ),
            )
            .switch(TempMatch)
            .join(
                SdssIdStacked,
                join_type=JOIN.LEFT_OUTER,
                on=(TempMatch.catalogidx == SdssIdStacked.catalogid25),
            )
            .where((TempMatch.version_idx == 25) & (SdssIdStacked.catalogid25.is_null()))
        )

        v31_cid_query_x = (
            TempMatch.select(TempMatch.catalogidx.alias("catalogid31"))
            .join(
                Catalog,
                on=(
                    (Catalog.catalogid == TempMatch.catalogidx) & (Catalog.lead == TempMatch.lead)
                ),
            )
            .switch(TempMatch)
            .join(
                SdssIdStacked,
                join_type=JOIN.LEFT_OUTER,
                on=(TempMatch.catalogidx == SdssIdStacked.catalogid31),
            )
            .where((TempMatch.version_idx == 31) & (SdssIdStacked.catalogid31.is_null()))
        )

        v21_cid_query_y = (
            TempMatch.select(TempMatch.catalogidy.alias("catalogid21"))
            .join(
                Catalog,
                on=(
                    (Catalog.catalogid == TempMatch.catalogidx) & (Catalog.lead == TempMatch.lead)
                ),
            )
            .switch(TempMatch)
            .join(
                SdssIdStacked,
                join_type=JOIN.LEFT_OUTER,
                on=(TempMatch.catalogidy == SdssIdStacked.catalogid21),
            )
            .where((TempMatch.version_idy == 21) & (SdssIdStacked.catalogid21.is_null()))
        )

        v25_cid_query_y = (
            TempMatch.select(TempMatch.catalogidy.alias("catalogid25"))
            .join(
                Catalog,
                on=(
                    (Catalog.catalogid == TempMatch.catalogidx) & (Catalog.lead == TempMatch.lead)
                ),
            )
            .switch(TempMatch)
            .join(
                SdssIdStacked,
                join_type=JOIN.LEFT_OUTER,
                on=(TempMatch.catalogidy == SdssIdStacked.catalogid25),
            )
            .where((TempMatch.version_idy == 25) & (SdssIdStacked.catalogid25.is_null()))
        )

        v31_cid_query_y = (
            TempMatch.select(TempMatch.catalogidy.alias("catalogid31"))
            .join(
                Catalog,
                on=(
                    (Catalog.catalogid == TempMatch.catalogidx) & (Catalog.lead == TempMatch.lead)
                ),
            )
            .switch(TempMatch)
            .join(
                SdssIdStacked,
                join_type=JOIN.LEFT_OUTER,
                on=(TempMatch.catalogidy == SdssIdStacked.catalogid31),
            )
            .where((TempMatch.version_idy == 31) & (SdssIdStacked.catalogid31.is_null()))
        )

        for table in [TempCatalogidV21, TempCatalogidV25, TempCatalogidV31]:
            if table.table_exists():
                database.drop_tables([table])

        # self.database.create_tables([TempCatalogidV21, TempCatalogidV25, TempCatalogidV31])
        TempCatalogidV21.create_table()
        TempCatalogidV25.create_table()
        TempCatalogidV31.create_table()

        temp_catalodids_indexes = """ CREATE INDEX temp_catalogid_v21_catalogid21_idx ON
                                        sandbox.temp_catalogid_v21 (catalogid21);
                                      CREATE INDEX temp_catalogid_v25_catalogid25_idx ON
                                        sandbox.temp_catalogid_v25 (catalogid25);
                                      CREATE INDEX temp_catalogid_v31_catalogid31_idx ON
                                        sandbox.temp_catalogid_v31 (catalogid31); """
        self.database.execute_sql(temp_catalodids_indexes)

        TempCatalogidV21.insert_from(v21_cid_query_x, [TempCatalogidV21.catalogid21]).execute()
        TempCatalogidV21.insert_from(v21_cid_query_y, [TempCatalogidV21.catalogid21]).execute()
        TempCatalogidV25.insert_from(v25_cid_query_x, [TempCatalogidV25.catalogid25]).execute()
        TempCatalogidV25.insert_from(v25_cid_query_y, [TempCatalogidV25.catalogid25]).execute()
        TempCatalogidV31.insert_from(v31_cid_query_x, [TempCatalogidV31.catalogid31]).execute()
        TempCatalogidV31.insert_from(v31_cid_query_y, [TempCatalogidV31.catalogid31]).execute()

        if self.catalogid_list is not None:
            catid_input_query = (
                Catalog.select(Catalog.catalogid, Catalog.version_id)
                .join(
                    SdssIdFlat,
                    join_type=JOIN.LEFT_OUTER,
                    on=(SdssIdFlat.catalogid == Catalog.catalogid),
                )
                .where(Catalog.catalogid << self.catalogid_list)
                .where(SdssIdFlat.catalogid.is_null())
            )

            for row in catid_input_query:
                if row.version_id == 21:
                    TempCatalogidV21.insert(catalogid21=row.catalogid).execute()
                elif row.version_id == 25:
                    TempCatalogidV25.insert(catalogid25=row.catalogid).execute()
                elif row.version_id == 31:
                    TempCatalogidV31.insert(catalogid31=row.catalogid).execute()
        elif self.individual_table is not None:
            try:
                ind_table = self.database.models[self.individual_table]
            except:
                raise ValueError(
                    "Could not find the table in database model: " + self.individual_table
                )

            ind_table_input_query = (
                Catalog.select(Catalog.catalogid, Catalog.version_id)
                .join(ind_table, on=(Catalog.catalogid == ind_table.catalogid))
                .join(
                    SdssIdFlat,
                    join_type=JOIN.LEFT_OUTER,
                    on=(SdssIdFlat.catalogid == Catalog.catalogid),
                )
                .where(SdssIdFlat.catalogid.is_null())
            )

            for row in ind_table_input_query:
                if row.version_id == 21:
                    TempCatalogidV21.insert(catalogid21=row.catalogid).execute()
                elif row.version_id == 25:
                    TempCatalogidV25.insert(catalogid25=row.catalogid).execute()
                elif row.version_id == 31:
                    TempCatalogidV31.insert(catalogid31=row.catalogid).execute()
        else:
            target_input_query = (
                Catalog.select(Catalog.catalogid, Catalog.version_id)
                .join(Target, on=(Catalog.catalogid == Target.catalogid))
                .join(
                    SdssIdFlat,
                    join_type=JOIN.LEFT_OUTER,
                    on=(SdssIdFlat.catalogid == Catalog.catalogid),
                )
                .where(SdssIdFlat.catalogid.is_null())
            )

            for row in target_input_query:
                if row.version_id == 21:
                    TempCatalogidV21.insert(catalogid21=row.catalogid).execute()
                elif row.version_id == 25:
                    TempCatalogidV25.insert(catalogid25=row.catalogid).execute()
                elif row.version_id == 31:
                    TempCatalogidV31.insert(catalogid31=row.catalogid).execute()

    def create_sdss_id_stacked_addendum(self, database, output_name):
        """This method matches the previous catalogids with the 'sdss_id' process
        and places them in a temporary table in the sandbox (to be added later). At
        this point, the ra/dec_sdss_id are assigned, which are just the coordinates of
        the most recent catalogid in the match."""

        large_cte_query = f"""DROP TABLE IF EXISTS sandbox.sdss_id_stacked_addendum;
        CREATE TABLE sandbox.sdss_id_stacked_addendum as (
            with cat_ids21 as (select distinct tc.catalogid21, cat.lead from
                sandbox.temp_catalogid_v21 tc join catalog cat on cat.catalogid=tc.catalogid21),
            cat_ids25 as (select distinct tc.catalogid25, cat.lead from
                sandbox.temp_catalogid_v25 tc join catalog cat on cat.catalogid=tc.catalogid25),
            cat_ids31 as (select distinct tc.catalogid31, cat.lead from
                sandbox.temp_catalogid_v31 tc join catalog cat on cat.catalogid=tc.catalogid31),
            sq21_25 as (select cc.catalogidx, cc.catalogidy from sandbox.{output_name} cc
                    join catalog cat on cat.catalogid=cc.catalogidx and cat.lead=cc.lead
                    where cc.version_idx=21 and cc.version_idy=25),
            sq25_31 as (select cc.catalogidx, cc.catalogidy from sandbox.{output_name} cc
                    join catalog cat on cat.catalogid=cc.catalogidx and cat.lead=cc.lead
                    where cc.version_idx=25 and cc.version_idy=31),
            left21_to_25 as (select cat_ids21.catalogid21, sq21_25.catalogidy as catalogid25 from
                cat_ids21 left join sq21_25 on cat_ids21.catalogid21=sq21_25.catalogidx),
            add_outer25 as (select left21_to_25.catalogid21,
                cat_ids25.catalogid25 as catalogid25 from
                left21_to_25 full outer join cat_ids25 on
                left21_to_25.catalogid25=cat_ids25.catalogid25),
            left_25_to_31 as (select add_outer25.catalogid21,
                add_outer25.catalogid25, sq25_31.catalogidy as catalogid31 from
                add_outer25 left join sq25_31 on add_outer25.catalogid25=sq25_31.catalogidx),
            add_outer31 as (select left_25_to_31.catalogid21, left_25_to_31.catalogid25,
                cat_ids31.catalogid31 from
                left_25_to_31 full outer join cat_ids31 on
                left_25_to_31.catalogid31=cat_ids31.catalogid31)
            select * from add_outer31);"""
        self.database.execute_sql(large_cte_query)

        add_ra_dec_columns = """ ALTER TABLE sandbox.sdss_id_stacked_addendum
                                    ADD COLUMN ra_sdss_id double precision;
                                 ALTER TABLE sandbox.sdss_id_stacked_addendum
                                    ADD COLUMN dec_sdss_id double precision; """
        self.database.execute_sql(add_ra_dec_columns)

        add_last_updated_column = f""" ALTER TABLE sandbox.sdss_id_stacked_addendum
                                           ADD COLUMN last_updated DATE;
                                       UPDATE sandbox.sdss_id_stacked_addendum SET
                                           last_updated = '{str(date.today())}'"""
        self.database.execute_sql(add_last_updated_column)

        create_indexes = """ CREATE INDEX sdss_id_stacked_addendum_catalogid21_idx ON
                                sandbox.sdss_id_stacked_addendum (catalogid21);
                             CREATE INDEX sdss_id_stacked_addendum_catalogid25_idx ON
                                sandbox.sdss_id_stacked_addendum (catalogid25);
                             CREATE INDEX sdss_id_stacked_addendum_catalogid31_idx ON
                                sandbox.sdss_id_stacked_addendum (catalogid31); """
        self.database.execute_sql(create_indexes)

        SdssIdStackedAddendum._meta.table_name = "sdss_id_stacked_addendum"

        ra_dec_update31 = (
            SdssIdStackedAddendum.update(ra_sdss_id=Catalog.ra, dec_sdss_id=Catalog.dec)
            .from_(Catalog)
            .where(SdssIdStackedAddendum.catalogid31 == Catalog.catalogid)
            .where(SdssIdStackedAddendum.catalogid31.is_null(False))
        )

        ra_dec_update25 = (
            SdssIdStackedAddendum.update(ra_sdss_id=Catalog.ra, dec_sdss_id=Catalog.dec)
            .from_(Catalog)
            .where(SdssIdStackedAddendum.catalogid25 == Catalog.catalogid)
            .where(SdssIdStackedAddendum.catalogid31.is_null())
            .where(SdssIdStackedAddendum.catalogid25.is_null(False))
        )

        ra_dec_update21 = (
            SdssIdStackedAddendum.update(ra_sdss_id=Catalog.ra, dec_sdss_id=Catalog.dec)
            .from_(Catalog)
            .where(SdssIdStackedAddendum.catalogid21 == Catalog.catalogid)
            .where(SdssIdStackedAddendum.catalogid31.is_null())
            .where(SdssIdStackedAddendum.catalogid25.is_null())
            .where(SdssIdStackedAddendum.catalogid21.is_null(False))
        )

        ra_dec_update31.execute()
        ra_dec_update25.execute()
        ra_dec_update21.execute()

    def add_to_sdss_id_stacked(self, database):
        """This method adds the previously made sdss_id_stached_addendum to sdss_id_stacked"""

        sid_stacked_f = [
            SdssIdStacked.catalogid21,
            SdssIdStacked.catalogid25,
            SdssIdStacked.catalogid31,
            SdssIdStacked.ra_sdss_id,
            SdssIdStacked.dec_sdss_id,
            SdssIdStacked.last_updated,
        ]

        to_add = SdssIdStackedAddendum.select(
            SdssIdStackedAddendum.catalogid21,
            SdssIdStackedAddendum.catalogid25,
            SdssIdStackedAddendum.catalogid31,
            SdssIdStackedAddendum.ra_sdss_id,
            SdssIdStackedAddendum.dec_sdss_id,
            SdssIdStackedAddendum.last_updated,
        ).tuples()
        with self.database.atomic():
            SdssIdStacked.insert_many(to_add, fields=sid_stacked_f).execute()

    def create_sdss_id_flat_addendum(self, database):
        """This method takes the new sdss_id rows in sdss_id_stacked and
        adds them to a temporary table in the sandbox in the sdss_id_flat format.
        Here, the ra/dec_catalogid are populated with the specific values associated
        with the catalogid."""

        if SdssIdFlatAddendum.table_exists():
            self.database.drop_tables([SdssIdFlatAddendum])
            # self.database.create_tables([SdssIdFlatAddendum])
            SdssIdFlatAddendum.create_table()
        else:
            # self.database.create_tables([SdssIdFlatAddendum])
            SdssIdFlatAddendum.create_table()

        sid_flat_fields = [
            SdssIdFlatAddendum.sdss_id,
            SdssIdFlatAddendum.catalogid,
            SdssIdFlatAddendum.version_id,
            SdssIdFlatAddendum.ra_sdss_id,
            SdssIdFlatAddendum.dec_sdss_id,
        ]

        max_sdss_id = SdssIdFlat.select(fn.MAX(SdssIdFlat.sdss_id)).scalar()

        sid_flat_add_v21 = (
            SdssIdStacked.select(
                SdssIdStacked.sdss_id,
                SdssIdStacked.catalogid21,
                peewee.Value(21),
                SdssIdStacked.ra_sdss_id,
                SdssIdStacked.dec_sdss_id,
            )
            .where(SdssIdStacked.catalogid21.is_null(False))
            .where(SdssIdStacked.sdss_id > max_sdss_id)
            .tuples()
        )

        with self.database.atomic():
            SdssIdFlatAddendum.insert_many(sid_flat_add_v21, fields=sid_flat_fields).execute()

        sid_flat_add_v25 = (
            SdssIdStacked.select(
                SdssIdStacked.sdss_id,
                SdssIdStacked.catalogid25,
                peewee.Value(25),
                SdssIdStacked.ra_sdss_id,
                SdssIdStacked.dec_sdss_id,
            )
            .where(SdssIdStacked.catalogid25.is_null(False))
            .where(SdssIdStacked.sdss_id > max_sdss_id)
            .tuples()
        )

        with self.database.atomic():
            SdssIdFlatAddendum.insert_many(sid_flat_add_v25, fields=sid_flat_fields).execute()

        sid_flat_add_v31 = (
            SdssIdStacked.select(
                SdssIdStacked.sdss_id,
                SdssIdStacked.catalogid31,
                peewee.Value(31),
                SdssIdStacked.ra_sdss_id,
                SdssIdStacked.dec_sdss_id,
            )
            .where(SdssIdStacked.catalogid31.is_null(False))
            .where(SdssIdStacked.sdss_id > max_sdss_id)
            .tuples()
        )
        with self.database.atomic():
            SdssIdFlatAddendum.insert_many(sid_flat_add_v31, fields=sid_flat_fields).execute()

        cte = (
            SdssIdFlatAddendum.select(SdssIdFlatAddendum.catalogid, fn.COUNT("*").alias("ct"))
            .group_by(SdssIdFlatAddendum.catalogid)
            .cte("catid_n_associated", columns=("catalogid", "ct"))
        )

        (
            SdssIdFlatAddendum.update(n_associated=cte.c.ct)
            .from_(cte)
            .where(SdssIdFlatAddendum.catalogid == cte.c.catalogid)
            .with_cte(cte)
            .execute()
        )

        (
            SdssIdFlatAddendum.update(ra_catalogid=Catalog.ra, dec_catalogid=Catalog.dec)
            .from_(Catalog)
            .where(SdssIdFlatAddendum.catalogid == Catalog.catalogid)
            .execute()
        )

        ranked_values = SdssIdFlatAddendum.select(
            SdssIdFlatAddendum.catalogid,
            SdssIdFlatAddendum.sdss_id,
            fn.RANK()
            .over(
                partition_by=[SdssIdFlatAddendum.catalogid], order_by=[SdssIdFlatAddendum.sdss_id]
            )
            .alias("rank"),
        )

        (
            SdssIdFlatAddendum.update(rank=ranked_values.c.rank)
            .from_(ranked_values)
            .where(
                (SdssIdFlatAddendum.catalogid == ranked_values.c.catalogid)
                & (SdssIdFlatAddendum.sdss_id == ranked_values.c.sdss_id)
            )
            .execute()
        )

    def add_to_sdss_id_flat(self, database):
        """This method adds sdss_id_flat_addendum to sdss_id_flat"""

        sid_flat_f = [
            SdssIdFlat.sdss_id,
            SdssIdFlat.catalogid,
            SdssIdFlat.version_id,
            SdssIdFlat.ra_sdss_id,
            SdssIdFlat.dec_sdss_id,
            SdssIdFlat.n_associated,
            SdssIdFlat.ra_catalogid,
            SdssIdFlat.dec_catalogid,
            SdssIdFlat.rank,
        ]

        to_add = SdssIdFlatAddendum.select(
            SdssIdFlatAddendum.sdss_id,
            SdssIdFlatAddendum.catalogid,
            SdssIdFlatAddendum.version_id,
            SdssIdFlatAddendum.ra_sdss_id,
            SdssIdFlatAddendum.dec_sdss_id,
            SdssIdFlatAddendum.n_associated,
            SdssIdFlatAddendum.ra_catalogid,
            SdssIdFlatAddendum.dec_catalogid,
            SdssIdFlatAddendum.rank,
        ).tuples()

        with self.database.atomic():
            SdssIdFlat.insert_many(to_add, fields=sid_flat_f).execute()


"""
EXAMPLE OF ADDING A SERIES OF CATALOGIDS
foo = AppendToTables(database, catalogid_list=[list, of, integer, catalogids])
output_name = foo.run_MetaXMatch(database)
foo.create_temp_catalogid_lists(database, output_name)
foo.create_sdss_id_stacked_addendum(database, output_name)
foo.add_to_sdss_id_stacked(database)
foo.create_sdss_id_flat_addendum(database)
foo.add_to_sdss_id_flat(database)
"""

"""
EXAMPLE OF ADDING MISSING CATALOGIDS FROM A TABLE
foo = AppendToTables(database, individual_table='catalogdb.catalog_to_too_target')
output_name = foo.run_MetaXMatch(database)
foo.create_temp_catalogid_lists(database, output_name)
foo.create_sdss_id_stacked_addendum(database, output_name)
foo.add_to_sdss_id_stacked(database)
foo.create_sdss_id_flat_addendum(database)
foo.add_to_sdss_id_flat(database)
"""

"""
EXAMPLE OF ADDING MISSING CATALOGIDS FROM TARGETDB
foo = AppendToTables(database)
output_name = foo.run_MetaXMatch(database)
foo.create_temp_catalogid_lists(database, output_name)
foo.create_sdss_id_stacked_addendum(database, output_name)
foo.add_to_sdss_id_stacked(database)
foo.create_sdss_id_flat_addendum(database)
foo.add_to_sdss_id_flat(database)
"""
