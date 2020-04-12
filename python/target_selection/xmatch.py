#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-04-06
# @Filename: xmatch.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import inspect
import itertools
import os
import re

import networkx
import numpy
import peewee
import yaml
from peewee import fn

from sdssdb.connection import PeeweeDatabaseConnection
from sdssdb.utils import get_row_count

import target_selection
from target_selection.utils import (Timer, copy_pandas,
                                    sql_apply_pm, sql_iauname)


EPOCH = 'J2015.5'


class Catalog(peewee.Model):
    """Model for the output table."""

    catalogid = peewee.BigIntegerField(null=False)
    iauname = peewee.TextField(null=True)
    ra = peewee.DoubleField(null=False)
    dec = peewee.DoubleField(null=False)
    pmra = peewee.FloatField(null=True)
    pmdec = peewee.FloatField(null=True)
    parallax = peewee.FloatField(null=True)
    version = peewee.TextField(null=False)

    class Meta:
        primary_key = False


def get_relational_model(model):

    assert model._meta.primary_key != '__composite_key__', \
        'invalid pk for model {model__name__!r}.'

    class Meta:
        database = model._meta.database
        table_name = 'catalog_to_' + model._meta.table_name
        schema = model._meta.schema
        primary_key = False

    RelationalModel = type(
        'CatalogTo' + model.__name__, (peewee.Model,),
        {
            'catalogid': peewee.BigIntegerField(),
            model._meta.table_name + '_id': model._meta.primary_key.__class__(),
            'best': peewee.BooleanField(null=True, default=True),
            'Meta': Meta
        })

    return RelationalModel


def XMatchModel(Model, resolution=None, ra_column=None, dec_column=None,
                pmra_column=None, pmdec_column=None, is_pmra_cos=True,
                parallax_column=None, epoch_column=None, epoch=None,
                has_duplicates=False, skip_unlinked=False):
    """Expands the model `peewee:Metadata` with cross-matching parameters.

    The parameters defined can be accessed with the same name as
    ``Model._meta.<parameter>`` (e.g., ``Model._meta.has_duplicates``).

    Parameters
    ----------
    resolution : float
        The spatial resolution of the catalogue, in arcsec.
    ra_column : str
        The name of the right ascension column. If not provided, an attempt
        will be made to recover it from the Q3C index, if it exists. Assumed
        to be in degrees.
    dec_column : str
        As ``ra``, for the declination column.
    pmra_column : str
        The RA proper motion column, assumed to be in milliarcseconds per year.
    pmdec_column : str
        As ``pmra_column`` for the declination proper motion.
    is_pmra_cos : bool
        Whether ``pmra_column`` provides the RA proper motion corrected from
        declination(``pmra * cos(dec)``) or not.
    parallax_column : str
        The column containing the parallax, assumed to be in arcsec.
    epoch_column : str
        The column containing the epoch of the target coordinates, assumed to
        be in Julian date.
    epoch : float
        The Julian date of the targets. Applies to all the records in the
        table. ``epoch`` and ``epoch_column`` are mutually exclusive. If
        neither ``epoch_column`` nor ``epoch`` are defined, assumes that the
        epoch is ``J2015.5``.
    has_duplicates : bool
        Whether the table contains duplicates.
    skip_unlinked : bool
        If True, any row that cannot be directly linked to a record in the
        output table will be skipped. This option is ignored if the table
        is the first one to be processed.

    Returns
    -------
    Model
        The same input model with the additional cross-matching parameters
        added to the metadata.

    """

    meta = Model._meta

    meta.resolution = resolution or numpy.nan

    if not ra_column or not dec_column:
        indexes = meta.database.get_indexes(meta.table_name, meta.schema)
        for index in indexes:
            if 'q3c' in index.sql.lower():
                match = re.match(r'.+q3c_ang2ipix\("*(\w+)"*, "*(\w+)"*\).+',
                                 index.sql)
                if match:
                    ra_column, dec_column = match.groups()

    meta.ra_column = ra_column
    meta.dec_column = dec_column

    meta.pmra_column = pmra_column
    meta.pmdec_column = pmdec_column
    meta.is_pmra_cos = is_pmra_cos
    meta.parallax_column = parallax_column

    epoch = 'J2015.5' if (not epoch) & (not epoch_column) else epoch

    assert (((epoch is None) & (epoch_column is None)) or
            ((epoch is not None) ^ (epoch_column is not None))), \
        'epoch and epoch_column are mutually exclusive.'

    meta.epoch = epoch
    meta.epoch_column = epoch_column

    meta.has_duplicates = has_duplicates
    meta.skip_unlinked = skip_unlinked

    meta.approximate_row_count = get_row_count(meta.database, meta.table_name,
                                               schema=meta.schema, approximate=True)

    return Model


class XMatchPlanner(object):
    """Prepares and runs the catalogue cross-matching.

    Parameters
    ----------
    database : ~sdssdb.connection.PeeweeDatabaseConnection
        A ~sdssdb.connection.PeeweeDatabaseConnection to the database
        containing ``catalogdb``.
    models : list
        The list of `.XMatchModel` classes to be cross-matched.
    version : str
        The cross-matching version.

    """

    def __init__(self, database, models, version, order='hierarchical',
                 key='row_count', start_node=None, output_table='catalogdb.catalog',
                 log_path='./xmatch_{version}.log', debug=False):

        self.log = target_selection.log

        if log_path:
            if self.log.fh:
                self.log.removeHandler(self.log.fh)
                self.log.fh = None
            self.log.start_file_logger(log_path.format(version=version))

        if debug is True:
            self.log.sh.setLevel(0)
        elif debug is False:
            self.log.sh.setLevel(100)
        else:
            self.log.sh.setLevel(debug)

        self.database = database
        assert self.database.connected, 'database is not connected.'

        self.models = {model._meta.table_name: model for model in models}
        self.version = version

        # Sets the metadata of the Catalog table.
        if '.' in output_table:
            Catalog._meta.schema, Catalog._meta.table_name = output_table.split('.')
        else:
            Catalog._meta.table_name = output_table

        self.database.bind([Catalog])

        if Catalog.table_exists() and self.database.get_primary_keys(Catalog._meta.table_name,
                                                                     Catalog._meta.schema):
            self.log.warning('output table has a primary key. This is not needed and will '
                             'slow down insertions. Consider dropping the pk and '
                             'recreating it later.')

        self.relational_models = []

        self._check_models()
        self.update_model_graph()

        self.process_order = []
        self._set_process_order(order=order, key=key, start_node=start_node)

    @classmethod
    def read(cls, base_model, version, config_file=None, **kwargs):
        """Instantiates `.XMatchPlanner` from a configuration file.

        Parameters
        ----------
        base_model
            A PeeWee model from which all the models to use subclass, or a
            `~sdssdb.connection.PeeweeDatabaseConnection` to the database
            containing the models. In the latter case, the models must have
            been imported so that they are available in the ``models``
            attribute.
        version : str
            The cross-matching version.
        config_file : str
            The path to the configuration file to use. Defaults to
            ``config/xmatch.yml``. The file must contain a hash with the
            cross-match version.
        kwargs : dict
            User arguments that will override the configuration file values.

        """

        # TODO: I don't like too much this way of defining the models.
        if inspect.isclass(base_model) and issubclass(base_model, peewee.Model):
            database = base_model._meta.database
            models = set(base_model.__subclasses__())
            while True:
                old_models = models.copy()
                for model in old_models:
                    models |= set(model.__subclasses__())
                if models == old_models:
                    break
        elif isinstance(base_model, PeeweeDatabaseConnection):
            database = base_model
            models = database.models
        else:
            raise TypeError(f'invalid input of type {type(base_model)!r}')

        assert database.connected, 'database is not connected.'

        if config_file is None:
            config_file = os.path.dirname(target_selection.__file__) + '/config/xmatch.yml'

        config = yaml.load(open(config_file, 'r'), Loader=yaml.SafeLoader)

        assert version in config, f'version {version!r} not found in configuration.'
        config = config[version]

        table_config = config.pop('tables', {}) or {}
        include = config.pop('include', [])
        exclude = config.pop('exclude', [])

        assert 'schema' in config, 'schema is required in configuration.'
        schema = config.pop('schema')
        all_models = [model for model in database.models
                      if model._meta.schema == schema]

        if include:
            all_models = [model for model in all_models
                          if model._meta.table_name in include]
        if exclude:
            all_models = [model for model in all_models
                          if model._meta.table_name not in exclude]

        xmatch_models = []
        for model in all_models:
            table_name = model._meta.table_name
            table_params = table_config[table_name] if table_name in table_config else {}
            xmatch_models.append(XMatchModel(model, **table_params))

        xmodel_table_names = set([model._meta.table_name for model in xmatch_models])
        missing_tables = set(table_config.keys()) - xmodel_table_names
        if len(missing_tables) > 0:
            raise ValueError('some tables for which configuration was provided '
                             f'are not in the list of models: {missing_tables}')

        config.update(kwargs)

        return cls(database, xmatch_models, version, **config)

    def run(self):
        """Runs the cross-matching process."""

        # Make sure the output table exists.
        Catalog.create_table()

        for table_name in self.process_order:
            self.process_model(self.models[table_name])
            self.update_model_graph(silent=True)
            break

    def process_model(self, model):
        """Processes a model, loading it into the output table."""

        table_name = model._meta.table_name

        self.log.info(f'[{table_name.upper()}] Processing table {table_name}.')

        if table_name == self.process_order[0]:
            assert not model._meta.has_duplicates, 'first model to ingest cannot have duplicates.'
            is_first_model = True

        if is_first_model:
            return self._ingest(*self._build_select(model))

    def _build_select(self, model):
        """Returns the ModelSelect and associated Catalog fields."""

        meta = model._meta
        fields = meta.fields

        model_fields = []
        catalog_fields = []

        ra_field = fields[meta.ra_column]
        dec_field = fields[meta.dec_column]

        if meta.pmra_column:

            pmra_field = fields[meta.pmra_column]
            pmdec_field = fields[meta.pmdec_column]

            model_fields.extend([pmra_field, pmdec_field])
            catalog_fields.extend([Catalog.pmra, Catalog.pmdec])

            if (meta.epoch and meta.epoch != EPOCH) or meta.epoch_column:

                if meta.epoch:
                    # TODO: Is this correct?
                    delta_years = float(EPOCH[1:]) - float(meta.epoch[1:])
                else:
                    epoch_model = fn.ltrim(fields[meta.epoch_column], 'J').cast('REAL')
                    delta_years = float(EPOCH[1:]) - epoch_model

                ra_field, dec_field = sql_apply_pm(ra_field, dec_field,
                                                   pmra_field, pmdec_field,
                                                   delta_years,
                                                   is_pmra_cos=meta.is_pmra_cos)

        model_fields.extend([ra_field, dec_field])
        catalog_fields.extend([Catalog.ra, Catalog.dec])

        if meta.parallax_column:
            model_fields.append(fields[meta.parallax_column])
            catalog_fields.append(Catalog.parallax)

        model_fields.append(sql_iauname(ra_field, dec_field))
        catalog_fields.append(Catalog.iauname)

        return model.select(*model_fields), catalog_fields

    def _ingest(self, query, fields, chunk_size=1000000):
        """Ingests a query into the output table with coordinate conversion."""

        table_name = query.model._meta.table_name
        pk = query.model._meta.primary_key

        header = f'[{table_name.upper()}] '

        with Timer() as timer:

            with self.database.atomic():

                # Get the max catalogid currently in the table.
                max_id = Catalog.select(fn.MAX(Catalog.catalogid)).scalar() or 0
                self.log.debug(header + f'catalogid will start at {max_id+1}.')

                # Add the version and a row_number() function that sets the serial
                # value for catalogid. The ordering of the query happens here so no need
                # to order the select itself.
                query = query.select_extend(
                    self.version,
                    fn.row_number().over(order_by=pk) + max_id)
                fields.extend([Catalog.version, Catalog.catalogid])

                insert_query = Catalog.insert_many(query, fields).returning(Catalog.catalogid)
                self.log.debug(header + f'Running INSERT query: {insert_query}')

                catalogids = insert_query.tuples().execute(self.database)

        self.log.debug(header + f'Inserted {len(catalogids)} rows '
                       f'in {timer.interval:.3f} s.')

        RelationalModel = get_relational_model(query.model)
        RelationalModel.create_table()

        self.log.debug(header + 'Created relational table '
                       f'{RelationalModel._meta.table_name!r}')

        self.log.debug(header + f'Populating {RelationalModel._meta.table_name!r}')

        with Timer() as timer:

            # Keep query but now return only the pk of the queried table.
            query._returning = (pk,)
            query = query.order_by(pk)
            pks = query.tuples()

            for range0 in range(0, len(catalogids), chunk_size):

                pks_chunk = itertools.islice(pks, range0, range0 + chunk_size)
                cids_chunk = itertools.islice(catalogids, range0, range0 + chunk_size)

                data = {'cids': (cid[0] for cid in cids_chunk),
                        'pks': (pk[0] for pk in pks_chunk),
                        'best': True}

                copy_pandas(data, self.database,
                            RelationalModel._meta.table_name,
                            schema=RelationalModel._meta.schema,
                            columns=['catalogid',
                                     table_name + '_id',
                                     'best'])

        self.log.debug(header + f'Table {RelationalModel._meta.table_name!r} loading '
                       f'took {timer.interval:.3f} s.')

        return

    def _check_models(self):
        """Checks the input models."""

        non_existing_models = []
        for model in self.models.values():
            if not model.table_exists():
                self.log.warning(f'table {model._meta.table_name!r} does not exist.')
                non_existing_models.append(model)
            elif model._meta.table_name == Catalog._meta.table_name:
                self.log.debug(f'skipping output table {Catalog._meta.table_name!r}')
                non_existing_models.append(model)
            else:
                meta = model._meta
                self.log.debug(f'added table {meta.schema}.{meta.table_name} '
                               f'with resolution={meta.resolution:.2f} arcsec and '
                               f'row count={meta.approximate_row_count:.0f} records.')

        [self.models.pop(model._meta.table_name) for model in non_existing_models]

    def update_model_graph(self, silent=False):
        """Updates the model graph using models as nodes and fks as edges."""

        self.model_graph = networkx.Graph()

        all_models = list(self.models.values()) + [Catalog] + self.relational_models

        for model in all_models:

            table_name = model._meta.table_name
            schema = model._meta.schema

            self.model_graph.add_node(table_name)

            fks = self.database.get_foreign_keys(table_name, schema)
            for fk in fks:
                self.model_graph.add_edge(table_name, fk.dest_table)

        return self.model_graph

    def _set_process_order(self, order='hierarchical', key='row_count',
                           start_node=None):
        """Sets and returns the order in which tables will be processed.

        See `.XMatchPlanner` for details on how the order is decided depending
        on the input parameters.

        """

        assert order in ['hierarchical', 'global'], f'invalid order value {order!r}.'
        self.log.info(f'processing order mode is {order!r}')

        if isinstance(order, (list, tuple)):
            self.log.info(f'processing order: {order}')
            return order

        assert key in ['row_count', 'resolution'], f'invalid order key {key}.'
        self.log.info(f'ordering key is {key!r}.')

        ordered_tables = []

        graph = self.model_graph.copy()

        if order == 'hierarchical':
            subgraphs = networkx.connected_components(graph)
        else:
            subgraphs = [node for node in subgraphs.nodes]

        for model in [Catalog] + self.relational_models:
            graph.remove_node(model._meta.table_name)

        subgraphs_ext = []
        for sg in subgraphs:
            if start_node and start_node in sg:
                # Last item in record is 0 for initial table, 1 for other.
                # This prioritises the initial table in a sort without
                # having to reverse.
                subgraphs_ext.append((sg, numpy.nan, 0))
            else:
                if key == 'row_count':
                    total_row_count = sum([self.models[tn]._meta.approximate_row_count
                                           for tn in sg])
                    # Use -total_row_count to avoid needing reverse order.
                    subgraphs_ext.append((sg, -total_row_count, 1))
                elif key == 'resolution':
                    resolution = [self.models[tn]._meta.resolution for tn in sg]
                    if all(numpy.isnan(resolution)):
                        min_resolution = numpy.nan
                    else:
                        min_resolution = numpy.nanmin(resolution)
                    subgraphs_ext.append((sg, -(min_resolution or numpy.nan), 1))

        subgraphs_ordered = list(zip(*sorted(subgraphs_ext,
                                             key=lambda x: (x[2], x[1], x[0]))))[0]

        if order == 'global':
            ordered_tables = subgraphs_ordered
        else:
            for sgo in subgraphs_ordered:
                sg_ext = []
                for table_name in sgo:
                    if start_node and start_node == table_name:
                        sg_ext.append((table_name, numpy.nan, 0))
                        continue
                    if key == 'row_count':
                        row_count = self.models[table_name]._meta.approximate_row_count
                        sg_ext.append((table_name, -row_count, 1))
                    elif key == 'resolution':
                        resolution = self.models[table_name]._meta.resolution
                        sg_ext.append((table_name, resolution, 1))
                # Use table name as second sorting order to use alphabetic order in
                # case of draw.
                sg_ordered = list(zip(*sorted(sg_ext, key=lambda x: (x[2], x[1], x[0]))))[0]
                ordered_tables.extend(sg_ordered)

        self.log.info(f'processing order: {ordered_tables}')
        self.process_order = ordered_tables

        return ordered_tables
