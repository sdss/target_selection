#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-04-06
# @Filename: xmatch.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import inspect
import os
import re

import networkx
import numpy
import peewee
import yaml
from peewee import fn
from playhouse.migrate import PostgresqlMigrator, migrate

from sdssdb.connection import PeeweeDatabaseConnection
from sdssdb.utils import get_row_count
from sdsstools.color_print import color_text

import target_selection
from target_selection.utils import Timer, copy_pandas, get_epoch, sql_apply_pm


EPOCH = 2015.5


class Catalog(peewee.Model):
    """Model for the output table."""

    catalogid = peewee.BigIntegerField(primary_key=True, null=False)
    iauname = peewee.TextField(null=True)
    ra = peewee.DoubleField(null=False)
    dec = peewee.DoubleField(null=False)
    pmra = peewee.FloatField(null=True)
    pmdec = peewee.FloatField(null=True)
    parallax = peewee.FloatField(null=True)
    lead = peewee.TextField(null=False)
    version = peewee.TextField(null=False)

    class Meta:
        primary_key = False


def get_relational_model(model, prefix='catalog_to_'):

    meta = model._meta

    assert meta.primary_key != '__composite_key__', \
        f'composite pk fount for model {model.__name__!r}.'

    # Auto/Serial are automatically PKs. Convert them to integers to avoid
    # having two pks in the relational table.
    if meta.primary_key.__class__.field_type == 'AUTO':
        model_pk_class = peewee.IntegerField
    elif meta.primary_key.__class__.field_type == 'BIGAUTO':
        model_pk_class = peewee.BigIntegerField
    else:
        model_pk_class = meta.primary_key.__class__

    class BaseRelationalModel(peewee.Model):

        id = peewee.BigAutoField(primary_key=True)
        catalogid = peewee.BigIntegerField()
        target_id = model_pk_class()
        distance = peewee.DoubleField(null=True)
        best = peewee.BooleanField(null=False, default=True)

        class Meta:
            database = meta.database
            schema = meta.schema

    model_prefix = ''.join(x.capitalize() or '_' for x in prefix.rstrip().split('_'))

    RelationalModel = type(model_prefix + model.__name__, (BaseRelationalModel,), {})
    RelationalModel._meta.table_name = prefix + meta.table_name

    return RelationalModel


def add_fks(database, RelationalModel, related_model):
    """Adds foreign keys and indexes to a relational model class."""

    rtschema = RelationalModel._meta.schema
    rtname = RelationalModel._meta.table_name

    migrator = PostgresqlMigrator(database)

    with database.atomic():

        migrate(
            migrator.set_search_path(rtschema),
            migrator.add_index(rtname, ('catalogid',), False),
            migrator.add_constraint(
                rtname, 'catalogid',
                peewee.SQL('FOREIGN KEY (catalogid) '
                           f'REFERENCES {Catalog._meta.schema}.{Catalog._meta.table_name} '
                           '(catalogid)')),
            migrator.add_index(rtname, ('target_id',), False),
            migrator.add_constraint(
                rtname, 'target_id',
                peewee.SQL(
                    f'FOREIGN KEY (target_id) '
                    f'REFERENCES {related_model._meta.schema}.{related_model._meta.table_name} '
                    f'({related_model._meta.primary_key.column_name})'))
        )

    migrator.database.close()


def XMatchModel(Model, resolution=None, ra_column=None, dec_column=None,
                pmra_column=None, pmdec_column=None, is_pmra_cos=True,
                parallax_column=None, epoch_column=None, epoch=None,
                epoch_format='jyear', has_duplicates=False,
                skip_unlinked=False, skip=False):
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
        The column containing the epoch of the target coordinates.
    epoch : float
        The epoch of the targets that applies to all the records in the
        table. ``epoch`` and ``epoch_column`` are mutually exclusive. If
        neither ``epoch_column`` nor ``epoch`` are defined, assumes that the
        epoch is ``2015.5``.
    epoch_format : str
        The format of the epoch. Either Julian year (``'jyear'``) or Julian
        date (``'jd'``).
    has_duplicates : bool
        Whether the table contains duplicates.
    skip : bool
        If `True`, the table will be used as a join node but will not be
        cross-matched.
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

    class XMatchMeta:
        pass

    meta = Model._meta
    meta.xmatch = XMatchMeta()

    meta.xmatch.resolution = resolution or numpy.nan

    if not ra_column or not dec_column:
        indexes = meta.database.get_indexes(meta.table_name, meta.schema)
        for index in indexes:
            if 'q3c' in index.sql.lower():
                match = re.match(r'.+q3c_ang2ipix\("*(\w+)"*, "*(\w+)"*\).+',
                                 index.sql)
                if match:
                    ra_column, dec_column = match.groups()

    meta.xmatch.ra_column = ra_column
    meta.xmatch.dec_column = dec_column

    meta.xmatch.pmra_column = pmra_column
    meta.xmatch.pmdec_column = pmdec_column
    meta.xmatch.is_pmra_cos = is_pmra_cos
    meta.xmatch.parallax_column = parallax_column

    epoch = EPOCH if (not epoch) & (not epoch_column) else epoch

    assert (((epoch is None) & (epoch_column is None)) or
            ((epoch is not None) ^ (epoch_column is not None))), \
        'epoch and epoch_column are mutually exclusive.'

    meta.xmatch.epoch = epoch
    meta.xmatch.epoch_column = epoch_column
    meta.xmatch.epoch_format = epoch_format

    meta.xmatch.has_duplicates = has_duplicates
    meta.xmatch.skip_unlinked = skip_unlinked
    meta.xmatch.skip = skip

    meta.xmatch.row_count = get_row_count(meta.database, meta.table_name,
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

    def __init__(self, database, models, version, extra_nodes=[],
                 order='hierarchical', key='row_count', epoch=EPOCH,
                 query_radius=None, start_node=None, schema='catalogdb',
                 output_table='catalog', order_by='pk',
                 log_path='./xmatch_{version}.log', debug=False,
                 show_sql=False, sample_region=None):

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

        self.schema = schema
        self.database = database
        assert self.database.connected, 'database is not connected.'

        self.version = version

        # Sets the metadata of the Catalog table.
        Catalog._meta.schema = schema
        Catalog._meta.table_name = output_table
        self.database.bind([Catalog])

        self.models = {model._meta.table_name: model for model in models}

        self.extra_nodes = {model._meta.table_name: model for model in extra_nodes}

        self._check_models()
        self.update_model_graph()

        self.process_order = []
        self._set_process_order(order=order, key=key, start_node=start_node)

        self._options = {'order_by': order_by,
                         'query_radius': query_radius or 1 / 3600.,
                         'show_sql': show_sql,
                         'sample_region': sample_region,
                         'epoch': epoch}

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
            models = database.models.values()
        else:
            raise TypeError(f'invalid input of type {type(base_model)!r}')

        assert database.connected, 'database is not connected.'

        if config_file is None:
            config_file = os.path.dirname(target_selection.__file__) + '/config/xmatch.yml'

        config = yaml.load(open(config_file, 'r'), Loader=yaml.SafeLoader)

        assert version in config, f'version {version!r} not found in configuration.'
        config = config[version]

        table_config = config.pop('tables', {}) or {}
        exclude = config.pop('exclude', []) or []

        assert 'schema' in config, 'schema is required in configuration.'
        schema = config['schema']

        models = {model._meta.table_name: model for model in models
                  if (model.table_exists() and
                      model._meta.schema == schema and
                      model not in exclude)}

        xmatch_models = {}
        for table_name in table_config:
            table_params = table_config[table_name]
            xmatch_models[table_name] = XMatchModel(models[table_name],
                                                    **table_params)

        extra_nodes = [models[table_name] for table_name in models
                       if table_name not in xmatch_models]

        config.update(kwargs)
        return cls(database, xmatch_models.values(), version,
                   extra_nodes=extra_nodes, **config)

    def _check_models(self):
        """Checks the input models."""

        remove_models = []

        for model in self.models.values():

            meta = model._meta
            table_name = meta.table_name

            if not model.table_exists():
                self.log.warning(f'table {table_name!r} does not exist.')
            elif table_name == Catalog._meta.table_name:
                self.log.debug(f'skipping output table {Catalog._meta.table_name!r}.')
                self.extra_nodes[Catalog._meta.table_name] = Catalog
            elif table_name.startswith(Catalog._meta.table_name + '_to_'):
                self.log.debug(f'skipping relational table {meta.table_name!r} '
                               'and adding it as an extra node.')
                self.extra_nodes[table_name] = model
            else:
                self.log.debug(f'added table {meta.schema}.{table_name} '
                               f'with resolution={meta.xmatch.resolution:.2f} arcsec and '
                               f'row count={meta.xmatch.row_count:.0f} records.')
                continue

            remove_models.append(model)

        [self.models.pop(model._meta.table_name) for model in remove_models]
        [self.extra_nodes.pop(y) for y in [x for x in self.extra_nodes
                                           if not self.extra_nodes[x].table_exists()]]

        if len(self.models) == 0:
            raise RuntimeError('no models to cross-match.')

    def update_model_graph(self, silent=False):
        """Updates the model graph using models as nodes and fks as edges."""

        self.model_graph = networkx.Graph()

        all_models = list(self.models.values()) + list(self.extra_nodes.values())

        for model in all_models:
            table_name = model._meta.table_name
            self.model_graph.add_node(table_name, model=model)

        for model in all_models:
            table_name = model._meta.table_name
            for fk_model in model._meta.model_refs:
                ref_table_name = fk_model._meta.table_name
                if ref_table_name not in self.model_graph.nodes:
                    continue
                self.model_graph.add_edge(table_name, ref_table_name)

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

        for model in self.extra_nodes.values():
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
                    total_row_count = sum([self.models[tn]._meta.xmatch.row_count
                                           for tn in sg
                                           if not self.models[tn]._meta.xmatch.skip])
                    # Use -total_row_count to avoid needing reverse order.
                    subgraphs_ext.append((sg, -total_row_count, 1))
                elif key == 'resolution':
                    resolution = [self.models[tn]._meta.xmatch.resolution
                                  for tn in sg
                                  if not self.models[tn]._meta.xmatch.skip]
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
                    model = self.models[table_name]
                    if start_node and start_node == table_name:
                        sg_ext.append((table_name, numpy.nan, 0))
                        continue
                    elif model._meta.xmatch.skip:
                        continue
                    if key == 'row_count':
                        row_count = model._meta.xmatch.row_count
                        sg_ext.append((table_name, -row_count, 1))
                    elif key == 'resolution':
                        resolution = model._meta.xmatch.resolution
                        sg_ext.append((table_name, resolution, 1))
                # Use table name as second sorting order to use alphabetic order in
                # case of draw.
                sg_ordered = list(zip(*sorted(sg_ext, key=lambda x: (x[2], x[1], x[0]))))[0]
                ordered_tables.extend(sg_ordered)

        self.log.info(f'processing order: {ordered_tables}')
        self.process_order = ordered_tables

        return ordered_tables

    def get_simple_paths(self, source, dest, allow_relational_table=True,
                         return_models=False):
        """Determines all possible join path between two tables.

        Instead of returning all possible simple paths between the ``source``
        and ``dest`` tables, this algorithm follows a "scorched earth" approach
        in which once an edge has been used for a join it cannot be used again.
        This produces only distint joins between the two nodes.

        Parameters
        ----------
        source : str
            The initial table for the path.
        dest : str
            The final table for the path.
        allow_relational_tables : bool, optional
            If `False`, the path that contains only the initial and final
            tables joined via their relational table is not included.
        return_models : bool
            If `True`, returns each path as a list of models. Otherwise returns
            the table names.

        Returns
        -------
        list
            A list in which each item is a path with the table names or models
            joining ``source`` to ``dest``. The list is sorted in order of
            increasing path length.

        """

        graph = self.model_graph.copy()

        if not allow_relational_table:
            rel_table_name = Catalog._meta.table_name + '_to_' + source
            if rel_table_name in graph:
                graph.remove_node(rel_table_name)

        paths = []

        while True:
            try:
                shortest_path = networkx.algorithms.shortest_path(graph, source, dest)
                paths.append(shortest_path)
                for node in shortest_path[1:-1]:
                    graph.remove_node(node)
            except networkx.NetworkXNoPath:
                break

        if len(paths) == 0:
            return []

        if return_models:
            paths = [[self.model_graph.nodes[node]['model'] for node in path]
                     for path in paths]

        return paths

    def run(self):
        """Runs the cross-matching process."""

        # Make sure the output table exists.
        if not Catalog.table_exists():
            Catalog.create_table()
            self.log.info(f'Created table {Catalog._meta.table_name!r}.')

        self.extra_nodes[Catalog._meta.table_name] = Catalog
        self.update_model_graph(silent=True)

        if self._options['sample_region']:
            sample_region = self._options['sample_region']
            self.log.warning(f'Using sample region {sample_region!r}.')

        with Timer() as timer:
            for table_name in self.process_order:
                model = self.models[table_name]
                self.process_model(model)

        self.log.info(f'Cross-matching completed in {timer.interval:.3f} s.')

    def process_model(self, model):
        """Processes a model, loading it into the output table."""

        table_name = model._meta.table_name
        self.log.header = f'[{table_name.upper()}] '

        self.log.info(f'Processing table {table_name}.')

        if model._meta.xmatch.has_duplicates:
            raise ValueError('handling of tables with duplicates '
                             'is not implemented.')

        # Check if there are already records in catalog for this version.
        if Catalog.select().where(Catalog.version == self.version).limit(1).count() == 0:
            is_first_model = True
        else:
            is_first_model = False

        rel_model, rel_model_created = self._create_relational_model(model)
        rel_table_name = rel_model._meta.table_name
        if rel_model_created:
            self.log.debug(f'Created relational table {rel_table_name!r}.')
        else:
            self.log.debug(f'Relational table {rel_table_name!r} already existed.')

        if is_first_model:
            self._load_catalog(model, rel_model)
        else:
            self._run_phase_1(model, rel_model)
            self._run_phase_2(model, rel_model)
            self._run_phase_3(model, rel_model)

        if rel_model_created:
            self.log.debug(f'Adding indexes and FKs to {rel_model._meta.table_name!r}.')
            add_fks(self.database, rel_model, model)

        self.log.header = ''

    def _build_select(self, model):
        """Returns the ModelSelect with all the fields to populate Catalog."""

        meta = model._meta
        xmatch = meta.xmatch
        fields = meta.fields

        model_fields = [meta.primary_key.alias('id')]

        ra_field = fields[xmatch.ra_column]
        dec_field = fields[xmatch.dec_column]

        to_epoch = self._options['epoch']

        if xmatch.pmra_column:

            pmra_field = fields[xmatch.pmra_column]
            pmdec_field = fields[xmatch.pmdec_column]
            model_fields.extend([pmra_field.alias('pmra'),
                                 pmdec_field.alias('pmdec')])

            if (xmatch.epoch and xmatch.epoch != to_epoch) or xmatch.epoch_column:
                delta_years = to_epoch - get_epoch(model)
                ra_field, dec_field = sql_apply_pm(ra_field, dec_field,
                                                   pmra_field, pmdec_field,
                                                   delta_years,
                                                   is_pmra_cos=xmatch.is_pmra_cos)
                if not xmatch.is_pmra_cos:
                    pmra_field *= fn.cos(fn.radians(dec_field))

        else:

            model_fields.extend([peewee.Value(None).cast('REAL').alias('pmra'),
                                 peewee.Value(None).cast('REAL').alias('pmdec')])

        model_fields.extend([ra_field.alias('ra'), dec_field.alias('dec')])

        if xmatch.parallax_column:
            model_fields.append(fields[xmatch.parallax_column].alias('parallax'))
        else:
            model_fields.append(peewee.Value(None).cast('REAL').alias('parallax'))

        # model_fields.append(sql_iauname(ra_field, dec_field).alias('iauname'))

        query = model.select(*model_fields)

        if self._options['order_by'] == 'q3c':
            query = query.order_by(fn.q3c_ang2ipix(ra_field, dec_field))
        else:
            query = query.order_by(meta.primary_key.asc())

        if self._options['sample_region']:
            sample_region = self._options['sample_region']
            query = query.where(fn.q3c_radial_query(fields[xmatch.ra_column],
                                                    fields[xmatch.dec_column],
                                                    sample_region[0],
                                                    sample_region[1],
                                                    sample_region[2]))

        return query

    def _load_catalog(self, model, rel_model, query=None):
        """Ingests a query into the output and relational tables."""

        meta = model._meta
        table_name = meta.table_name

        query = query or self._build_select(model)

        with Timer() as timer:

            with self.database.atomic():

                # Get the max catalogid currently in the table.
                max_cid = Catalog.select(fn.MAX(Catalog.catalogid)).scalar() or 0

                x = query.cte('x')

                y = peewee.CTE('y',
                               rel_model.insert_from(
                                   x.select(fn.row_number().over() + max_cid, x.c.id,
                                            peewee.SQL('true'))
                                   .where(x.c.ra.is_null(False),  # Make sure RA and Dec exist.
                                          x.c.dec.is_null(False)),
                                   [rel_model.catalogid,
                                    rel_model.target_id,
                                    rel_model.best])
                               .returning(rel_model.catalogid,
                                          rel_model.target_id))

                insert_query = Catalog.insert_from(
                    y.select(y.c.catalogid,
                             x.c.ra,
                             x.c.dec,
                             x.c.pmra,
                             x.c.pmdec,
                             x.c.parallax,
                             peewee.Value(table_name),
                             peewee.Value(self.version))
                     .join(x, on=(x.c.id == y.c.target_id)),
                    [Catalog.catalogid,
                     Catalog.ra,
                     Catalog.dec,
                     Catalog.pmra,
                     Catalog.pmdec,
                     Catalog.parallax,
                     Catalog.lead,
                     Catalog.version]
                ).returning().with_cte(x, y)

                self.log.debug(f'Running INSERT query{self._get_sql(insert_query)}')
                n_rows = insert_query.execute()

        self.log.debug(f'Inserted {n_rows} rows in {timer.interval:.3f} s.')

        return

    def _create_relational_model(self, model):
        """Creates a relational table for a given model.

        Returns the relational model and `True` if the table was created or
        `False` if it already existed.

        """

        RelationalModel = get_relational_model(model, prefix=Catalog._meta.table_name + '_to_')
        rtname = RelationalModel._meta.table_name

        if not RelationalModel.table_exists():
            RelationalModel.create_table()
            created = True
        else:
            if (RelationalModel
                    .select()
                    .join(Catalog, on=(Catalog.catalogid == RelationalModel.catalogid))
                    .where(Catalog.version == self.version)
                    .limit(1).count()) > 0:
                raise RecursionError(f'Relational table {rtname!r} contains records '
                                     f'for this version of cross-matching ({self.version}). '
                                     'This problem needs to be solved manually.')

            if RelationalModel.select().count() > 0:
                self.log.warning(f'Relational table {rtname!r} is not empty!')

            created = False

        # Add foreign key fields here. We want to avoid Peewee creating them
        # as constraints and indexes if the table is created because that would
        # slow down inserts. We'll created them manually with add_fks.
        RelationalModel._meta.add_field('catalog',
                                        peewee.ForeignKeyField(Catalog,
                                                               column_name='catalogid',
                                                               backref='+'))
        RelationalModel._meta.add_field('target',
                                        peewee.ForeignKeyField(model,
                                                               column_name='target_id',
                                                               backref='+'))

        self.extra_nodes[rtname] = RelationalModel
        self.update_model_graph(silent=True)

        return RelationalModel, created

    def _build_join(self, path):
        """Returns a build query for a given join path."""

        model = path[0]
        meta = model._meta
        xmatch = meta.xmatch

        query = model.select()
        for node in path[1:]:
            query = query.join(node)

        if self._options['sample_region']:
            sample_region = self._options['sample_region']
            query = query.where(fn.q3c_radial_query(meta.fields[xmatch.ra_column],
                                                    meta.fields[xmatch.dec_column],
                                                    sample_region[0],
                                                    sample_region[1],
                                                    sample_region[2]))

        return query

    def _load_relational_table(self, rel_model, catalogids, target_ids):
        """Inserts rows into the relational table associated with a model.

        Parameters
        ----------
        rel_model : peewee:Model
            The relational model to load.
        catalogids : tuple
            The catalogids to be associated.
        target_ids : tuple or peewee:ModelSelect
            Either a tuple with the target_ids from the related model or a
            query that will return them.
        chunk_size : int
            How many records to insert at a time.

        Returns
        -------
        n_rows : int
            The number of rows inserted.

        """

        meta = rel_model._meta

        if isinstance(target_ids, peewee.ModelSelect):
            # Keep query but now return only the pk of the queried table.
            target_ids._returning = (target_ids.model._meta.primary_key,)
            target_ids = (tid[0] for tid in target_ids.tuples())

        data = {'cids': catalogids, 'tids': target_ids, 'best': True}

        df = copy_pandas(data, self.database, meta.table_name, schema=meta.schema,
                         columns=['catalogid', 'target_id', 'best'])

        return len(df)

    def _run_phase_1(self, model, rel_model):
        """Runs the linking against matched catalogids stage."""

        table_name = model._meta.table_name
        rel_table_name = rel_model._meta.table_name

        self.log.info('Phase 1: linking existing targets.')

        join_paths = self.get_simple_paths(table_name, Catalog._meta.table_name,
                                           allow_relational_table=False)

        if len(join_paths) == 0:
            self.log.debug(f'No paths found between {table_name!r} and output table.')
            return False

        self.log.debug(f'Found {len(join_paths)} paths between '
                       f'{table_name!r} and output table.')

        # Go through the list following the same order as process_order to ensure
        # that the cross-matching has already happened for those tables.
        paths_processed = []

        for through_table in self.process_order:
            # Loop over he paths that contain the through table. Because of the
            # scorched earth mode of the path-finding algorithm, there can only
            # be 0 or 1 paths that contain through_table.
            for path in [pth for pth in join_paths if (through_table in pth and
                                                       pth not in paths_processed)]:

                join_models = [self.model_graph.nodes[node]['model'] for node in path]

                with Timer() as timer:

                    with self.database.atomic():

                        query = self._build_join(join_models).select(
                            model._meta.primary_key, Catalog.catalogid,
                            peewee.Value(False))

                        insert_query = rel_model.insert_from(
                            query, fields=[rel_model.target_id,
                                           rel_model.catalogid,
                                           rel_model.best]).returning()

                        self.log.debug(f'Inserting linked targets into '
                                       f'{rel_table_name!r} with join path {path}'
                                       f'{self._get_sql(insert_query)}')

                        nids = insert_query.execute()

                self.log.debug(f'Inserted {nids} records in {timer.interval:.3f} s.')

                paths_processed.append(path)

        if len(paths_processed) == 0:
            self.log.warning('No join paths were found using table that '
                             'had already been cross-matched. You should '
                             'review your processing order.')
            return False

        return True

    def _run_phase_2(self, model, rel_model):
        """Associates existing targets in Catalog with entries in the model."""

        self.log.info('Phase 2: cross-matching against existing targets.')

        meta = model._meta
        xmatch = meta.xmatch

        table_name = meta.table_name
        rel_table_name = rel_model._meta.table_name

        model_pk = meta.primary_key
        model_ra = meta.fields[xmatch.ra_column]
        model_dec = meta.fields[xmatch.dec_column]

        model_epoch = get_epoch(model)
        catalog_epoch = self._options['epoch']

        self.log.debug('Determining maximum delta epoch between tables.')
        max_delta_epoch = peewee.Value(
            model.select(fn.MAX(fn.ABS(model_epoch - catalog_epoch))).scalar())

        # Determine what targets in the model haven't been matched in phase 1.
        unmatched = (Catalog.select().join(rel_model, peewee.JOIN.LEFT_OUTER)
                     .where(rel_model.catalogid.is_null()).cte('unmatched'))

        # Create a subquery that returns the entries in the model that are
        # within query_radius of the Catalog target.
        subq = (model.select(model_pk.alias('target_id'),
                             fn.q3c_dist_pm(unmatched.c.ra,
                                            unmatched.c.dec,
                                            unmatched.c.pmra,
                                            unmatched.c.pmdec,
                                            1,  # Catalog pmra is pmra*cos_delta by definition.
                                            catalog_epoch,
                                            model_ra,
                                            model_dec,
                                            model_epoch).alias('distance'))
                .where(fn.q3c_join_pm(unmatched.c.ra,
                                      unmatched.c.dec,
                                      unmatched.c.pmra,
                                      unmatched.c.pmdec,
                                      1,
                                      catalog_epoch,
                                      model_ra,
                                      model_dec,
                                      model_epoch,
                                      max_delta_epoch,
                                      self._options['query_radius'])))

        # Do a lateral join (for each loop). Partition over each group of
        # targets that are within radius of a catalogid. Mark the one with
        # the smallest distance to the Catalog target as best.
        partition = fn.first_value(subq.c.target_id).over(partition_by=[unmatched.c.catalogid],
                                                          order_by=[subq.c.distance.asc()])
        xmatches = (unmatched.select(unmatched.c.catalogid,
                                     subq.c.target_id,
                                     subq.c.distance,
                                     peewee.Value(partition == subq.c.target_id).alias('best'))
                    .join(subq, 'CROSS JOIN LATERAL')
                    .order_by(unmatched.c.catalogid)
                    .with_cte(unmatched))

        # Insert as a CTE because we want to gather some stats from the returned values.
        insert_cte = peewee.CTE('insert_cte',
                                rel_model.insert_from(
                                    xmatches,
                                    fields=[rel_model.catalogid,
                                            rel_model.target_id,
                                            rel_model.distance,
                                            rel_model.best])
                                .returning(rel_model.catalogid,
                                           rel_model.target_id))

        insert_query = (insert_cte
                        .select(
                            fn.count(insert_cte.c.catalogid.distinct()),
                            fn.count(insert_cte.c.target_id.distinct()))
                        .with_cte(insert_cte))

        with Timer() as timer:
            with self.database.atomic():

                self.log.debug('Inserting cross-matched targets into '
                               f'{rel_table_name!r}{self._get_sql(insert_query)}')

                n_catalogid, n_target_id = list(insert_query
                                                .tuples()
                                                .execute(self.database))[0]

        self.log.debug(f'Inserted {n_catalogid} records in {rel_table_name!r} '
                       f'associated with {n_target_id} targets in {table_name!r}, '
                       f'in {timer.interval:.3f} s.')

    def _run_phase_3(self, model, rel_model):
        """Add non-matched targets to Catalog and the relational table."""

        self.log.info('Phase 3: adding non cross-matched targets.')

        # Build the initial query.
        query = self._build_select(model)

        # Limit the query to targets not matched.
        query = (query
                 .join(rel_model, peewee.JOIN.LEFT_OUTER)
                 .where(rel_model.catalogid.is_null()))

        # Insert the records
        self._load_catalog(model, rel_model, query=query)

    def _get_sql(self, query):
        """Returns coulourised SQL text for logging."""

        if self._options['show_sql']:
            return f': {color_text(str(query), "darkgrey")}'
        else:
            return '.'
