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
from networkx.algorithms import shortest_path
from peewee import fn
from playhouse.migrate import PostgresqlMigrator, migrate

from sdssdb.connection import PeeweeDatabaseConnection
from sdssdb.utils import get_row_count
from sdsstools.color_print import color_text

import target_selection
from target_selection.exceptions import (TargetSelectionNotImplemented,
                                         XMatchError)
from target_selection.utils import (Timer, get_epoch, set_config_parameter,
                                    sql_apply_pm)


EPOCH = 2015.5


class Catalog(peewee.Model):
    """Model for the output table."""

    catalogid = peewee.BigIntegerField(null=False)
    iauname = peewee.TextField(null=True)
    ra = peewee.DoubleField(null=False)
    dec = peewee.DoubleField(null=False)
    pmra = peewee.FloatField(null=True)
    pmdec = peewee.FloatField(null=True)
    parallax = peewee.FloatField(null=True)
    lead = peewee.TextField(null=False)
    version = peewee.TextField(null=False)

    class Meta:
        primary_key = peewee.CompositeKey('catalogid', 'version')


def get_relational_model(model, prefix='catalog_to_'):

    meta = model._meta
    pk = meta.primary_key

    if isinstance(pk, str) and pk == '__composite_key__':
        raise XMatchError(f'composite pk found for model {model.__name__!r}.')

    # Auto/Serial are automatically PKs. Convert them to integers
    # to avoid having two pks in the relational table.
    if pk.__class__.field_type == 'AUTO':
        model_pk_class = peewee.IntegerField
    elif pk.__class__.field_type == 'BIGAUTO':
        model_pk_class = peewee.BigIntegerField
    else:
        model_pk_class = pk.__class__

    class BaseModel(peewee.Model):

        id = peewee.BigAutoField(primary_key=True)
        catalogid = peewee.BigIntegerField()
        target_id = model_pk_class()
        distance = peewee.DoubleField(null=True)
        best = peewee.BooleanField(null=False, default=True)

        class Meta:
            database = meta.database
            schema = meta.schema

    model_prefix = ''.join(x.capitalize() or '_'
                           for x in prefix.rstrip().split('_'))

    RelationalModel = type(model_prefix + model.__name__, (BaseModel,), {})
    RelationalModel._meta.table_name = prefix + meta.table_name

    return RelationalModel


def XMatchModel(Model, resolution=None, ra_column=None, dec_column=None,
                pmra_column=None, pmdec_column=None, is_pmra_cos=True,
                parallax_column=None, epoch_column=None, epoch=None,
                epoch_format='jyear', has_duplicates=False,
                skip_unlinked=False, skip=False):
    """Expands the model `peewee:Metadata` with cross-matching parameters.

    The parameters defined can be accessed with the same name as
    ``Model._meta.xmatch.<parameter>`` (e.g.,
    ``Model._meta.xmatch.has_duplicates``).

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
    :obj:`peewee:Model`
        The same input model with the additional cross-matching parameters
        added to the metadata namespace ``xmatch``.

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
    """Prepares and runs catalogue cross-matching.

    This class prepares the execution of a cross-matching between multiple
    catalogues, the result being an output table of unique targets linked to
    each input catalog via a relational table. Target coordinates are
    propagated to a common epoch when proper motions are available.
    Instantiating the class only prepares the process and sets up the
    processing order; cross-matching itself happens by calling the `.run`
    method.

    The output table contains a sequential integer identifier for each unique
    target (``catalogid``), along with the following columns: ``ra``, ``dec``,
    ``pmra``, ``pmdec``, ``parallax``, ``lead``, and ``version``. ``version``
    matches the version string for the cross-matching process. ``lead``
    indicates from which one of the input catalogues the coordinates were
    obtained.

    The output table is related to each input catalogue via a many-to-many
    table with the format ``<output_table>_to_<catalogue_table>`` where
    ``<output_table>`` is the table name of the output table and
    ``<catalog_table>`` is the table name of the referred catalogue. Each
    relational table contains a unique identifier column, ``id``,
    ``catalogid``, ``target_id`` pointing to the primary key of the referred
    catalogue, and a boolean ``best`` indicating whether it is the best
    (closest) possible match when a target in the output table correspond to
    multiple targets in the input catalogue.

    The cross-matching process roughly follows the following process:

    |

    .. image:: _static/Catalogdb_Crossmatch.png
        :scale: 75 %
        :align: center

    |

    In practice the cross-matching process begins by creating a graph of all
    nodes (tables to be processed and additional tables, ``extra_nodes``, in
    the database) and edges (foreign keys relating two tables). This graph is
    used to determine join conditions and to establish the order in which the
    input models will be processed. The processing order is determined by the
    ``order`` and ``key`` input parameters. When ``key='row_count'``, tables
    are sorted by number of decreasing rows so that tables with more targets
    are processed first (note that to speed things up the row count is always
    the latest known approximate determined by ``ANALYZE``); if
    ``key='resolution'`` the associated spatial resolution for a catalogue is
    used to process catalogues with high resolution first. If
    ``order='hierarchical'``, all the tables are divided into as many
    disconnected subgraphs as exist; then for each subgraph the maximum row
    count or minim resolution is calculated (depending on the value of
    ``key``). Subgraphs are sorted based on this result and then tables
    belonging to each subgraph are sorted by key. If ``order='global'`` the
    ``key`` ordering is applied to all tables without taking into account
    subgraphs.

    Once the order has been determined and when `.run` is called, each table
    model is processed in order. The first model is just ingested completely
    into the output table and its associated relational table is created.

    For each additional model the following three stages are applied:

    - In *phase 1* we determine what targets in the input model have an
      existing cross-match to targets already ingested into the output table.
      To do that we build all possible joins between the model and the output
      table. If multiple joins are possible via a given table only the shortest
      is used (see `.get_simple_paths`). For all matched targets we insert
      entries in the relational table. The *lead* of the original entries is
      not changed.

    - In *phase 2* we perform the actual cross-match between targets in the
      output table and the ones in the input catalogue. Currently the only
      cross-matching method available is a spatial cone query with radius
      ``query_radius``. All matched targets are added to the relational table
      and the one with the smallest distance is defined as *best*.

    - In *phase 3* we determine any target in the input catalogue that has not
      been cross-matched at this point and insert them into the output table as
      new entries. The *lead* is set to the input catalogue and the one-to-one
      match is added to the relational table.

    After each model has been processed the relational table is modified to
    add foreign keys and indexes to speed up future queries. The model graph
    is updated to reflect the newly created relational table and foreign keys.

    In addition to the limitations of the spatial cone query method, the
    following caveats are known:

    - Input tables with duplicate targets are not currently supported.

    Parameters
    ----------
    database : peewee:PostgresqlDatabase
        A `peewee:PostgresqlDatabase` to the database the tables to
        cross-match.
    models : list
        The list of `.XMatchModel` classes to be cross-matched. If the model
        correspond to a non-existing table it will be silently ignored.
    version : str
        The cross-matching version.
    extra_nodes : list
        List of PeeWee models to be used as extra nodes for joins (i.e.,
        already established cross-matches between catalogues). This models
        are not processed or inserted into the output table.
    order : str or list
        The type of sorting to be applies to the input models to decide in
        what order to process them. Currently allowed values are
        ``'hierarchical'`` and ``'global'`` (refer to the description above).
        The order can also be a list of table names, in which case that order
        is used without any further sorting.
    key : str
        The key to be used while sorting. Can be ``'row_count'`` or
        ``'resolution'``.
    epoch : float
        The epoch to which to convert all the target coordinates as they are
        inserted into the output table.
    start_node : str
        If specified, the name of the table that will be inserted first
        regarding of the above sorting process.
    query_radius : float
        The radius, in degrees, for cross-matching between existing targets.
        Used in phase 2. Defaults to 1 arcsec.
    schema : str
        The schema in which all the tables to cross-match live (multiple
        schemas are not supported), and the schema in which the output table
        will be created.
    output_table : str
        The name of the output table. Defaults to ``catalog``.
    allow_existing : bool
        By default instantiation will fail if the output table already exists
        and contains entries associated with this cross-matching version. This
        option allows to continue but should be used only for testing.
    log_path : str
        The path to which to log or `False` to disable file logging.
    debug : bool or int
        Controls the level to which to log to the screen. If `False` no logging
        is done to ``stdout``. If `True` the logging level is set to ``debug``
        (all messages). It's also possible specify a numerical value for the
        `logging level <logging>`.
    show_sql : bool
        Whether to log the full SQL queries being run.
    sample_region : tuple
        Allows to specify a 3-element tuple with the
        (``ra``, ``dec``, ``radius``) of the region to which to limit the
        cross-match. All values must be in degrees.
    disable_seqscan : bool
        Whether to disable sequential scans completely. Depending on the tables
        to process and the database server configuration, this may sometimes
        be needed to force the query planner to use the Q3C indexes.

    """

    def __init__(self, database, models, version, extra_nodes=[],
                 order='hierarchical', key='row_count', epoch=EPOCH,
                 start_node=None, query_radius=None, schema='catalogdb',
                 output_table='catalog', allow_existing=False,
                 log_path='./xmatch_{version}.log', debug=False,
                 show_sql=False, sample_region=None, disable_seqscan=False):

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

        self._options = {'query_radius': query_radius or 1 / 3600.,
                         'allow_existing': allow_existing,
                         'show_sql': show_sql,
                         'sample_region': sample_region,
                         'epoch': epoch,
                         'disable_seqscan': disable_seqscan}

        self._check_models()

        self.model_graph = None
        self.update_model_graph()

        self.process_order = []
        self.set_process_order(order=order, key=key, start_node=start_node)

    @classmethod
    def read(cls, in_models, version, config_file=None, **kwargs):
        """Instantiates `.XMatchPlanner` from a configuration file.

        The YAML configuration file must organised by version string (multiple
        versions can live in the same file). Any parameter that
        `.XMatchPlanner` accepts can be passed via the configuration file.
        Additionally, the configuration file accepts the two extra parameters:
        ``exclude``, a list of table names that will be ignored (this is useful
        if you pass a datbase or base class as ``in_models`` and want to ignore
        some of the models), and ``tables``, a dictionary of table names with
        parameters to be passed to `.XMatchModel` for its corresponding model.
        An example of a valid configuration files is:

        .. code-block:: yaml

            '0.1.0':
                order: hierarchical
                key: resolution
                query_radius: 0.000277778  # 1 arcsec
                schema: catalogdb
                output_table: catalog
                start_node: tic_v8
                debug: true
                log_path: false
                exclude: ['catwise']
                tables:
                    tic_v8:
                        ra_column: ra
                        dec_column: dec
                        pmra_column: pmra
                        pmdec_column: pmdec
                        is_pmra_cos: true
                        parallax_column: plx
                        epoch: 2015.5
                    gaia_dr2_source:
                        ra_column: ra
                        dec_column: dec
                        pmra_column: pmra
                        pmdec_column: pmdec
                        is_pmra_cos: true
                        parallax_column: parallax
                        epoch: 2015.5
                        skip: true

        Note that only models that match the table names in ``tables`` will be
        passed to `.XMatchPlanner` to be processed; any other table will be
        used as an extra joining node unless it's listed in ``exclude``, in
        which case it will be ignored completely. It's possible to set the
        ``skip`` option for a table; this has the same effect as removing the
        entry in ``table``.

        Parameters
        ----------
        in_models
            The models to cross-match. Can be a list or tuple of PeeWee
            :obj:`peewee:Model` instances, a base class from which all the
            models to use subclass, or a
            `~sdssdb.connection.PeeweeDatabaseConnection` to the database
            containing the models. In the latter case, the models must have
            been imported so that they are available via the ``models``
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
        if isinstance(in_models, (list, tuple)):
            models = in_models
        elif inspect.isclass(in_models) and issubclass(in_models, peewee.Model):
            database = in_models._meta.database
            models = set(in_models.__subclasses__())
            while True:
                old_models = models.copy()
                for model in old_models:
                    models |= set(model.__subclasses__())
                if models == old_models:
                    break
        elif isinstance(in_models, PeeweeDatabaseConnection):
            database = in_models
            models = database.models.values()
        else:
            raise TypeError(f'invalid input of type {type(in_models)!r}')

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
                      model._meta.schema == schema and model not in exclude)}

        xmatch_models = {}
        for table_name in table_config:
            if table_name not in models:
                continue
            table_params = table_config[table_name] or {}
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
                self.extra_nodes[Catalog._meta.table_name] = Catalog
            elif table_name.startswith(Catalog._meta.table_name + '_to_'):
                self.extra_nodes[table_name] = model
            elif meta.xmatch.skip:
                self.extra_nodes[table_name] = model
            else:
                continue

            remove_models.append(model)

        [self.models.pop(model._meta.table_name) for model in remove_models]
        [self.extra_nodes.pop(y) for y in [x for x in self.extra_nodes
                                           if not self.extra_nodes[x].table_exists()]]

        if len(self.models) == 0:
            raise XMatchError('no models to cross-match.')

        error_msg = (f'{Catalog._meta.table_name!r} contains '
                     'records for this version of cross-matching '
                     f'({self.version}).')

        # Check if Catalog already has entries for this xmatch version.
        if Catalog.table_exists() and (Catalog.select()
                                       .where(Catalog.version == self.version)
                                       .limit(1)
                                       .count() > 0):
            if self._options['allow_existing']:
                self.log.warning(error_msg)
            else:
                raise XMatchError(error_msg + ' This problem needs to be '
                                  'solved manually.')

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

            # A bit of a hack here. There is no FK between catalog_to_table
            # and catalog because it depends on catalogid and version. So
            # We add an edge manually and we'll deal with the join later.
            if table_name.startswith(Catalog._meta.table_name + '_to_'):
                ref_table_name = Catalog._meta.table_name
                self.model_graph.add_edge(table_name, ref_table_name)

        return self.model_graph

    def set_process_order(self, order='hierarchical', key='row_count',
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

        graph = self.model_graph.copy()

        for model in self.extra_nodes.values():
            graph.remove_node(model._meta.table_name)

        if order == 'hierarchical':
            subgraphs = networkx.connected_components(graph)
        else:
            subgraphs = [node for node in subgraphs.nodes]

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
                                           for tn in sg])
                    # Use -total_row_count to avoid needing reverse order.
                    subgraphs_ext.append((sg, -total_row_count, 1))
                elif key == 'resolution':
                    resolution = [self.models[tn]._meta.xmatch.resolution
                                  for tn in sg]
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
            ordered_tables = []
            for sgo in subgraphs_ordered:
                sg_ext = []
                for table_name in sgo:
                    model = self.models[table_name]
                    if start_node and start_node == table_name:
                        sg_ext.append((table_name, numpy.nan, 0))
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
        `list`
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
                spath = shortest_path(graph, source, dest)
                paths.append(spath)
                for node in spath[1:-1]:
                    graph.remove_node(node)
            except networkx.NetworkXNoPath:
                break

        if len(paths) == 0:
            return []

        if return_models:
            nodes = self.model_graph.nodes
            paths = [[nodes[node]['model'] for node in path] for path in paths]

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
            raise TargetSelectionNotImplemented(
                'handling of tables with duplicates is not implemented.')

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
            self.log.debug(f'Adding indexes and FKs to '
                           f'{rel_model._meta.table_name!r}.')
            self.add_fks(rel_model, model)

        self.log.header = ''

    def _get_model_fields(self, model):
        """Returns the model fields needed to populate Catalog."""

        meta = model._meta
        xmatch = meta.xmatch
        fields = meta.fields

        # List of fields that will become the SELECT clause.
        model_fields = []

        ra_field = fields[xmatch.ra_column]
        dec_field = fields[xmatch.dec_column]

        to_epoch = self._options['epoch']

        # RA, Dec, and pro[er motion fields.
        if xmatch.pmra_column:

            pmra_field = fields[xmatch.pmra_column]
            pmdec_field = fields[xmatch.pmdec_column]

            if (xmatch.epoch and xmatch.epoch != to_epoch) or xmatch.epoch_column:
                delta_years = to_epoch - get_epoch(model)
                ra_field, dec_field = sql_apply_pm(ra_field, dec_field,
                                                   pmra_field, pmdec_field,
                                                   delta_years,
                                                   xmatch.is_pmra_cos)
                if not xmatch.is_pmra_cos:
                    pmra_field *= fn.cos(fn.radians(dec_field))

        else:

            pmra_field = peewee.Value(None).cast('REAL')
            pmdec_field = peewee.Value(None).cast('REAL')

        # Actually add the RA/Dec fields as defined above
        model_fields.extend([ra_field, dec_field])
        model_fields.extend([pmra_field, pmdec_field])

        # Parallax
        if xmatch.parallax_column:
            model_fields.append(fields[xmatch.parallax_column])
        else:
            model_fields.append(peewee.Value(None).cast('REAL'))

        return model_fields

    def _load_catalog(self, model, rel_model, query=None):
        """Ingests a query into the output and relational tables."""

        meta = model._meta
        table_name = meta.table_name

        model_fields = self._get_model_fields(model)

        pk = model._meta.primary_key
        ra_field = meta.fields[meta.xmatch.ra_column]
        dec_field = meta.fields[meta.xmatch.dec_column]

        disable_seqscan = self._options['disable_seqscan']

        # Get the max catalogid currently in the table for the version.
        max_cid = (Catalog
                   .select(fn.MAX(Catalog.catalogid))
                   .where(Catalog.version == self.version)
                   .scalar() or 0)

        if query is None:
            query = model.select()

        # Build the query and exclude rows without RA/Dec. Note that this
        # resets any SELECT field passed to the query but keeps the WHEREs.
        query = (query
                 .select(fn.row_number().over() + max_cid,
                         model._meta.primary_key,
                         peewee.SQL('true'))
                 .where(ra_field.is_null(False),
                        dec_field.is_null(False))
                 .order_by(meta.primary_key.asc()))

        # Apply any restriction to the sampled area.
        if self._options['sample_region']:
            sample_region = self._options['sample_region']
            query = query.where(fn.q3c_radial_query(ra_field,
                                                    dec_field,
                                                    sample_region[0],
                                                    sample_region[1],
                                                    sample_region[2]))

        with Timer() as timer:

            with self.database.atomic():

                with set_config_parameter(self.database, 'enable_seqscan',
                                          'OFF' if disable_seqscan else 'ON',
                                          log=self.log):

                    # CTE to populate the relational table.
                    rel_insert = peewee.CTE('rel_insert',
                                            rel_model.insert_from(
                                                query,
                                                [rel_model.catalogid,
                                                 rel_model.target_id,
                                                 rel_model.best])
                                            .returning(rel_model.catalogid,
                                                       rel_model.target_id))

                    # INSERT INTO catalog using the catalogids and targetids
                    # returned by y. Return only the number of rows inserted.
                    insert_query = Catalog.insert_from(
                        rel_insert.select(rel_insert.c.catalogid,
                                          *model_fields,
                                          peewee.Value(table_name),
                                          peewee.Value(self.version))
                        .join(model, on=(pk == rel_insert.c.target_id)),
                        [Catalog.catalogid,
                         Catalog.ra,
                         Catalog.dec,
                         Catalog.pmra,
                         Catalog.pmdec,
                         Catalog.parallax,
                         Catalog.lead,
                         Catalog.version]
                    ).returning().with_cte(rel_insert)

                    self.log.debug(f'Running INSERT query'
                                   f'{self._get_sql(insert_query)}')

                    n_rows = insert_query.execute()

        self.log.debug(f'Inserted {n_rows} rows in {timer.interval:.3f} s.')

        return

    def _create_relational_model(self, model):
        """Creates a relational table for a given model.

        Returns the relational model and `True` if the table was created or
        `False` if it already existed.

        """

        cat_table = Catalog._meta.table_name

        RelationalModel = get_relational_model(model, prefix=cat_table + '_to_')
        rtname = RelationalModel._meta.table_name

        if not RelationalModel.table_exists():
            RelationalModel.create_table()
            created = True
        else:
            if RelationalModel.select().count() > 0:
                self.log.warning(f'Relational table {rtname!r} is not empty!')

            created = False

        # Add foreign key field here. We want to avoid Peewee creating it
        # as a constraint and index if the table is created because that would
        # slow down inserts. We'll created them manually with add_fks.
        # Note that we do not create an FK between the relational model and
        # Catalog because the relationship is only unique on (catalogid, version).
        RelationalModel._meta.add_field(
            'target',
            peewee.ForeignKeyField(model, column_name='target_id', backref='+'))

        self.extra_nodes[rtname] = RelationalModel
        self.update_model_graph(silent=True)

        return RelationalModel, created

    def _build_join(self, path):
        """Returns a build query for a given join path."""

        model = path[0]
        meta = model._meta
        xmatch = meta.xmatch

        query = model.select()
        for inode in range(1, len(path)):
            if path[inode] is Catalog:
                query = query.join(Catalog,
                                   on=((Catalog.catalogid == path[inode - 1].catalogid) &
                                       (Catalog.version == self.version)))
            else:
                query = query.join(path[inode])

        if self._options['sample_region']:
            sample_region = self._options['sample_region']
            query = query.where(fn.q3c_radial_query(meta.fields[xmatch.ra_column],
                                                    meta.fields[xmatch.dec_column],
                                                    sample_region[0],
                                                    sample_region[1],
                                                    sample_region[2]))

        return query

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

                        query = (self._build_join(join_models)
                                 .select(model._meta.primary_key,
                                         Catalog.catalogid,
                                         peewee.Value(True))
                                 .where(Catalog.version == self.version))

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
        disable_seqscan = self._options['disable_seqscan']

        # TODO: find a way to get this delta that is not too time consuming.
        # self.log.debug('Determining maximum delta epoch between tables.')
        # max_delta_epoch = peewee.Value(
        #     model.select(fn.MAX(fn.ABS(model_epoch - catalog_epoch))).scalar())
        max_delta_epoch = 50.

        # Create a subquery that returns the entries in the model that are
        # within query_radius of the Catalog target.
        subq = (model.select(model_pk.alias('target_id'),
                             fn.q3c_dist_pm(Catalog.ra,
                                            Catalog.dec,
                                            Catalog.pmra,
                                            Catalog.pmdec,
                                            1,  # Catalog pmra is pmra*cos by definition.
                                            catalog_epoch,
                                            model_ra,
                                            model_dec,
                                            model_epoch).alias('distance'))
                .where(fn.q3c_join_pm(Catalog.ra,
                                      Catalog.dec,
                                      Catalog.pmra,
                                      Catalog.pmdec,
                                      1,
                                      catalog_epoch,
                                      model_ra,
                                      model_dec,
                                      model_epoch,
                                      max_delta_epoch,
                                      self._options['query_radius'])))

        # We'll partition over each group of targets that are within
        # radius of a catalogid and mark the one with the smallest
        # distance to the Catalog target as best.
        partition = (fn.first_value(subq.c.target_id)
                     .over(partition_by=[Catalog.catalogid],
                           order_by=[subq.c.distance.asc()]))
        best = peewee.Value(partition == subq.c.target_id)

        # Determine what targets in the model haven't been matched in phase 1.
        # We do this starting from the Catalog side and determining what targets
        # don't have an entry in the relational table.
        unmatched = (Catalog.select(Catalog.catalogid,
                                    subq.c.target_id,
                                    subq.c.distance,
                                    best.alias('best'))
                     .join(rel_model, peewee.JOIN.LEFT_OUTER,
                           on=((Catalog.catalogid == rel_model.catalogid) &
                               (Catalog.version == self.version)))
                     .where(rel_model.catalogid.is_null(),
                            Catalog.version == self.version))

        # Limit the query region
        if self._options['sample_region']:
            racen, deccen, radius = self._options['sample_region']
            unmatched = unmatched.where(fn.q3c_radial_query(Catalog.ra,
                                                            Catalog.dec,
                                                            racen,
                                                            deccen,
                                                            radius))

        # Do a lateral join (for each loop). Remove pairs in which
        # target_id is null (i.e., where there is no target within the query
        # radius).
        xmatches = (unmatched
                    .join(subq, 'LEFT JOIN LATERAL', on=peewee.SQL('true'))
                    .order_by(Catalog.catalogid)
                    .where(subq.c.target_id.is_null(False)))

        with Timer() as timer:

            with self.database.atomic():

                with set_config_parameter(self.database, 'enable_seqscan',
                                          'OFF' if disable_seqscan else 'ON',
                                          log=self.log):

                    self.log.debug('Inserting cross-matched targets into '
                                   f'{rel_table_name!r}{self._get_sql(xmatches)}')

                    # Insert as CTE to gather some statistincs
                    insert_cte = peewee.CTE(
                        'insert_cte',
                        rel_model.insert_from(xmatches,
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

                    n_catalogid, n_target_id = list(insert_query
                                                    .tuples()
                                                    .execute(self.database))[0]

        ratio_msg = (f' (ratio {n_target_id / n_catalogid:.2f}). '
                     if n_catalogid > 0 else '. ')
        self.log.debug(f'Cross-matched {n_catalogid} catalogids '
                       f'associated with {n_target_id} targets in {table_name!r}' +
                       ratio_msg + f'Run in {timer.interval:.3f} s.')

    def _run_phase_3(self, model, rel_model):
        """Add non-matched targets to Catalog and the relational table."""

        self.log.info('Phase 3: adding non cross-matched targets.')

        # Limit the query to targets not matched.
        query = (model
                 .select()
                 .where(~fn.EXISTS(rel_model.select(rel_model.target_id)
                                   .join(Catalog,
                                         on=((Catalog.catalogid == rel_model.catalogid) &
                                             (Catalog.version == self.version)))
                                   .where(rel_model.target_id == model._meta.primary_key))))

        # Insert the records
        self._load_catalog(model, rel_model, query=query)

    def _get_sql(self, query):
        """Returns coulourised SQL text for logging."""

        if self._options['show_sql']:
            return f': {color_text(str(query), "darkgrey")}'
        else:
            return '.'

    def add_fks(self, relational_model, to_model):
        """Adds foreign keys and indexes to a relational model class."""

        rtschema = relational_model._meta.schema
        rtname = relational_model._meta.table_name

        to_schema = to_model._meta.schema
        to_table_name = to_model._meta.table_name

        migrator = PostgresqlMigrator(self.database)

        with self.database.atomic():

            migrate(
                migrator.set_search_path(rtschema),
                migrator.add_index(rtname, ('catalogid',), False),
                migrator.add_index(rtname, ('target_id',), False),
                migrator.add_constraint(
                    rtname, 'target_id',
                    peewee.SQL(
                        f'FOREIGN KEY (target_id) '
                        f'REFERENCES {to_schema}.{to_table_name} '
                        f'({to_model._meta.primary_key.column_name})')))

        migrator.database.close()
