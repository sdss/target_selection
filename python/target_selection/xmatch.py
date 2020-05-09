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
import types
import warnings

import networkx
import numpy
import peewee
import yaml
from networkx.algorithms import shortest_path
from peewee import SQL, Select, fn

from sdssdb.connection import PeeweeDatabaseConnection
from sdssdb.utils import get_row_count
from sdsstools import merge_config
from sdsstools.color_print import color_text

import target_selection
from target_selection.exceptions import (TargetSelectionNotImplemented,
                                         TargetSelectionUserWarning,
                                         XMatchError)
from target_selection.utils import (Timer, get_configuration_values, get_epoch,
                                    sql_apply_pm, vacuum_outputs, vacuum_table)


EPOCH = 2015.5
QUERY_RADIUS = 1.


class Version(peewee.Model):
    """Model for the version table."""

    id = peewee.AutoField()
    version = peewee.TextField()
    tag = peewee.TextField()

    class Meta:
        table_name = 'version'


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
    version_id = peewee.IntegerField(null=False, index=True)


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

        id = peewee.BigAutoField()
        catalogid = peewee.BigIntegerField(null=False, index=True)
        target_id = model_pk_class(null=False)
        version_id = peewee.SmallIntegerField(null=False)
        distance = peewee.DoubleField(null=True)
        best = peewee.BooleanField(null=False)

        class Meta:
            database = meta.database
            schema = meta.schema
            indexes = [(('version_id', 'target_id'), False)]

    model_prefix = ''.join(x.capitalize() or '_'
                           for x in prefix.rstrip().split('_'))

    RelationalModel = type(model_prefix + model.__name__, (BaseModel,), {})
    RelationalModel._meta.table_name = prefix + meta.table_name

    return RelationalModel


def XMatchModel(Model, resolution=None, ra_column=None, dec_column=None,
                pmra_column=None, pmdec_column=None, is_pmra_cos=True,
                parallax_column=None, epoch_column=None, epoch=None,
                epoch_format='jyear', has_duplicates=False,
                has_missing_coordinates=False, skip=False,
                skip_phases=None, query_radius=None, join_weight=1):
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
    has_missing_coordinates : bool
        Whether the catalogue contains rows in which the RA/Dec are null.
    skip : bool
        If `True`, the table will be used as a join node but will not be
        cross-matched. This is useful for testing and also if the table is
        in a previous version of the configuration file and you are using the
        ``base_version`` option but want to remove that table. It can also
        be used when setting a ``join_path`` for a model but otherwise don't
        want the table to be processed.
    skip_phases : list
        A list of cross-matching phases to be skipped for this model. Refer
        to the `.XMatchPlanner` documentation for definitions on what each
        phase does.
    query_radius : float
        The radius, in arcsec, to use in the radial query for cross-matching.
        If not provided defaults to the `.XMatchPlanner` value.
    join_weight : float
        The weight used by `.XMatchPlanner.get_join_paths` to determine
        the cost of using this table as a join. Lower weights translate to
        better chances of that join path to be selected.

    Returns
    -------
    :obj:`peewee:Model`
        The same input model with the additional cross-matching parameters
        added to the metadata namespace ``xmatch``.

    """

    meta = Model._meta
    meta.xmatch = types.SimpleNamespace()

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
    meta.xmatch.has_missing_coordinates = has_missing_coordinates
    meta.xmatch.skip = skip
    meta.xmatch.skip_phases = skip_phases or []
    meta.xmatch.query_radius = query_radius

    meta.xmatch.row_count = int(get_row_count(meta.database,
                                              meta.table_name,
                                              schema=meta.schema,
                                              approximate=True))

    meta.xmatch.join_weight = join_weight

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
    ``pmra``, ``pmdec``, ``parallax``, ``lead``, and ``version_id``.
    ``version_id`` relates the record with the ``version`` table which contains
    the cross-matching version and the version of the code used when it was
    run. ``lead`` indicates from which one of the input catalogues the
    coordinates were obtained.

    The output table is related to each input catalogue via a many-to-many
    table with the format ``<output_table>_to_<catalogue_table>`` where
    ``<output_table>`` is the table name of the output table and
    ``<catalog_table>`` is the table name of the referred catalogue. Each
    relational table contains the unique identifier column, ``catalogid``,
    ``target_id`` pointing to the primary key of the referred catalogue,
    ``version_id``, and a boolean ``best`` indicating whether it is the best
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
    into the output table and its associated relational table is created
    (the equivalent of *phase 3* below).

    For each additional model the following three stages are applied:

    - In *phase 1* we determine what targets in the input model have an
      existing cross-match to targets already ingested into the output table.
      To do that we build all possible joins between the model and the output
      table. If multiple joins are possible via a given table only the shortest
      is used (see `.get_join_paths`). For all matched targets we insert
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

    - In *phase 2* there is no current measure in place for the same target
      to be associated with more than one catalogid (i.e., each cross-match
      is performed independently).

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
        The radius, in arcsec, for cross-matching between existing targets.
        Used in phase 2. Defaults to 1 arcsec.
    schema : str
        The schema in which all the tables to cross-match live (multiple
        schemas are not supported), and the schema in which the output table
        will be created.
    output_table : str
        The name of the output table. Defaults to ``catalog``.
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
        cross-match. All values must be in degrees. It can also be a list
        of tuples, in which case the union of all the regions will be sampled.

    """

    def __init__(self, database, models, version, extra_nodes=[],
                 order='hierarchical', key='row_count', epoch=EPOCH,
                 start_node=None, query_radius=None, schema='catalogdb',
                 output_table='catalog', log_path='./xmatch_{version}.log',
                 debug=False, show_sql=False, sample_region=None):

        self.log = target_selection.log
        self.log.header = ''

        if log_path:
            log_path = os.path.realpath(log_path)
            if self.log.fh:
                self.log.removeHandler(self.log.fh)
                self.log.fh = None
            self.log.start_file_logger(log_path.format(version=version),
                                       rotating=False, mode='a')

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
        self.log.info(f'version = {self.version!r}; '
                      f'tag = {target_selection.__version__!r}.')

        # Sets the metadata of the Catalog and Version tables.
        Catalog._meta.schema = schema
        Catalog._meta.table_name = output_table
        Version._meta.schema = schema
        self.database.bind([Catalog, Version])

        self.models = {model._meta.table_name: model for model in models}
        self.extra_nodes = {model._meta.table_name: model for model in extra_nodes}

        self._options = {'query_radius': query_radius or QUERY_RADIUS,
                         'show_sql': show_sql,
                         'sample_region': sample_region,
                         'epoch': epoch}

        if self._options['sample_region']:
            sample_region = self._options['sample_region']
            self.log.warning(f'Using sample region {sample_region!r}.')

        self._check_models()

        self._log_db_configuration()

        self.model_graph = None
        self.update_model_graph()

        self.process_order = []
        self.set_process_order(order=order, key=key, start_node=start_node)

        self._version_id = None
        self._max_cid = None
        self._catalog_count = 0

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
                query_radius: 1.
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

        It is also possible to use a parameter ``base_version`` pointing to a
        previous version string. In that case the previous version
        configuration will be used as base and the new values will be merged
        (the update happens recursively as with normal Python dictionaries).

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

        config = XMatchPlanner._read_config(config_file, version)

        table_config = config.pop('tables', {}) or {}
        exclude = config.pop('exclude', []) or []

        assert 'schema' in config, 'schema is required in configuration.'
        schema = config['schema']

        models = {model._meta.table_name: model for model in models
                  if (model.table_exists() and
                      model._meta.schema == schema and
                      model._meta.table_name not in exclude)}

        xmatch_models = {}
        for table_name in table_config:
            if table_name not in models:
                continue
            table_params = table_config[table_name] or {}
            xmatch_models[table_name] = XMatchModel(models[table_name], **table_params)

        extra_nodes = [models[table_name] for table_name in models
                       if table_name not in xmatch_models]

        config.update(kwargs)

        signature = inspect.signature(XMatchPlanner)

        valid_kw = {}
        for kw in config:
            if kw not in signature.parameters:
                warnings.warn(f'ignoring invalid configuration value {kw!r}.',
                              TargetSelectionUserWarning)
                continue
            valid_kw[kw] = config[kw]

        return cls(database, xmatch_models.values(), version,
                   extra_nodes=extra_nodes, **valid_kw)

    @staticmethod
    def _read_config(file_, version):
        """Reads the configuration file, recursively."""

        config = yaml.load(open(file_, 'r'), Loader=yaml.SafeLoader)

        assert version in config, f'version {version!r} not found in configuration.'

        base_version = config[version].pop('base_version', None)
        if base_version:
            config = merge_config(config[version],
                                  XMatchPlanner._read_config(file_, base_version))
        else:

            config = config[version]

        return config

    def _check_models(self):
        """Checks the input models."""

        catalog_tname = Catalog._meta.table_name

        remove_models = []

        for model in self.models.values():

            meta = model._meta
            table_name = meta.table_name

            if not model.table_exists():
                self.log.warning(f'table {table_name!r} does not exist.')
            elif table_name == catalog_tname:
                self.extra_nodes[catalog_tname] = Catalog
            elif table_name.startswith(catalog_tname + '_to_'):
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

        # Check if Catalog already has entries for this xmatch version.
        if Catalog.table_exists():
            version_id = (Version.select()
                          .where(Version.version == self.version)
                          .scalar())
            exists = (peewee.Select(
                columns=[fn.EXISTS(Catalog
                                   .select(SQL('1'))
                                   .where(Catalog.version_id == version_id))])
                      .tuples()
                      .execute(self.database))[0][0]
            if version_id and exists:
                self.log.warning(f'{catalog_tname!r} contains '
                                 'records for this version of cross-matching '
                                 f'({self.version}).')

    def _log_db_configuration(self):
        """Logs some key database configuration parameters."""

        parameters = ['shared_buffers',
                      'effective_cache_size',
                      'wal_buffers',
                      'effective_io_concurrency',
                      'work_mem',
                      'max_worker_processes',
                      'random_page_cost',
                      'seq_page_cost',
                      'cpu_index_tuple_cost',
                      'cpu_operator_cost',
                      'default_statistics_target',
                      'temp_buffers',
                      'plan_cache_mode',
                      'geqo_effort',
                      'force_parallel_mode',
                      'enable_seqscan',
                      'enable_nestloop']

        values = get_configuration_values(self.database, parameters)

        self.log.debug('Current database configuration parameters.')
        for parameter in values:
            self.log.debug(f'{parameter} = {values[parameter]}')

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

                # Determines the join weight as the average of the
                # weights of the joined nodes.
                if hasattr(model._meta, 'xmatch'):
                    model_weight = model._meta.xmatch.join_weight
                else:
                    model_weight = 1

                if hasattr(fk_model._meta, 'xmatch'):
                    fk_model_weight = fk_model._meta.xmatch.join_weight
                else:
                    fk_model_weight = 1

                join_weight = 0.5 * (model_weight + fk_model_weight)

                self.model_graph.add_edge(table_name, ref_table_name,
                                          join_weight=join_weight)

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

        if isinstance(order, (list, tuple)):
            self.log.info(f'processing order: {order}')
            self.process_order = order
            return order

        assert order in ['hierarchical', 'global'], f'invalid order value {order!r}.'
        self.log.info(f'processing order mode is {order!r}')

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

    def get_join_paths(self, source, return_models=False):
        """Determines all possible join path between two tables.

        Instead of returning all possible simple paths between the ``source``
        and output table, this algorithm follows a "scorched earth" approach
        in which once an edge has been used for a join it cannot be used again.
        This produces only distint joins between the two nodes. Paths that
        include only the source and destination through their relational table
        are ignored, as are paths in which the last node before the output
        table has not yet been processed.

        Weights can be defined by setting the ``join_weight`` value in
        `.XMatchModel`. The weight for each edge is the average of the the
        ``join_weight`` of the two nodes joined. The weight defaults to 1.
        Lower weights translate to better chances of that join path to be
        selected.

        Parameters
        ----------
        source : str
            The initial table for the path.
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

        dest = Catalog._meta.table_name
        porder = self.process_order

        rel_table_name = Catalog._meta.table_name + '_to_' + source
        if rel_table_name in graph:
            graph.remove_node(rel_table_name)

        paths = []

        while True:
            try:
                spath = shortest_path(graph, source, dest, weight='join_weight')
                if (spath[-3] in porder and
                        porder.index(spath[-3]) < porder.index(source)):
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

    def run(self, vacuum=False, analyze=False, from_=None):
        """Runs the cross-matching process.

        Parameters
        ----------
        vacuum : bool
            Vacuum all output tables before processing new catalogues.
        analyze : bool
            Analyze all output tables before processing new catalogues.
        from_ : str
            Table from which to start running the process. Useful in reruns
            to skip tables already processed.

        """

        # If from_, we don't know how many targets are in the output table
        # so we set self._catalog_count to None and _get_larger_table
        # will take care of it.
        if from_:
            self._catalog_count = None

        if vacuum or analyze:
            cmd = ' '.join(('VACUUM' if vacuum else '',
                            'ANALYZE' if analyze else '')).strip()
            self.log.info(f'Running {cmd} on output tables.')
            vacuum_outputs(self.database, vacuum=vacuum, analyze=analyze,
                           schema=Catalog._meta.schema,
                           table=Catalog._meta.table_name)

        # Make sure the output table exists.
        if not Catalog.table_exists():
            # Add Q3C index for Catalog
            Catalog.index(SQL(f'CREATE INDEX ON '
                              f'{Catalog._meta.schema}.{Catalog._meta.table_name} '
                              f'(q3c_ang2ipix(ra, dec))'))
            self.database.create_tables([Catalog, Version])
            self.log.info(f'Created tables {Catalog._meta.table_name!r} and '
                          f'{Version._meta.table_name!r}.')

        version, vcreated = Version.get_or_create(version=self.version,
                                                  tag=target_selection.__version__)
        self._version_id = version.id
        if vcreated:
            vmsg = 'Added version record '
        else:
            vmsg = 'Using version record '
        self.log.info(vmsg + f'({self._version_id}, {self.version}, '
                             f'{target_selection.__version__}).')

        self.extra_nodes[Catalog._meta.table_name] = Catalog
        self.update_model_graph(silent=True)

        if vcreated or from_ is None:
            self._catalog_count = 0
        else:
            self._get_catalog_row_count()

        with Timer() as timer:
            for table_name in self.process_order:
                if (from_ and
                        self.process_order.index(table_name) < self.process_order.index(from_)):
                    self.log.warning(f'Skipping table {table_name}.')
                    continue
                model = self.models[table_name]
                self.process_model(model)

        self.log.info(f'Cross-matching completed in {timer.interval:.3f} s.')

    def process_model(self, model):
        """Processes a model, loading it into the output table."""

        table_name = model._meta.table_name
        self.log.header = f'[{table_name.upper()}] '

        self.log.info(f'Processing table {table_name}.')
        self._log_table_configuration(model)

        if model._meta.xmatch.has_duplicates:
            raise TargetSelectionNotImplemented(
                'handling of tables with duplicates is not implemented.')

        # Check if there are already records in catalog for this version.
        if self.process_order.index(model._meta.table_name) == 0:
            is_first_model = True
        else:
            is_first_model = False

        rel_model, rel_model_created = self._create_relational_model(model)
        rel_table_name = rel_model._meta.table_name
        if rel_model_created:
            self.log.debug(f'Created relational table {rel_table_name}.')
        else:
            self.log.debug(f'Relational table {rel_table_name} already exists.')

        if is_first_model:
            self._run_phase_3(model, rel_model)
        else:
            self._run_phase_1(model, rel_model)
            self._run_phase_2(model, rel_model)
            self._run_phase_3(model, rel_model)

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

            pmra_field = peewee.SQL('null')
            pmdec_field = peewee.SQL('null')

        # Actually add the RA/Dec fields as defined above
        model_fields.extend([ra_field, dec_field])
        model_fields.extend([pmra_field, pmdec_field])

        # Parallax
        if xmatch.parallax_column:
            model_fields.append(fields[xmatch.parallax_column])
        else:
            model_fields.append(peewee.SQL('null'))

        return model_fields

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

        query = model.select()
        for inode in range(1, len(path)):
            if path[inode] is Catalog:
                query = query.join(Catalog,
                                   on=(Catalog.catalogid == path[inode - 1].catalogid))
            else:
                query = query.join(path[inode])

        return query

    def _run_phase_1(self, model, rel_model):
        """Runs the linking against matched catalogids stage."""

        table_name = model._meta.table_name
        rel_table_name = rel_model._meta.table_name

        model_pk = model._meta.primary_key

        self.log.info('Phase 1: linking existing targets.')

        if 1 in model._meta.xmatch.skip_phases:
            self.log.warning('Skipping due to configuration.')
            return

        join_paths = self.get_join_paths(table_name)

        if len(join_paths) == 0:
            self.log.debug(f'No paths found between {table_name} and output table.')
            return False

        self.log.debug(f'Found {len(join_paths)} paths between '
                       f'{table_name} and output table.')

        for n_path, path in enumerate(join_paths):

            join_models = [self.model_graph.nodes[node]['model'] for node in path]

            with Timer() as timer:

                with self.database.atomic():

                    query = (self._build_join(join_models)
                             .select(model_pk,
                                     Catalog.catalogid,
                                     peewee.Value(self._version_id),
                                     peewee.SQL('true'))
                             .where(Catalog.version_id == self._version_id)
                             .where(~fn.EXISTS(rel_model
                                               .select(SQL('1'))
                                               .where(rel_model.version_id == self._version_id,
                                                      rel_model.target_id == model_pk))))

                    # In query we do not include a Q3C where for the sample region
                    # because Catalog for this version should already be sample
                    # region limited.

                    insert_query = rel_model.insert_from(
                        query, fields=[rel_model.target_id,
                                       rel_model.catalogid,
                                       rel_model.version_id,
                                       rel_model.best]).returning()

                    self.log.debug(f'Inserting linked targets into '
                                   f'{rel_table_name} with join path {path}'
                                   f'{self._get_sql(insert_query)}')

                    nids = insert_query.execute()

            self.log.debug(f'Linked {nids:,} records in {timer.interval:.3f} s.')

        self._analyze(rel_model)

    def _run_phase_2(self, model, rel_model):
        """Associates existing targets in Catalog with entries in the model."""

        meta = model._meta
        xmatch = meta.xmatch

        self.log.info('Phase 2: cross-matching against existing targets.')

        if 2 in xmatch.skip_phases:
            self.log.warning('Skipping due to configuration.')
            return

        table_name = meta.table_name
        rel_table_name = rel_model._meta.table_name

        model_pk = meta.primary_key
        model_ra = meta.fields[xmatch.ra_column]
        model_dec = meta.fields[xmatch.dec_column]

        catalog_epoch = self._options['epoch']

        query_radius = xmatch.query_radius or self._options['query_radius']

        # TODO: find a way to get this delta that is not too time consuming.
        # self.log.debug('Determining maximum delta epoch between tables.')
        # max_delta_epoch = peewee.Value(
        #     model.select(fn.MAX(fn.ABS(model_epoch - catalog_epoch))).scalar())
        max_delta_epoch = 50.

        # Should we use proper motions?
        model_epoch = get_epoch(model)
        use_pm = (isinstance(model_epoch, (peewee.Expression, peewee.Function)) or
                  model_epoch != catalog_epoch)

        # Determine which of the two tables is smaller. Q3C really wants the
        # larger table last.
        idx = self._get_larger_table(Catalog, model)

        if idx == 0:  # Catalog is larger

            self.log.debug(f'Cross-matching model against {Catalog._meta.table_name}.')

            if use_pm:

                model_pmra = meta.fields[xmatch.pmra_column]
                model_pmdec = meta.fields[xmatch.pmdec_column]
                model_is_pmra_cos = int(xmatch.is_pmra_cos)

                q3c_dist = fn.q3c_dist_pm(model_ra, model_dec,
                                          model_pmra, model_pmdec,
                                          model_is_pmra_cos, model_epoch,
                                          Catalog.ra, Catalog.dec, EPOCH)
                q3c_join = fn.q3c_join_pm(model_ra, model_dec,
                                          model_pmra, model_pmdec,
                                          model_is_pmra_cos, model_epoch,
                                          Catalog.ra, Catalog.dec,
                                          EPOCH, max_delta_epoch,
                                          query_radius / 3600.)
            else:

                q3c_dist = fn.q3c_dist(model_ra, model_dec, Catalog.ra, Catalog.dec)
                q3c_join = fn.q3c_join(model_ra, model_dec, Catalog.ra, Catalog.dec,
                                       query_radius / 3600.)

        else:  # Model is larger

            self.log.debug(f'Cross-matching {Catalog._meta.table_name} against model.')

            if use_pm:
                q3c_dist = fn.q3c_dist_pm(Catalog.ra, Catalog.dec,
                                          Catalog.pmra, Catalog.pmdec,
                                          1, catalog_epoch,
                                          model_ra, model_dec, model_epoch)
                q3c_join = fn.q3c_join_pm(Catalog.ra, Catalog.dec,
                                          Catalog.pmra, Catalog.pmdec,
                                          1, catalog_epoch,
                                          model_ra, model_dec,
                                          model_epoch, max_delta_epoch,
                                          query_radius / 3600.)
            else:
                q3c_dist = fn.q3c_dist(Catalog.ra, Catalog.dec, model_ra, model_dec)
                q3c_join = fn.q3c_join(Catalog.ra, Catalog.dec, model_ra, model_dec,
                                       query_radius / 3600.)

        # Get the cross-matched catalogid and model target pk (target_id),
        # and their distance. We need to use a materialised CTE here to
        # make sure the Q3C index is actually applied.
        xmatched = (Catalog
                    .select(Catalog.catalogid,
                            model_pk.alias('target_id'),
                            q3c_dist.alias('distance'),
                            Catalog.version_id)
                    .join(model, peewee.JOIN.CROSS)
                    .where(q3c_join))

        # This may break the use of the index but I think it's needed if
        # the model is the second table in q3c_join and has empty RA/Dec.
        if xmatch.has_missing_coordinates and use_pm:
            xmatched = xmatched.where(model_ra.is_null(False),
                                      model_dec.is_null(False))

        xmatched = xmatched.cte('xmatched', materialized=True)

        # We'll partition over each group of targets that match the same
        # catalogid and mark the one with the smallest distance to it as best.
        partition = (fn.first_value(xmatched.c.target_id)
                     .over(partition_by=[xmatched.c.catalogid],
                           order_by=[xmatched.c.distance.asc()]))
        best = peewee.Value(partition == xmatched.c.target_id)

        # Select the values to insert. Remove target_ids that were already
        # present in the relational table due to phase 1.
        insert_query = (xmatched
                        .select(xmatched.c.catalogid,
                                xmatched.c.target_id,
                                peewee.Value(self._version_id),
                                xmatched.c.distance,
                                best.alias('best'))
                        .where(xmatched.c.version_id == self._version_id)
                        # We do two NOT EXISTS instead of a single one to fully
                        # take advantage of the indexes on catalogid and
                        # (version_id, target_id).
                        .where(~fn.EXISTS(rel_model
                                          .select(SQL('1'))
                                          .where(rel_model.catalogid == xmatched.c.catalogid)))
                        .where(~fn.EXISTS(rel_model
                                          .select(SQL('1'))
                                          .where(rel_model.version_id == self._version_id,
                                                 rel_model.target_id == xmatched.c.target_id))))

        # TODO: this is probably non-optimal because we are doing a cross-match
        # of all of catalog against the model and catalog contains multiple
        # versions that we don't filter out until insert_query. But we cannot
        # add the filter in the CTE without breaking the index. A solution may
        # be to create a temporary table temp_catalog at the beginning of the
        # processing with a Q3C index and do everything there. Then we insert
        # all the rows in the real catalog table. This allows us to remove all
        # the filters on version_id because in the temporary table there is
        # only the current version being processed.

        # Insert into relational model.
        insert_query = rel_model.insert_from(
            insert_query, fields=[rel_model.catalogid,
                                  rel_model.target_id,
                                  rel_model.version_id,
                                  rel_model.distance,
                                  rel_model.best]).returning().with_cte(xmatched)

        with Timer() as timer:

            with self.database.atomic():

                self.log.debug('Inserting cross-matched targets into '
                               f'{rel_table_name}{self._get_sql(insert_query)}')

                n_catalogid = insert_query.execute(self.database)

        self.log.debug(f'Cross-matched {Catalog._meta.table_name} with '
                       f'{n_catalogid:,} targets in {table_name}. '
                       f'Run in {timer.interval:.3f} s.')

        self._analyze(rel_model)

    def _run_phase_3(self, model, rel_model):
        """Add non-matched targets to Catalog and the relational table."""

        meta = model._meta
        xmatch = meta.xmatch

        self.log.info('Phase 3: adding non cross-matched targets.')

        if 3 in xmatch.skip_phases:
            self.log.warning('Skipping due to configuration.')
            return

        table_name = meta.table_name

        rel_table_name = rel_model._meta.table_name
        catalog_table_name = Catalog._meta.table_name

        model_fields = self._get_model_fields(model)

        model_pk = meta.primary_key
        model_ra = meta.fields[xmatch.ra_column]
        model_dec = meta.fields[xmatch.dec_column]

        # Create a CTE for the sample region.
        sample = (model
                  .select(model_pk.alias('target_id'))
                  .where(self._get_sample_where(model_ra, model_dec)))

        if xmatch.has_missing_coordinates:
            sample = sample.where(model_ra.is_null(False), model_dec.is_null(False))

        sample = sample.cte('sample')

        # Get the max catalogid currently in the table for the version.
        if not self._max_cid:
            self.log.debug('Getting max. catalogid.')
            self._max_cid = Catalog.select(fn.MAX(Catalog.catalogid)).scalar() or 0

        unmatched = (sample
                     .select(fn.row_number().over() + self._max_cid,
                             sample.c.target_id,
                             peewee.Value(self._version_id),
                             peewee.SQL('true'))
                     .where(~fn.EXISTS(rel_model
                            .select(SQL('1'))
                            .where(rel_model.target_id == sample.c.target_id,
                                   rel_model.version_id == self._version_id))))

        with Timer() as timer:
            with self.database.atomic():

                rel_insert = rel_model.insert_from(
                    unmatched,
                    [rel_model.catalogid,
                     rel_model.target_id,
                     rel_model.version_id,
                     rel_model.best]).returning().with_cte(sample)

                self.log.debug(f'Inserting data into {rel_table_name}'
                               f'{self._get_sql(rel_insert)}')

                n_rows = rel_insert.execute()
                self.log.debug(f'Insertion into {rel_table_name} completed '
                               f'with {n_rows:,} rows in {timer.elapsed:.3f} s.')

                subq = (rel_model
                        .select(rel_model.catalogid,
                                *model_fields,
                                peewee.Value(table_name),
                                peewee.Value(self._version_id))
                        .where(~fn.EXISTS(Catalog.select(SQL('1'))
                                          .where(Catalog.catalogid == rel_model.catalogid)))
                        .join(model, on=(rel_model.target_id == model_pk)))

                insert_query = Catalog.insert_from(
                    subq,
                    [Catalog.catalogid,
                     Catalog.ra,
                     Catalog.dec,
                     Catalog.pmra,
                     Catalog.pmdec,
                     Catalog.parallax,
                     Catalog.lead,
                     Catalog.version_id]).returning()

                self.log.debug(f'Running INSERT query into {catalog_table_name}'
                               f'{self._get_sql(insert_query)}')
                n_rows = insert_query.execute()

        self._max_cid += n_rows  # Avoid having to calculate max_cid again
        self._catalog_count += n_rows

        self.log.debug(f'Inserted {n_rows:,} rows in {timer.elapsed:.3f} s.')

        self._analyze(rel_model, catalog=True)

    def _get_sql(self, query):
        """Returns coulourised SQL text for logging."""

        if self._options['show_sql']:
            return f': {color_text(str(query), "darkgrey")}'
        else:
            return '.'

    def _get_sample_where(self, ra_field, dec_field):
        """Returns the list of conditions to sample a model."""

        sample_region = self._options['sample_region']

        if sample_region is None:
            return True

        if len(sample_region) == 3 and not isinstance(sample_region[0], (list, tuple)):
            return fn.q3c_radial_query(ra_field,
                                       dec_field,
                                       sample_region[0],
                                       sample_region[1],
                                       sample_region[2])

        sample_conds = peewee.SQL('false')

        for ra, dec, radius in sample_region:
            sample_conds |= fn.q3c_radial_query(ra_field, dec_field, ra, dec, radius)

        return sample_conds

    def show_join_paths(self):
        """Prints all the available joint paths.

        This is useful before call `.run` to make sure the join paths to be
        used are correct or adjust the table weights otherwise. Note that
        the paths starting from the first model to be processed are ignored.

        """

        if len(self.process_order) == 1:
            return

        for table in self.process_order[1:]:
            paths = self.get_join_paths(table)
            if paths:
                for path in paths:
                    print(path)

    def _analyze(self, rel_model, catalog=False):
        """Analyses a relational model after insertion."""

        schema = rel_model._meta.schema
        table_name = rel_model._meta.table_name

        self.log.debug(f'Running ANALYZE on {table_name}.')
        vacuum_table(self.database, f'{schema}.{table_name}', vacuum=False, analyze=True)

        if catalog:
            cat_schema = Catalog._meta.schema
            cat_table_name = Catalog._meta.table_name
            self.log.debug(f'Running ANALYZE on {cat_table_name}.')
            vacuum_table(self.database, f'{cat_schema}.{cat_table_name}',
                         vacuum=False, analyze=True)

    def _log_table_configuration(self, model):
        """Logs the configuration used to cross-match a table."""

        xmatch = model._meta.xmatch

        parameters = ['ra_column', 'dec_column', 'pmra_column', 'pmdec_column',
                      'is_pmra_cos', 'parallax_column', 'epoch', 'epoch_column',
                      'epoch_format', 'has_duplicates', 'skip', 'skip_phases',
                      'query_radius', 'row_count', 'resolution', 'join_weight']

        self.log.debug('Table cross-matching parameters:')

        for parameter in parameters:

            value = getattr(xmatch, parameter, None)

            if parameter == 'query_radius':
                value = value or QUERY_RADIUS

            if value is True or value is False or value is None:
                value = str(value).lower()

            if parameter == 'row_count':
                self.log.debug(f'{parameter}: {value:,}')
            else:
                self.log.debug(f'{parameter}: {value}')

    def _get_larger_table(self, Catalog, Model):
        """Determines which table is larger."""

        # Use approximate count
        model_count = Model._meta.xmatch.row_count

        return 0 if self._catalog_count > model_count else 1

    def _get_catalog_row_count(self):
        """Determines catalog row count for a version."""

        # This requires the following function to have been defined in the
        # database.
        #
        # CREATE FUNCTION row_estimator(query text) RETURNS bigint
        #    LANGUAGE plpgsql AS
        # $$DECLARE
        #    plan jsonb;
        # BEGIN
        #    EXECUTE 'EXPLAIN (FORMAT JSON) ' || query INTO plan;
        #
        #    RETURN (plan->0->'Plan'->>'Plan Rows')::bigint;
        # END;$$;
        #
        #    RETURN (plan->0->'Plan'->>'Plan Rows')::bigint;
        # END;$$;

        self.log.debug(f'Determining row count for {Catalog._meta.table_name}.')

        self._catalog_count = (Select(columns=[
            fn.row_estimator(str(Catalog
                                 .select(SQL('*'))
                                 .where(Catalog.version_id == self._version_id)))
        ]).tuples().execute(self.database)[0][0])

        self.log.debug(f'{Catalog._meta.table_name} contains approximately '
                       f'{self._catalog_count} rows for '
                       f'version_id={self._version_id}.')
