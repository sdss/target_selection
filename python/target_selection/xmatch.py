#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-04-06
# @Filename: xmatch.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import copy
import hashlib
import inspect
import os
import re
import types
import uuid
import warnings

import networkx
import numpy
import peewee
import rich.markup
import yaml
from networkx.algorithms import shortest_path
from peewee import SQL, Case, Model, fn

from sdssdb.connection import PeeweeDatabaseConnection
from sdssdb.utils.internals import get_row_count
from sdsstools import merge_config
from sdsstools._vendor.color_print import color_text

import target_selection
from target_selection.exceptions import (
    TargetSelectionNotImplemented,
    TargetSelectionUserWarning,
    XMatchError,
)
from target_selection.utils import (
    Timer,
    get_configuration_values,
    get_epoch,
    is_view,
    sql_apply_pm,
    vacuum_outputs,
    vacuum_table,
)


EPOCH = 2016.0
QUERY_RADIUS = 1.0

#: Reserve last 11 bits for the run id.
RUN_ID_BIT_SHIFT = 64 - 11

TEMP_SCHEMA = "sandbox"


class Version(peewee.Model):
    """Model for the version table."""

    id = peewee.AutoField()
    plan = peewee.TextField()
    tag = peewee.TextField()

    class Meta:
        table_name = "version"


class Catalog(peewee.Model):
    """Model for the output table."""

    catalogid = peewee.BigIntegerField(index=True, null=False)
    iauname = peewee.TextField(null=True)
    ra = peewee.DoubleField(null=False)
    dec = peewee.DoubleField(null=False)
    pmra = peewee.FloatField(null=True)
    pmdec = peewee.FloatField(null=True)
    parallax = peewee.FloatField(null=True)
    lead = peewee.TextField(null=False)
    version_id = peewee.IntegerField(null=False, index=True)


class TempCatalog(Catalog):
    """Temporary output table."""

    catalogid = peewee.BigIntegerField(index=True, null=False)
    version_id = peewee.IntegerField(index=False)

    class Meta:
        primary_key = False


def XMatchModel(
    Model,
    resolution=None,
    ra_column=None,
    dec_column=None,
    pmra_column=None,
    pmdec_column=None,
    is_pmra_cos=True,
    parallax_column=None,
    epoch_column=None,
    epoch=None,
    epoch_format="jyear",
    relational_table=None,
    has_duplicates=False,
    has_missing_coordinates=False,
    skip=False,
    skip_phases=None,
    query_radius=None,
    join_weight=1,
    database_options=None,
):
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
    table_name : str
        Overrides the default model table name. This can be useful sometimes
        if, for example, a view has been created that contains only the
        columns from the main table needed for cross-matching.
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
    database_options : dict
        A dictionary of database configuration parameters to be set locally
        for this model for each processing phase transaction, temporarily
        overriding the default database configuration. Keys must be the
        database parameter to modify. The value can be a simple string with
        the value to set, or a dictionary that more accurately defines
        when the parameter will be applied. See `.XMatchPlanner` for more
        information.

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
            if "q3c" in index.sql.lower():
                match = re.match(r'.+q3c_ang2ipix\("*(\w+)"*, "*(\w+)"*\).+', index.sql)
                if match:
                    ra_column, dec_column = match.groups()

    meta.xmatch.ra_column = ra_column
    meta.xmatch.dec_column = dec_column

    meta.xmatch.pmra_column = pmra_column
    meta.xmatch.pmdec_column = pmdec_column
    meta.xmatch.is_pmra_cos = is_pmra_cos
    meta.xmatch.parallax_column = parallax_column

    assert ((epoch is None) & (epoch_column is None)) or (
        (epoch is not None) ^ (epoch_column is not None)
    ), "epoch and epoch_column are mutually exclusive."

    meta.xmatch.epoch = epoch
    meta.xmatch.epoch_column = epoch_column
    meta.xmatch.epoch_format = epoch_format

    meta.xmatch.relational_table = relational_table

    meta.xmatch.has_duplicates = has_duplicates
    meta.xmatch.has_missing_coordinates = has_missing_coordinates
    meta.xmatch.skip = skip
    meta.xmatch.skip_phases = skip_phases or []
    meta.xmatch.query_radius = query_radius

    meta.xmatch.row_count = int(
        get_row_count(
            meta.database,
            meta.table_name,
            schema=meta.schema,
            approximate=True,
        )
    )

    meta.xmatch.join_weight = join_weight

    meta.xmatch.database_options = database_options or {}

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
    ``version_id`` relates the record in the ``version`` table which contains
    the cross-matching plan and the tag of the code used when it was
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

    To speed things up unique targets are initially inserted into a temporary
    table ``<output_table>_<uid>`` where ``<output_table>`` is the name of
    the output table and ``<uid>`` is a unique identifier based on the version.

    Once the order has been determined and when `.run` is called, each table
    model is processed in order. The first model is just ingested completely
    into the temporary table and its associated relational table is created if
    it does not exist (the equivalent of *phase 3* below).

    For each additional model the following three stages are applied:

    - In *phase 1* we determine what targets in the input model have an
      existing cross-match to targets already ingested into the temporary
      table. To do that we build all possible joins between the model and the
      temporary table. If multiple joins are possible via a given table only
      the shortest is used (see `.get_join_paths`). For all matched targets we
      insert entries in the relational table. The *lead* of the original
      entries is not changed.

    - In *phase 2* we perform the actual cross-match between targets in the
      temporary table and the ones in the input catalogue. Currently the only
      cross-matching method available is a spatial cone query with radius
      ``query_radius``. All matched targets are added to the relational table
      and the one with the smallest distance is defined as *best*.

    - In *phase 3* we determine any target in the input catalogue that has not
      been cross-matched at this point and insert them into the temporary table
      as new entries. The *lead* is set to the input catalogue and the
      one-to-one match is added to the relational table.

    In phases 1 and 3 the queries are initially stored as a Postgresql
    temporary table for efficiency, and then copied to the relational table.
    After all the tables have been processed the output temporary table is
    inserted in bulk into the output table and dropped.

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
    plan : str
        The cross-matching plan version.
    run_id : int
        An integer to identify this run of cross-matching. The ID is bit
        shifted `.RUN_ID_BIT_SHIFT` positions and added to the catalogid.
        This allows to quickly associate a ``catalogid`` with a run without
        having to query ``catalogdb``. A ``run_id`` cannot be used if there
        are targets already in ``catalog`` using that same ``run_id``.
    version_id
        The ``catalogdb.version.id`` cross-match version to use. Normally this
        will be `None`, in which case a new version id will be created. For
        "addendum" runs, this should be set to the run to which to append.
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
        schemas are not supported), and the schema in which the output tables
        will be created.
    output_table : str
        The name of the output table. Defaults to ``catalog``.
    temp_schema
        The schema where the temporary ``catalog`` table will be initially created.
    log
        A logger to which to log messages. If not provided the ``target_selection``
        logger is used.
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
    path_mode : str
        The mode to look for join paths to link tables to output catalog.
        ``original`` mode uses a scorched earth algorithm that iteratively
        looks the shortest path and removes the first node trasspased and
        starts again until no shortest paths are available. ``full`` mode
        retrieves all the paths that are not a subsample of another path.
        ``config_list`` mode takes the list of paths from ``join_paths``
        parameter in the .yml configuration file.
    join_paths : list
        When using path_mode=``config_list`` is the list of paths to link
        tables to output catalog table in phase_1.
    database_options : dict
        A dictionary of database configuration parameters to be set locally
        during each phase transaction, temporarily overriding the default
        database configuration. Keys must be the database parameter to modify.
        The value can be a simple string with the value to set, or a dictionary
        that more accurately defines when the parameter will be applied. For
        example ::

            database_options:
                work_mem: '2GB'
                temp_buffers : {value: '500MB', phases: [3]}

        In this case the ``temp_buffers='500MB'`` option will only be set for
        phase 3. Configuration options to be used for a specific table can
        be set up when defining the `.XMatchModel`.

    """

    def __init__(
        self,
        database,
        models,
        plan,
        run_id,
        version_id=None,
        extra_nodes=[],
        order="hierarchical",
        key="row_count",
        epoch=EPOCH,
        start_node=None,
        query_radius=None,
        schema="catalogdb",
        temp_schema=TEMP_SCHEMA,
        output_table="catalog",
        log=None,
        log_path="./xmatch_{plan}.log",
        debug=False,
        show_sql=False,
        sample_region=None,
        database_options=None,
        path_mode="full",
        join_paths=None,
    ):
        self.log = log or target_selection.log
        self.log.header = ""

        if log_path:
            log_path = os.path.realpath(log_path)
            if self.log.fh:
                self.log.removeHandler(self.log.fh)
                self.log.fh = None
            self.log.start_file_logger(log_path.format(plan=plan), rotating=False, mode="a")

        if debug is True:
            self.log.sh.setLevel(0)
        elif debug is False:
            self.log.sh.setLevel(100)
        else:
            self.log.sh.setLevel(debug)

        self.schema = schema
        self.temp_schema = temp_schema

        self.output_table = output_table
        self.md5 = hashlib.md5(plan.encode()).hexdigest()[0:16]
        self._temp_table = self.output_table + "_" + self.md5

        self.database = database
        assert self.database.connected, "database is not connected."

        self.plan = plan
        self.run_id = run_id

        self.version_id = version_id
        self.is_addendum = version_id is not None

        self.tag = target_selection.__version__
        self.log.info(f"plan = {self.plan!r}; run_id = {self.run_id}; tag = {self.tag!r}.")
        self.log.info(f"Reference Epoch = {epoch}")
        self.models = {model._meta.table_name: model for model in models}
        self.extra_nodes = {model._meta.table_name: model for model in extra_nodes}
        self._check_models()

        self._options = {
            "query_radius": query_radius or QUERY_RADIUS,
            "show_sql": show_sql,
            "sample_region": sample_region,
            "epoch": epoch,
            "database_options": database_options or None,
        }

        qrd = self._options["query_radius"]
        self.log.info(f"Query radius = {qrd}")

        if self._options["sample_region"]:
            sample_region = self._options["sample_region"]
            self.log.warning(f"Using sample region {sample_region!r}.")

        self._log_db_configuration()

        self.model_graph = None
        self.update_model_graph()

        self.process_order = []
        self.set_process_order(order=order, key=key, start_node=start_node)

        self._max_cid = self.run_id << RUN_ID_BIT_SHIFT
        self.path_mode = path_mode

        if path_mode == "config_list":
            join_paths_warning = "join_paths needed for path_mode=config_list"
            assert join_paths is not None, join_paths_warning
            self.join_paths = join_paths

    @classmethod
    def read(cls, in_models, plan, config_file=None, **kwargs):
        """Instantiates `.XMatchPlanner` from a configuration file.

        The YAML configuration file must organised by plan string (multiple
        plans can live in the same file). Any parameter that
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

        It is also possible to use a parameter ``base_plan`` pointing to a
        previous plan string. In that case the previous plan configuration
        will be used as base and the new values will be merged (the update
        happens recursively as with normal Python dictionaries).

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
        plan : str
            The cross-matching plan.
        config_file : str or dict
            The path to the configuration file to use. Defaults to
            ``config/xmatch.yml``. The file must contain a hash with the
            cross-match plan.
        kwargs : dict
            User arguments that will override the configuration file values.

        """

        # HACK: this ensures that the catalogdb.models are populated. In principle thos would not
        # work if schema != catalogdb but anyway many other things would fail in that case.
        from sdssdb.peewee.sdss5db import catalogdb  # noqa

        if isinstance(in_models, (list, tuple)):
            models = in_models
        elif inspect.isclass(in_models) and issubclass(in_models, Model):
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
            raise TypeError(f"invalid input of type {type(in_models)!r}")

        assert database.connected, "database is not connected."

        if config_file is None:
            config_file = os.path.dirname(target_selection.__file__) + "/config/xmatch.yml"

        config = XMatchPlanner._read_config(config_file, plan)
        table_config = config.pop("tables", {}) or {}
        exclude_nodes = config.pop("exclude_nodes", []) or []
        assert "schema" in config, "schema is required in configuration."
        schema = config["schema"]

        models = {
            model._meta.table_name: model for model in models if model._meta.schema == schema
        }

        xmatch_models = {}
        for table_name in table_config:
            if table_name not in models:
                continue
            table_params = table_config[table_name] or {}
            xmatch_models[table_name] = XMatchModel(models[table_name], **table_params)

        extra_nodes_config = config.pop("extra_nodes", [])
        for table_name in extra_nodes_config:
            if table_name in models:
                xmatch_models[table_name] = XMatchModel(models[table_name])
            else:
                warnings.warn(
                    f"Cannot find model for extra node {table_name!r}.",
                    TargetSelectionUserWarning,
                )

        extra_nodes = [
            models[table_name]
            for table_name in models
            if table_name not in xmatch_models and table_name not in exclude_nodes
        ]

        config.update(kwargs)

        signature = inspect.signature(XMatchPlanner)

        valid_kw = {}
        for kw in config:
            if kw not in signature.parameters:
                warnings.warn(
                    f"ignoring invalid configuration value {kw!r}.",
                    TargetSelectionUserWarning,
                )
                continue
            valid_kw[kw] = config[kw]

        return cls(database, xmatch_models.values(), plan, extra_nodes=extra_nodes, **valid_kw)

    @staticmethod
    def _read_config(file_, plan):
        """Reads the configuration file, recursively."""

        if isinstance(file_, dict):
            config = copy.deepcopy(file_)
        else:
            config = yaml.load(open(file_, "r"), Loader=yaml.SafeLoader)

        assert plan in config, f"plan {plan!r} not found in configuration."

        base_plan = config[plan].pop("base_plan", None)
        if base_plan:
            config = merge_config(config[plan], XMatchPlanner._read_config(file_, base_plan))
        else:
            config = config[plan]

        return config

    def _check_models(self):
        """Checks the input analyse models."""

        catalog_tname = self.output_table

        # Remove extra nodes that are relational tables because we'll use
        # temporary relational tables instead.
        self.extra_nodes = {
            tname: self.extra_nodes[tname]
            for tname in self.extra_nodes
            if (
                not tname.startswith(catalog_tname + "_to_")
                and not tname == catalog_tname
                and self.extra_nodes[tname].table_exists()
            )
        }

        for tname in list(self.models.keys()):
            model = self.models[tname]
            meta = model._meta

            view_exists = any(
                [is_view(self.database, tname, self.schema, mat) for mat in [True, False]]
            )

            if not model.table_exists() and not view_exists:
                self.log.warning(f"table {tname!r} does not exist.")
            elif tname == catalog_tname:
                pass
            elif tname.startswith(catalog_tname + "_to_"):
                pass
            elif meta.xmatch.skip:
                self.extra_nodes[tname] = model
            else:
                continue

            self.models.pop(tname)

        if len(self.models) == 0:
            raise XMatchError("no models to cross-match.")

    def _log_db_configuration(self):
        """Logs some key database configuration parameters."""

        parameters = [
            "shared_buffers",
            "effective_cache_size",
            "wal_buffers",
            "effective_io_concurrency",
            "work_mem",
            "max_worker_processes",
            "random_page_cost",
            "seq_page_cost",
            "cpu_index_tuple_cost",
            "cpu_operator_cost",
            "default_statistics_target",
            "temp_buffers",
            "plan_cache_mode",
            "geqo_effort",
            "force_parallel_mode",
            "enable_seqscan",
            "enable_nestloop",
        ]

        values = get_configuration_values(self.database, parameters)

        self.log.debug("Current database configuration parameters.")
        for parameter in values:
            log_str = f"{parameter} = {values[parameter]}"
            self.log.debug(log_str)

    def update_model_graph(self):
        """Updates the model graph using models as nodes and fks as edges."""

        self.model_graph = networkx.Graph()
        self.model_graph.add_node(self._temp_table, model=TempCatalog)

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
                if hasattr(model._meta, "xmatch"):
                    model_weight = model._meta.xmatch.join_weight
                else:
                    model_weight = 1

                if hasattr(fk_model._meta, "xmatch"):
                    fk_model_weight = fk_model._meta.xmatch.join_weight
                else:
                    fk_model_weight = 1

                join_weight = 0.5 * (model_weight + fk_model_weight)

                self.model_graph.add_edge(table_name, ref_table_name, join_weight=join_weight)

            if model in self.models.values():
                rel_model = self.get_relational_model(model, sandboxed=False)
                rel_model._meta.schema = self.schema
                rel_model_tname = rel_model._meta.table_name
                self.model_graph.add_node(rel_model_tname, model=rel_model)
                self.model_graph.add_edge(self._temp_table, rel_model_tname)
                self.model_graph.add_edge(table_name, rel_model_tname)

        return self.model_graph

    def set_process_order(self, order="hierarchical", key="row_count", start_node=None):
        """Sets and returns the order in which tables will be processed.

        See `.XMatchPlanner` for details on how the order is decided depending
        on the input parameters.

        """

        if isinstance(order, (list, tuple)):
            self.log.info(f"processing order: {order}")
            self.process_order = order
            return order

        assert order in ["hierarchical", "global"], f"invalid order {order!r}."
        self.log.info(f"processing order mode is {order!r}")

        assert key in ["row_count", "resolution"], f"invalid key {key}."
        self.log.info(f"ordering key is {key!r}.")

        graph = self.model_graph.copy()

        for model in self.extra_nodes.values():
            graph.remove_node(model._meta.table_name)

        if order == "hierarchical":
            subgraphs = networkx.connected_components(graph)
        else:
            subgraphs = [node for node in graph.nodes]

        subgraphs_ext = []
        for sg in subgraphs:
            if start_node and start_node in sg:
                # Last item in record is 0 for initial table, 1 for other.
                # This prioritises the initial table in a sort without
                # having to reverse.
                subgraphs_ext.append((sg, numpy.nan, 0))
            else:
                if key == "row_count":
                    total_row_count = sum([self.models[tn]._meta.xmatch.row_count for tn in sg])
                    # Use -total_row_count to avoid needing reverse order.
                    subgraphs_ext.append((sg, -total_row_count, 1))
                elif key == "resolution":
                    resolution = [self.models[tn]._meta.xmatch.resolution for tn in sg]
                    if all(numpy.isnan(resolution)):
                        min_resolution = numpy.nan
                    else:
                        min_resolution = numpy.nanmin(resolution)
                    subgraphs_ext.append((sg, -(min_resolution or numpy.nan), 1))

        subgraphs_ordered = list(zip(*sorted(subgraphs_ext, key=lambda x: (x[2], x[1], x[0]))))[0]

        if order == "global":
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
                    if key == "row_count":
                        row_count = model._meta.xmatch.row_count
                        sg_ext.append((table_name, -row_count, 1))
                    elif key == "resolution":
                        resolution = model._meta.xmatch.resolution
                        sg_ext.append((table_name, resolution, 1))
                # Use table name as second sorting order
                #  to use alphabetic order in case of draw.
                sg_ordered = list(zip(*sorted(sg_ext, key=lambda x: (x[2], x[1], x[0]))))[0]
                ordered_tables.extend(sg_ordered)

        self.log.info(f"processing order: {ordered_tables}")
        self.process_order = ordered_tables

        return ordered_tables

    def get_join_paths(self, source, return_models=False, mode="full"):
        """Determines all possible join path between two tables.

        Mode ``original``  follows a scorched earth approach
        in which once an edge has been used for a join it cannot be used again.
        This produces only distinct joins between the two nodes.
        Mode ``full`` includes all possible paths that are not a subsample
        of another path.
        In both modes paths that include only the source and destination
        through their relational table are ignored, as are paths in which
        the last node before the output table has not yet been processed.
        Finally mode "config_list" is used to take the paths from a list
        indicated in the configuration file.

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
        mode : str
            The method used to obtain the paths. Possible values are
            ``original``, ``full``, and ``config_list``

        Returns
        -------
        `list`
            A list in which each item is a path with the table names or models
            joining ``source`` to ``dest``. The list is sorted in order of
            increasing path length.

        """

        graph = self.model_graph.copy()
        porder = self.process_order
        paths = []

        if mode == "original":
            dest = self._temp_table
            rel_table_name = self._temp_table + "_to_" + source
            if rel_table_name in graph:
                graph.remove_node(rel_table_name)

            while True:
                try:
                    spath = shortest_path(graph, source, dest, weight="join_weight")
                    if spath[-3] in porder and porder.index(spath[-3]) < porder.index(source):
                        paths.append(spath)
                    graph.remove_node(spath[1])
                except networkx.NetworkXNoPath:
                    break

        if mode == "full":
            ind = porder.index(source)
            for dest in porder[:ind]:
                paths1 = list(networkx.all_simple_paths(graph, source, dest, cutoff=1))
                if len(paths1) > 0:
                    paths += paths1
                    graph.remove_node(paths1[0][1])  # Only 1 possible path
            for dest in porder[:ind]:
                paths2orless = list(networkx.all_simple_paths(graph, source, dest, cutoff=2))
                paths2 = [el for el in paths2orless if len(el) == 3]
                if len(paths2) > 0:
                    paths += paths2
                    used2_edges = [(el[1], el[2]) for el in paths2]
                    graph.remove_edges_from(used2_edges)
            for dest in porder[:ind]:
                forb_couples = [set(el[1:]) for el in paths if len(el) == 3]
                paths3orless = list(networkx.all_simple_paths(graph, source, dest, cutoff=3))
                paths3 = [
                    el
                    for el in paths3orless
                    if len(el) == 4 and {el[1], el[3]} not in forb_couples
                ]
                if len(paths3) > 0:
                    paths += paths3

            full_paths = []
            for path in paths:
                full_path = path + ["catalog_to_" + path[-1], self._temp_table]
                full_paths.append(full_path)
            paths = full_paths

        if mode == "config_list":
            all_paths = self.join_paths
            paths = [path for path in all_paths if path[0] == source]

        if len(paths) == 0:
            return []

        if return_models:
            nodes = self.model_graph.nodes
            paths = [[nodes[node]["model"] for node in path] for path in paths]

        return paths

    def _prepare_models(self):
        """Prepare the Catalog, CatalogTemp, and Version models."""

        # Sets the metadata of the Catalog, TempCatalog, and Version tables.

        Catalog._meta.schema = self.schema
        Catalog._meta.table_name = self.output_table
        Catalog._meta.set_database(self.database)

        TempCatalog._meta.schema = self.temp_schema
        TempCatalog._meta.table_name = self._temp_table
        TempCatalog._meta.set_database(self.database)

        Version._meta.schema = self.schema
        Version._meta.set_database(self.database)

    def _check_version(self, model, force=False):
        """Checks if a model contains a plan version."""

        with self.database.atomic():
            self.database.execute_sql("SET LOCAL enable_seqscan = off;")
            vexists = (
                peewee.Select(
                    columns=[
                        fn.EXISTS(
                            model.select(SQL("1")).where(model.version_id == self.version_id)
                        )
                    ]
                )
                .tuples()
                .execute(self.database)
            )[0][0]

        if vexists:
            msg = (
                f"{model._meta.table_name!r} contains records for this "
                f"cross-matching plan ({self.plan})."
            )

            if force:
                self.log.warning(msg)
            else:
                raise XMatchError(msg)

    def _create_models(self, force=False):
        """Creates models and performs some checks."""

        self._prepare_models()

        # Bind models
        self.database.bind([Catalog, TempCatalog, Version])

        if not Version.table_exists():
            self.database.create_tables([Version])
            self.log.info(f"Created table {Version._meta.table_name}.")

        if self.version_id is None:
            version, vcreated = Version.get_or_create(plan=self.plan, tag=self.tag)
            self.version_id = version.id
        else:
            vcreated = False

        if vcreated:
            vmsg = "Added version record "
        else:
            vmsg = "Using version record "

        self.log.info(vmsg + f"({self.version_id}, {self.plan}, {self.tag}).")

        # Make sure the output table exists.
        if not Catalog.table_exists():
            # Add Q3C index for Catalog
            Catalog.add_index(
                SQL(
                    f"CREATE INDEX catalog_q3c_idx ON "
                    f"{self.schema}.{self.output_table} "
                    f"(q3c_ang2ipix(ra, dec))"
                )
            )

            self.database.create_tables([Catalog])

            self.log.info(f"Created table {self.output_table}.")

        # Check if Catalog already has entries for this xmatch version.
        if Catalog.table_exists():
            self._check_version(Catalog, force or self.is_addendum)

        if TempCatalog.table_exists():
            msg = f"Temporary table {self._temp_table} already exists."

            if force:
                self.log.warning(msg)
            else:
                raise XMatchError(msg)

            self._check_version(TempCatalog, force)

            try:
                self._temp_count = int(
                    get_row_count(
                        self.database,
                        self._temp_table,
                        schema=self.schema,
                        approximate=True,
                    )
                )
            except ValueError:
                self._temp_count = 0

        else:
            # Add Q3C index for TempCatalog
            TempCatalog.add_index(
                SQL(
                    f"CREATE INDEX IF NOT EXISTS "
                    f"{self._temp_table}_q3c_idx ON "
                    f"{self.temp_schema}.{self._temp_table} "
                    f"(q3c_ang2ipix(ra, dec))"
                )
            )

            self.database.create_tables([TempCatalog])

            self.log.info(f"Created table {self._temp_table}.")

            self._temp_count = 0

    def run(
        self, vacuum=False, analyze=False, from_=None, force=False, dry_run=False, keep_temp=False
    ):
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
        force : bool
            Allows to continue even if the temporary table exists or the
            output table contains records for this version. ``force=True`` is
            assumed if ``from_`` is defined.
        dry_run : bool
            If `False`, loads the temporary tables into ``Catalog`` and ``CatalogToXXX``.
            `True` implies ``keep_temp=True``; all the cross-matching steps will be run
            but the original tables won't be modified. A dry run can only be executed for
            a plan with a single catalogue since processing multiple catalogue requires the
            final tables to have been updated for successive catalogues.
        keep_temp : bool
            Whether to keep the temporary table or to drop it after the cross
            matching is done.

        """

        if len(self.process_order) > 1 and dry_run is True:
            raise RuntimeError("Cannot dry run with a plan that includes more than one catalogue.")

        if dry_run:
            keep_temp = True

        if vacuum or analyze:
            cmd = " ".join(("VACUUM" if vacuum else "", "ANALYZE" if analyze else "")).strip()
            self.log.info(f"Running {cmd} on output tables.")
            vacuum_outputs(
                self.database,
                vacuum=vacuum,
                analyze=analyze,
                schema=Catalog._meta.schema,
                table=Catalog._meta.table_name,
            )
            if TempCatalog.table_exists():
                vacuum_table(self.database, f"{self.schema}.{self._temp_table}")

        # Checks if there are multiple cross-matching plans running at the same
        # time. This is problematic because there can be catalogid collisions.
        # This should only be allowed if we define the starting catalogid
        # manually and make sure there cannot be collisions.
        # temp_tables = [table for table in self.database.get_tables(self.schema)
        #                if table.startswith(self.output_table + '_') and
        #                not table.startswith(self.output_table + '_to_') and
        #                table != self._temp_table]

        # if len(temp_tables) > 0:
        #     raise XMatchError('Another cross-match plan is currently running.')

        self._create_models(force or (from_ is not None))

        if from_ is not None:
            max_cid = TempCatalog.select(fn.MAX(TempCatalog.catalogid)).scalar()
            self._max_cid = max_cid + 1  # just to be sure numbers dont overlap

        with Timer() as timer:
            p_order = self.process_order
            for norder, table_name in enumerate(p_order):
                if dry_run and norder > 0:
                    raise RuntimeError("Cannot dry run more than one catalogue.")

                if from_ and p_order.index(table_name) < p_order.index(from_):
                    self.log.warning(f"Skipping table {table_name}.")
                    continue
                model = self.models[table_name]
                self.process_model(model, force=force)

                if not dry_run:
                    self.load_output_tables(model, keep_temp=keep_temp)

        self.log.info(f"Cross-matching completed in {timer.interval:.3f} s.")

    def process_model(self, model, force=False):
        """Processes a model, loading it into the output table."""

        table_name = model._meta.table_name
        self.log.header = f"[{table_name.upper()}] "

        self.log.info(f"Processing table {table_name}.")
        self._log_table_configuration(model)

        if model._meta.xmatch.has_duplicates:
            raise TargetSelectionNotImplemented(
                "handling of tables with duplicates is not implemented."
            )

        rel_model = self.get_relational_model(model, sandboxed=True, create=False)
        rel_model_table_name = rel_model._meta.table_name
        if rel_model.table_exists():
            if force is False:
                raise RuntimeError(
                    f"Sandboxed relational table {rel_model_table_name} exists. "
                    "Delete it manually before continuing."
                )
            else:
                self.log.warning(f"Sandboxed relational table {rel_model_table_name} exists.")

        # Check if there are already records in catalog for this version.
        if self.process_order.index(model._meta.table_name) == 0 and not self.is_addendum:
            is_first_model = True
        else:
            is_first_model = False

        self._phases_run = set()

        with Timer() as timer:
            if is_first_model:
                self._run_phase_3(model)
            else:
                self._run_phase_1(model)
                self._run_phase_2(model)
                self._run_phase_3(model)

        self.log.info(f"Fully processed {table_name} in {timer.elapsed:.0f} s.")

        self.update_model_graph()

        self.log.header = ""

    def _get_model_fields(self, model):
        """Returns the model fields needed to populate Catalog."""

        meta = model._meta
        xmatch = meta.xmatch
        fields = meta.fields

        # List of fields that will become the SELECT clause.
        model_fields = []

        ra_field = fields[xmatch.ra_column]
        dec_field = fields[xmatch.dec_column]

        to_epoch = self._options["epoch"]

        # RA, Dec, and proper motion fields.
        if model._meta.table_name == "tic_v8":
            # TODO: this should be handled in a way that can be opted-in from
            # the configuration, but for now I'll just hardcode it here.

            pmra_field = fields[xmatch.pmra_column]
            pmdec_field = fields[xmatch.pmdec_column]
            delta_years2000 = to_epoch - get_epoch(model)
            delta_years2015p5 = to_epoch - 2015.5

            racorr2000_field, deccorr2000_field = sql_apply_pm(
                ra_field,
                dec_field,
                pmra_field,
                pmdec_field,
                delta_years2000,
                xmatch.is_pmra_cos,
            )

            racorr2015p5_field, deccorr2015p5_field = sql_apply_pm(
                model.ra_orig,
                model.dec_orig,
                pmra_field,
                pmdec_field,
                delta_years2015p5,
                xmatch.is_pmra_cos,
            )

            ra_field = Case(
                None,
                [(model.posflag == "gaia2", racorr2015p5_field)],
                racorr2000_field,
            )
            dec_field = Case(
                None,
                [(model.posflag == "gaia2", deccorr2015p5_field)],
                deccorr2000_field,
            )

            model_fields.extend([ra_field.alias("ra"), dec_field.alias("dec")])
            model_fields.extend([pmra_field.alias("pmra"), pmdec_field.alias("pmdec")])

        elif xmatch.pmra_column:
            pmra_field = fields[xmatch.pmra_column]
            pmdec_field = fields[xmatch.pmdec_column]

            if (xmatch.epoch and xmatch.epoch != to_epoch) or xmatch.epoch_column:
                delta_years = to_epoch - get_epoch(model)
                ra_field, dec_field = sql_apply_pm(
                    ra_field,
                    dec_field,
                    pmra_field,
                    pmdec_field,
                    delta_years,
                    xmatch.is_pmra_cos,
                )
                if not xmatch.is_pmra_cos:
                    pmra_field *= fn.cos(fn.radians(dec_field))

            model_fields.extend([ra_field.alias("ra"), dec_field.alias("dec")])
            model_fields.extend([pmra_field.alias("pmra"), pmdec_field.alias("pmdec")])

        else:
            pmra_field = peewee.SQL("null")
            pmdec_field = peewee.SQL("null")

            model_fields.extend([ra_field.alias("ra"), dec_field.alias("dec")])

        # Parallax
        if xmatch.parallax_column:
            model_fields.append(fields[xmatch.parallax_column].alias("parallax"))

        return model_fields

    def get_output_model(self, temporary=False):
        """Returns the temporary or final output model (``catalog``)."""

        if temporary:
            return TempCatalog

        return Catalog

    def get_relational_model(self, model, sandboxed=False, temp=False, create=False):
        """Gets or creates a relational table for a given model.

        When the relational model is ``sandboxed``, the table is created in the
        temporary schema and suffixed with the same MD5 used for the run.

        """

        cat_table = Catalog._meta.table_name
        prefix = cat_table + "_to_"

        meta = model._meta
        pk = meta.primary_key

        if isinstance(pk, str) and pk == "__composite_key__":
            raise XMatchError(f"composite pk found for model {model.__name__!r}.")

        # Auto/Serial are automatically PKs. Convert them to integers
        # to avoid having two pks in the relational table.
        if pk.__class__.field_type == "AUTO":
            model_pk_class = peewee.IntegerField
        elif pk.__class__.field_type == "BIGAUTO":
            model_pk_class = peewee.BigIntegerField
        else:
            model_pk_class = pk.__class__

        class BaseModel(peewee.Model):
            catalogid = peewee.BigIntegerField(null=False, index=True)
            target_id = model_pk_class(null=False, index=True)
            version_id = peewee.SmallIntegerField(null=False, index=True)
            distance = peewee.DoubleField(null=True)
            best = peewee.BooleanField(null=False)
            plan_id = peewee.TextField(null=True)
            added_by_phase = peewee.SmallIntegerField(null=True)

            class Meta:
                database = meta.database
                schema = self.temp_schema if sandboxed else meta.schema
                primary_key = False

        model_prefix = "".join(x.capitalize() or "_" for x in prefix.rstrip().split("_"))

        RelationalModel = type(model_prefix + model.__name__, (BaseModel,), {})

        if temp:
            RelationalModel._meta.primary_key = False
            RelationalModel._meta.composite_key = False
            RelationalModel._meta.set_table_name(uuid.uuid4().hex[0:8])
            RelationalModel._meta.schema = None
            return RelationalModel

        table_name = prefix + meta.table_name
        if sandboxed:
            table_name += f"_{self.md5}"
        RelationalModel._meta.table_name = table_name

        if create and not RelationalModel.table_exists():
            RelationalModel.create_table()

        # Add foreign key field here. We want to avoid Peewee creating it
        # as a constraint and index if the table is created because that would
        # slow down inserts. We'll created them manually with add_fks.
        # Note that we do not create an FK between the relational model and
        # Catalog because the relationship is only unique on
        # (catalogid, version_id).
        RelationalModel._meta.remove_field("target_id")
        RelationalModel._meta.add_field(
            "target_id",
            peewee.ForeignKeyField(model, column_name="target_id", backref="+"),
        )

        return RelationalModel

    def _build_join(self, path):
        """Returns a build query for a given join path."""

        model = path[0]

        query = model.select()
        for inode in range(1, len(path)):
            if path[inode] is TempCatalog:
                query = query.join(
                    TempCatalog, on=(TempCatalog.catalogid == path[inode - 1].catalogid)
                )
            elif path[inode]._meta.table_name == "gaia_dr2_neighbourhood":
                query = (
                    query.join(path[inode])
                    .where(path[inode].angular_distance < 200)  # For 1-1
                    .distinct(path[inode].dr2_source_id)
                )  # To confirm
            else:
                query = query.join(path[inode])

        return query

    def _run_phase_1(self, model):
        """Runs the linking against matched catalogids stage."""

        table_name = model._meta.table_name
        rel_model_sb = self.get_relational_model(model, create=True, sandboxed=True)

        rel_model = self.get_relational_model(model, create=False, sandboxed=False)

        model_pk = model._meta.primary_key

        self.log.info("Phase 1: linking existing targets.")

        if 1 in model._meta.xmatch.skip_phases:
            self.log.warning("Skipping due to configuration.")
            return
        path_mode = self.path_mode
        join_paths = self.get_join_paths(table_name, mode=path_mode)

        if len(join_paths) == 0:
            self.log.debug(f"No paths found between {table_name} and temporary output table.")
            return False

        self.log.debug(
            f"Found {len(join_paths)} paths between {table_name} and temporary output table."
        )

        for n_path, path in enumerate(join_paths):
            # Remove the temporary catalog table at the end of the join path
            # because we only need catalogid and we can get that from the
            # relational model, saving us one join.
            path = path[0:-1]

            join_models = [self.model_graph.nodes[node]["model"] for node in path]

            # Get the relational model that leads to the temporary catalog
            # table in the join. We'll want to filter on version_id to avoid
            # a sequential scan.
            join_rel_model = join_models[-1]

            query = (
                self._build_join(join_models)
                .select(
                    model_pk.alias("target_id"),
                    join_rel_model.catalogid,
                    peewee.Value(True).alias("best"),
                )
                .where(join_rel_model.version_id == self.version_id, join_rel_model.best >> True)
                .where(
                    ~fn.EXISTS(
                        rel_model_sb.select(SQL("1")).where(
                            rel_model_sb.version_id == self.version_id,
                            (
                                (rel_model_sb.target_id == model_pk)
                                | (rel_model_sb.catalogid == join_rel_model.catalogid)
                            ),
                        )
                    )
                )
                # Select only one match per target in the catalogue with are cross-matching.
                .distinct(model_pk)
            )

            # Deal with duplicates in LS8
            if table_name == "legacy_survey_dr8":
                query = query.where(self._get_ls8_where(model))

            # Remove Duplicates and non-primary entries from LS10
            if table_name == "legacy_survey_dr10":
                query = query.where(
                    model.survey_primary >> True, fn.coalesce(model.ref_cat, "") != "T2"
                )

            # If the real relational model exists, exclude any matches that already exist there.
            if rel_model.table_exists():
                query = query.where(
                    ~fn.EXISTS(
                        rel_model.select(SQL("1")).where(
                            rel_model.version_id == self.version_id,
                            rel_model.target_id == model_pk,
                        )
                    )
                )

            # In query we do not include a Q3C where for the sample region because
            # TempCatalog for this plan should already be sample region limited.

            with Timer() as timer:
                with self.database.atomic():
                    temp_model = self.get_relational_model(model, temp=True, sandboxed=True)
                    temp_table = temp_model._meta.table_name

                    self._setup_transaction(model, phase=1)
                    self.log.debug(
                        f"Selecting linked targets into temporary "
                        f"table {temp_table!r} with join path "
                        f"{path}{self._get_sql(query)}"
                    )
                    query.create_table(temp_table, temporary=True)

                    self.log.debug(
                        f"Copying data into relational model "
                        f"{rel_model_sb._meta.table_name!r}."
                    )

                    fields = [
                        temp_model.target_id,
                        temp_model.catalogid,
                        temp_model.version_id,
                        temp_model.best,
                        temp_model.plan_id,
                        temp_model.added_by_phase,
                    ]

                    nids = (
                        rel_model_sb.insert_from(
                            temp_model.select(
                                temp_model.target_id,
                                temp_model.catalogid,
                                peewee.Value(self.version_id),
                                temp_model.best,
                                peewee.Value(self.plan) if self.is_addendum else None,
                                peewee.Value(1),
                            ),
                            fields,
                        )
                        .returning()
                        .execute()
                    )

            self.log.debug(f"Linked {nids.rowcount:,} records in {timer.interval:.3f} s.")
            self._phases_run.add(1)

            if nids.rowcount > 0:
                self._analyze(rel_model_sb)

    def _run_phase_2(self, model, source=TempCatalog):
        """Associates existing targets in Catalog with entries in the model.

        Here ``source`` is the catalogue with which we are spatially cross-matching.
        Normally this is the temporary catalog table which we are building for this
        cross-match run. But when we are doing an addendum, that table is going to be
        empty (at least of the first table in the addeundum), so we need to also have
        the option of using the real ``catalog`` table as the source. This method will
        call itself recursively with ``source=Catalog`` if this is an addendum run.

        """

        meta = model._meta
        xmatch = meta.xmatch

        table_name = meta.table_name

        s_table_name = source._meta.table_name
        self.log.info(f"Phase 2: cross-matching against existing targets ({s_table_name}).")

        if 2 in xmatch.skip_phases:
            self.log.warning("Skipping due to configuration.")
            return

        rel_model_sb = self.get_relational_model(model, create=True, sandboxed=True)
        rel_sb_table_name = rel_model_sb._meta.table_name

        rel_model = self.get_relational_model(model, create=False, sandboxed=False)

        model_pk = meta.primary_key
        model_ra = meta.fields[xmatch.ra_column]
        model_dec = meta.fields[xmatch.dec_column]

        catalog_epoch = self._options["epoch"]

        query_radius = xmatch.query_radius or self._options["query_radius"]

        # Should we use proper motions?
        model_epoch = get_epoch(model)

        is_model_expression = isinstance(model_epoch, (peewee.Expression, peewee.Function))
        use_pm = model_epoch and (is_model_expression or (model_epoch != catalog_epoch))

        if use_pm:
            self.log.debug("Determining maximum epoch delta " "between catalogues.")

            if isinstance(model_epoch, (int, float)):
                max_delta_epoch = float(abs(model_epoch - catalog_epoch))
            else:
                max_delta_epoch = float(
                    model.select(fn.MAX(fn.ABS(model_epoch - catalog_epoch)))
                    .where(self._get_sample_where(model_ra, model_dec))
                    .scalar()
                )

            max_delta_epoch += 0.1  # Add .1 yr to be sure it's an upper bound

            self.log.debug(f"Maximum epoch delta: {max_delta_epoch:.3f} (+ 0.1 year).")

        if use_pm:
            model_pmra = meta.fields[xmatch.pmra_column]
            model_pmdec = meta.fields[xmatch.pmdec_column]
            model_is_pmra_cos = int(xmatch.is_pmra_cos)

            q3c_dist = fn.q3c_dist_pm(
                model_ra,
                model_dec,
                model_pmra,
                model_pmdec,
                model_is_pmra_cos,
                model_epoch,
                source.ra,
                source.dec,
                catalog_epoch,
            )
            q3c_join = fn.q3c_join_pm(
                model_ra,
                model_dec,
                model_pmra,
                model_pmdec,
                model_is_pmra_cos,
                model_epoch,
                source.ra,
                source.dec,
                catalog_epoch,
                max_delta_epoch,
                query_radius / 3600.0,
            )
        else:
            q3c_dist = fn.q3c_dist(model_ra, model_dec, source.ra, source.dec)
            q3c_join = fn.q3c_join(
                model_ra, model_dec, source.ra, source.dec, query_radius / 3600.0
            )

        # Get the cross-matched catalogid and model target pk (target_id),
        # and their distance.
        xmatched = (
            source.select(
                source.catalogid,
                model_pk.alias("target_id"),
                q3c_dist.alias("distance"),
                source.version_id,
            )
            .join(model, peewee.JOIN.CROSS)
            .where(q3c_join)
            .where(self._get_sample_where(model_ra, model_dec))
        )

        if table_name == "legacy_survey_dr8":
            xmatched = xmatched.where(self._get_ls8_where(model))
        if table_name == "legacy_survey_dr10":
            xmatched = xmatched.where(
                model.survey_primary >> True, fn.coalesce(model.ref_cat, "") != "T2"
            )

        # This may break the use of the index but I think it's needed if
        # the model is the second table in q3c_join and has empty RA/Dec.
        if xmatch.has_missing_coordinates and use_pm:
            xmatched = xmatched.where(model_ra.is_null(False), model_dec.is_null(False))

        xmatched = xmatched.cte("xmatched", materialized=True)

        # We'll partition over each group of targets that match the same
        # catalogid and mark the one with the smallest distance to it as best.
        partition = fn.first_value(xmatched.c.target_id).over(
            partition_by=[xmatched.c.catalogid], order_by=[xmatched.c.distance.asc()]
        )
        best = peewee.Value(partition == xmatched.c.target_id)

        # Select the values to insert. Remove target_ids that were already
        # present in the relational table due to phase 1.
        # We separate the filter in two IF NOT EXISTS clauses to
        # be sure the query planner uses the indexes for each (it
        # won't necessarily do it if we do an OR). Also, make sure
        # we compare version_id, target_id and not the other way
        # because that's the order in which we defined the index.
        in_query = xmatched.select(
            xmatched.c.target_id,
            xmatched.c.catalogid,
            peewee.Value(self.version_id).alias("version_id"),
            xmatched.c.distance.alias("distance"),
            best.alias("best"),
            self.plan if self.is_addendum else None,
            peewee.Value(2).alias("added_by_phase"),
        )

        # We only need to care about already linked targets if phase 1 run.
        if 1 in self._phases_run:
            in_query = in_query.where(
                xmatched.c.version_id == self.version_id,
                ~fn.EXISTS(
                    rel_model_sb.select(SQL("1")).where(
                        (rel_model_sb.version_id == self.version_id)
                        & (
                            (rel_model_sb.catalogid == xmatched.c.catalogid)
                            | (rel_model_sb.target_id == xmatched.c.target_id)
                        )
                    )
                ),
            )

        if rel_model.table_exists():
            in_query = in_query.where(
                ~fn.EXISTS(
                    rel_model.select(SQL("1")).where(
                        (rel_model.version_id == self.version_id)
                        & (
                            (rel_model.catalogid == xmatched.c.catalogid)
                            | (rel_model.target_id == xmatched.c.target_id)
                        )
                    )
                )
            )

        with Timer() as timer:
            with self.database.atomic():
                # 1. Tweak database configuration for this transaction to
                #    ensure Q3C index is used.

                # May need to increase work_mem during this transaction to
                # make sure the Q3C index is used.
                self._setup_transaction(model, phase=2)

                # 2. Run cross-match and insert data into relational model.

                fields = [
                    rel_model_sb.target_id,
                    rel_model_sb.catalogid,
                    rel_model_sb.version_id,
                    rel_model_sb.distance,
                    rel_model_sb.best,
                    rel_model_sb.plan_id,
                    rel_model_sb.added_by_phase,
                ]

                in_query = rel_model_sb.insert_from(
                    in_query.with_cte(xmatched),
                    fields,
                ).returning()

                self.log.debug(
                    f"Running Q3C query and inserting cross-matched data into "
                    f"relational table {rel_sb_table_name!r}: "
                    f"{self._get_sql(in_query)}"
                )

                n_catalogid = in_query.execute().rowcount

        self.log.debug(
            f"Cross-matched {source._meta.table_name} with "
            f"{n_catalogid:,} targets in {table_name}. "
            f"Run in {timer.interval:.3f} s."
        )

        if n_catalogid > 0:
            self._phases_run.add(2)
            self._analyze(rel_model_sb)

        # For addenda it's not sufficient to cross-match with the temporary table, because that
        # does not contain all the cumulated targets from this cross-match version. We need to
        # also spatially cross-match with the real catalog table (but only for the targets with
        # version_id=<this-version-id>).
        if self.is_addendum and source != Catalog:
            self._run_phase_2(model, source=Catalog)

    def _run_phase_3(self, model):
        """Add non-matched targets to Catalog and the relational table."""

        meta = model._meta
        xmatch = meta.xmatch

        self.log.info("Phase 3: adding non cross-matched targets.")

        rel_model_sb = self.get_relational_model(model, create=True, sandboxed=True)
        rel_sb_table_name = rel_model_sb._meta.table_name

        rel_model = self.get_relational_model(model, create=False, sandboxed=False)

        if 3 in xmatch.skip_phases:
            self.log.warning("Skipping due to configuration.")
            return

        table_name = meta.table_name

        model_fields = self._get_model_fields(model)

        model_pk = meta.primary_key
        model_ra = meta.fields[xmatch.ra_column]
        model_dec = meta.fields[xmatch.dec_column]

        unmatched = model.select(
            (fn.row_number().over() + self._max_cid).alias("catalogid"),
            model_pk.alias("target_id"),
            *model_fields,
        ).where(self._get_sample_where(model_ra, model_dec))

        if 1 in self._phases_run or 2 in self._phases_run:
            unmatched = unmatched.where(
                ~fn.EXISTS(
                    rel_model_sb.select(SQL("1")).where(
                        rel_model_sb.version_id == self.version_id,
                        rel_model_sb.target_id == model_pk,
                        rel_model_sb.best >> True,
                    )
                )
            )

        if rel_model.table_exists():
            unmatched = unmatched.where(
                ~fn.EXISTS(
                    rel_model.select(SQL("1")).where(
                        rel_model.version_id == self.version_id,
                        rel_model.target_id == model_pk,
                        rel_model.best >> True,
                    )
                )
            )

        if xmatch.has_missing_coordinates:
            unmatched = unmatched.where(model_ra.is_null(False), model_dec.is_null(False))

        # TODO: this is horrible and should be moved to the configuration.
        if model._meta.table_name == "tic_v8":
            unmatched = unmatched.where(model.objtype != "EXTENDED")

        if table_name == "legacy_survey_dr8":
            unmatched = unmatched.where(self._get_ls8_where(model))

        if table_name == "legacy_survey_dr10a":
            unmatched = unmatched.where(model.ref_cat != "T2")

        if table_name == "legacy_survey_dr10":
            unmatched = unmatched.where(
                model.survey_primary >> True, fn.coalesce(model.ref_cat, "") != "T2"
            )

        with Timer() as timer:
            with self.database.atomic():
                # TODO: Not sure it's worth using a temporary table here.

                # 1. Run link query and create temporary table with results.

                self._setup_transaction(model, phase=3)

                temp_model = self.get_relational_model(model, temp=True, sandboxed=True)
                temp_model_name = temp_model._meta.table_name

                self.log.debug(
                    f"Selecting unique targets "
                    f"into temporary table "
                    f"{temp_model_name!r}{self._get_sql(unmatched)}"
                )

                unmatched.create_table(temp_model_name, temporary=True)

                # Analyze the temporary table to gather stats.
                # self.log.debug('Running ANALYZE on temporary table.')
                # self.database.execute_sql(f'ANALYZE "{temp_model_name}";')

                # 2. Copy data from temporary table to relational table. Add
                #    catalogid at this point.

                fields = [
                    temp_model.catalogid,
                    temp_model.target_id,
                    temp_model.version_id,
                    temp_model.best,
                    rel_model_sb.plan_id,
                    rel_model_sb.added_by_phase,
                ]

                rel_insert_query = rel_model_sb.insert_from(
                    temp_model.select(
                        temp_model.catalogid,
                        temp_model.target_id,
                        self.version_id,
                        peewee.SQL("true"),
                        self.plan if self.is_addendum else None,
                        peewee.Value(3).alias("added_by_phase"),
                    ),
                    fields,
                ).returning()

                self.log.debug(
                    f"Copying data into relational model "
                    f"{rel_sb_table_name!r}"
                    f"{self._get_sql(rel_insert_query)}"
                )

                cursor = rel_insert_query.execute()
                n_rows = cursor.rowcount

                self.log.debug(
                    f"Insertion into {rel_sb_table_name} completed "
                    f"with {n_rows:,} rows in "
                    f"{timer.elapsed:.3f} s."
                )

                # 3. Fill out the temporary catalog table with the information
                #    from the unique targets.

                temp_table = peewee.Table(temp_model_name)

                fields = [TempCatalog.catalogid, TempCatalog.lead, TempCatalog.version_id]
                select_columns = [
                    temp_table.c.catalogid,
                    peewee.Value(table_name),
                    self.version_id,
                ]
                for field in model_fields:
                    if field._alias == "ra":
                        fields.append(TempCatalog.ra)
                        select_columns.append(temp_table.c.ra)
                    elif field._alias == "dec":
                        fields.append(TempCatalog.dec)
                        select_columns.append(temp_table.c.dec)
                    elif field._alias == "pmra":
                        fields.append(TempCatalog.pmra)
                        select_columns.append(temp_table.c.pmra)
                    elif field._alias == "pmdec":
                        fields.append(TempCatalog.pmdec)
                        select_columns.append(temp_table.c.pmdec)
                    elif field._alias == "parallax":
                        fields.append(TempCatalog.parallax)
                        select_columns.append(temp_table.c.parallax)

                insert_query = TempCatalog.insert_from(
                    temp_table.select(*select_columns), fields
                ).returning()

                self.log.debug(
                    f"Running INSERT query into {self._temp_table}"
                    f"{self._get_sql(insert_query)}"
                )

                cursor = insert_query.execute()
                n_rows = cursor.rowcount

        # Avoid having to calculate max_cid again
        self._max_cid += n_rows
        self._temp_count += n_rows

        self.log.debug(f"Inserted {n_rows:,} rows. Total time: {timer.elapsed:.3f} s.")

        self._phases_run.add(3)

        if n_rows > 0.5 * self._temp_count:  # Cluster if > 50% of rows are new
            self.log.debug(f"Running CLUSTER on {self._temp_table} with q3c index.")
            self.database.execute_sql(
                f"CLUSTER {self.temp_schema}.{self._temp_table} "
                f"using {self._temp_table}_q3c_idx;"
            )
            self.log.debug(f"Running ANALYZE on {self._temp_table}.")
            self.database.execute_sql(f"ANALYZE {self.temp_schema}.{self._temp_table};")

        self._analyze(rel_model_sb, catalog=False)

    def load_output_tables(self, model, keep_temp=False):
        """Loads the temporary tables into the output tables."""

        self._load_output_table(TempCatalog, Catalog, keep_temp=keep_temp)

        rel_model_sb = self.get_relational_model(model, sandboxed=True, create=False)
        rel_model = self.get_relational_model(model, sandboxed=False, create=True)
        self._load_output_table(rel_model_sb, rel_model, keep_temp=keep_temp)

    def _load_output_table(self, from_model, to_model, keep_temp=False):
        """Copies the temporary table to the real output table."""

        to_table = f"{to_model._meta.schema}.{to_model._meta.table_name}"
        from_table = f"{from_model._meta.schema}.{from_model._meta.table_name}"

        self.log.info(f"Copying {from_table} to {to_table}.")

        with Timer() as timer:
            with self.database.atomic():
                self._setup_transaction()

                insert_query = to_model.insert_from(
                    from_model.select(), from_model.select()._returning
                ).returning()

                self.log.debug(
                    f"Running INSERT query into {to_table} {self._get_sql(insert_query)}"
                )

                cursor = insert_query.execute()
                n_rows = cursor.rowcount

        self.log.debug(f"Inserted {n_rows:,} rows in {timer.elapsed:.3f} s.")

        if not keep_temp:
            self.database.drop_tables([from_model])
            self.log.info(f"Dropped temporary table {from_table}.")

        if n_rows > 0:
            self.log.debug(f"Running VACUUM ANALYZE on {to_table}.")
            vacuum_table(self.database, to_table, vacuum=True, analyze=True)

    def _get_sql(self, query, return_string=False):
        """Returns colourised SQL text for logging."""

        query_str, query_params = query.sql()

        if query_params:
            for ind in range(len(query_params)):
                if isinstance(query_params[ind], str):
                    query_params[ind] = "'" + query_params[ind] + "'"

            query_str = query_str % tuple(query_params)

        query_str = query_str.replace("None", "Null")
        if return_string:
            return query_str
        elif self.log.rich_console:
            return f": {rich.markup.escape(query_str)}"
        elif self._options["show_sql"]:
            return f': {color_text(query_str, "blue")}'
        else:
            return "."

    def _setup_transaction(self, model=None, phase=None):
        """Sets database parameters for the transaction."""

        if not self._options["database_options"]:
            return

        options = self._options["database_options"].copy()
        if model and model._meta.xmatch.database_options:
            options.update(model._meta.xmatch.database_options)

        for param in options:
            if param == "maintenance_work_mem":
                continue
            param_config = options[param]
            if isinstance(param_config, dict):
                value = param_config.get("value")
                param_phases = param_config.get("phase", None)
                if param_phases and (not phase or phase not in param_phases):
                    continue
            else:
                value = param_config
            stm = f"SET LOCAL {param}={value!r};"
            self.database.execute_sql(stm)
            self.log.debug(stm)

    def _get_sample_where(self, ra_field, dec_field):
        """Returns the list of conditions to sample a model."""

        sample_region = self._options["sample_region"]

        if sample_region is None:
            return True

        if len(sample_region) == 3 and not isinstance(sample_region[0], (list, tuple)):
            return fn.q3c_radial_query(
                ra_field,
                dec_field,
                sample_region[0],
                sample_region[1],
                sample_region[2],
            )

        sample_conds = peewee.SQL("false")

        for ra, dec, radius in sample_region:
            sample_conds |= fn.q3c_radial_query(ra_field, dec_field, ra, dec, radius)

        return sample_conds

    def _get_ls8_where(self, model):
        """Removes duplicates from LS8 queries."""

        return ~(
            (
                (model.release == 8000)
                & (model.dec > 32.375)
                & (model.ra > 100.0)
                & (model.ra < 300.0)
            )
            | ((model.release == 8001) & (model.dec < 32.375))
        )

    def show_join_paths(self):
        """Prints all the available joint paths.

        This is useful before call `.run` to make sure the join paths to be
        used are correct or adjust the table weights otherwise. Note that
        the paths starting from the first model to be processed are ignored.

        """
        if len(self.process_order) == 1:
            return

        mode = self.path_mode
        for table in self.process_order[1:]:
            paths = self.get_join_paths(table, mode=mode)
            if paths:
                for path in paths:
                    print(path)

    def _analyze(self, rel_model, vacuum=False, catalog=False):
        """Analyses a relational model after insertion."""

        schema = rel_model._meta.schema
        table_name = rel_model._meta.table_name

        db_opts = self._options["database_options"]
        if db_opts:
            work_mem = db_opts.get("maintenance_work_mem", None)
            if work_mem:
                self.database.execute_sql(f"SET maintenance_work_mem = {work_mem!r}")

        self.log.debug(f"Running ANALYZE on {table_name}.")
        vacuum_table(self.database, f"{schema}.{table_name}", vacuum=vacuum, analyze=True)

        if catalog:
            self.log.debug(f"Running ANALYZE on {self._temp_table}.")
            vacuum_table(
                self.database,
                f"catalogdb.{self._temp_table}",
                vacuum=vacuum,
                analyze=True,
            )

    def _log_table_configuration(self, model):
        """Logs the configuration used to cross-match a table."""

        xmatch = model._meta.xmatch

        parameters = [
            "ra_column",
            "dec_column",
            "pmra_column",
            "pmdec_column",
            "is_pmra_cos",
            "parallax_column",
            "epoch",
            "epoch_column",
            "epoch_format",
            "has_duplicates",
            "skip",
            "skip_phases",
            "query_radius",
            "row_count",
            "resolution",
            "join_weight",
            "has_missing_coordinates",
        ]

        self.log.debug("Table cross-matching parameters:")

        for parameter in parameters:
            value = getattr(xmatch, parameter, None)

            if parameter == "query_radius":
                value = value or QUERY_RADIUS

            if value is True or value is False or value is None:
                value = str(value).lower()

            if parameter == "row_count":
                self.log.debug(f"{parameter}: {value:,}")
            else:
                self.log.debug(f"{parameter}: {value}")
