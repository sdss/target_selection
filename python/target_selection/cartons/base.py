#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-24
# @Filename: base.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import abc
import inspect
import warnings

import numpy
import peewee
from astropy import table

from sdssdb.peewee import BaseModel
from sdssdb.peewee.sdss5db import catalogdb as cdb
from sdssdb.peewee.sdss5db import targetdb as tdb
from sdsstools import read_yaml_file
from sdsstools.color_print import color_text

from target_selection import __version__, config, log
from target_selection.exceptions import (TargetSelectionError,
                                         TargetSelectionUserWarning)
from target_selection.utils import Timer


EPOCH = 2015.5


class BaseCarton(metaclass=abc.ABCMeta):
    """A base class for target cartons.

    This class is not intended for direct instantiation. Instead, it must be
    subclassed and the relevant class attributes overridden with the values
    corresponding to the carton.

    Parameters
    ----------
    targeting_plan : str
        The target selection plan version.
    config_file : str
        The path to the configuration file to use. If undefined, uses the
        internal ``target_selection.yml`` file.
    schema : str
        Schema in which the temporary table with the results of the
        query will be created. If `None`, tries to use the ``schema`` parameter
        from the configuration file for this plan of target selection. If
        the parameter is not set, defaults to ``'sandbox'``.
    table_name : str
        The name of the temporary table. Defaults to ``temp_<name>`` where
        ``<name>`` is the target class name.

    Attributes
    ----------
    name : str
        The name of the carton (required).
    cadence : str
        The label of the cadence rule for this carton.
    category : str
        The category of targets for this carton.
    mapper : str
        The mapper with which this carton is associated.
    orm : str
        The ORM library to be used, ``peewee`` or ``sqlalchemy``.
    tag : str
        The version of the ``target_selection`` code used.
    query_region : tuple
        A tuple defining the region over which the query should be performed,
        with the format ``(ra, dec, radius)`` in degrees. This will append a
        ``q3c_radial_query`` condition to the query.
    load_magnitudes : bool
        Whether to load target magnitudes. In general this must be `True`
        except for cartons for which it's known the magnitudes will not be
        used, e.g., skies.

    """

    name = None
    cadence = None
    category = None
    program = None
    mapper = None
    priority = None

    load_magnitudes = True

    query_region = None

    def __init__(self, targeting_plan, config_file=None,
                 schema=None, table_name=None):

        assert self.name, 'carton subclass must override name'
        assert self.category, 'carton subclass must override category'
        assert self.program, 'carton subclass must override program'

        self.plan = targeting_plan
        self.tag = __version__

        if config_file:
            this_config = read_yaml_file(config_file)
        else:
            this_config = config

        if self.plan not in this_config:
            raise TargetSelectionError(f'({self.name}): cannot find plan '
                                       f'{self.plan!r} in config.')

        self.config = this_config[self.plan]

        if 'parameters' in self.config:
            self.parameters = self.config['parameters'].get(self.name, None)
        else:
            self.parameters = None

        try:
            self.xmatch_plan = self.config['xmatch_plan']
        except KeyError:
            raise TargetSelectionError(f'({self.name}): xmatch_plan '
                                       'not found in config.')

        # Check the signature of build_query
        self._build_query_signature = inspect.signature(self.build_query)
        if 'version_id' not in self._build_query_signature.parameters:
            raise TargetSelectionError('build_query does not '
                                       'accept version_id')

        self.database = tdb.database
        assert self.database.connected, 'database is not connected.'

        self.schema = schema or self.config.get('schema', None) or 'sandbox'
        self.table_name = table_name or f'temp_{self.name}'

        if self.cadence:
            assert (tdb.Cadence
                    .select()
                    .where(tdb.Cadence.label == self.cadence)
                    .count() == 1), (f'{self.cadence!r} does not '
                                     'exist in targetdb.cadence.')

        self.log = log
        self.has_run = False

    @property
    def path(self):
        """The schema-qualified path to the output table."""

        return self.schema + '.' + self.table_name

    @abc.abstractmethod
    def build_query(self, version_id, query_region=None):
        """Builds and returns the query.

        The ORM query for the target class. Note that this must be the
        *un-executed* query, which will be executed in `.run`. The select
        statement must be run on ``catalogdb.catalog`` and needs to include,
        at least, ``catalogid``. Additional columns such as
        ``ra, dec, pmra, pmdec`` can be specified but are otherwise completed
        when calling `.run`. Magnitude columns can be included and propagated
        to ``targetdb`` but must be aliased as ``magnitude_<band>`` where
        ``<band>`` must be one of the columns in ``targetdb.magnitude``.

        The query returned must include a filter for ``version_id`` on
        ``catalogdb.catalog``. Normally this is applied to the ``Catalog``
        select as ``.where(Catalog.version_id == version_id)``.

        Parameters
        ----------
        version_id : int
            The id of the cross-match version to use.

        Returns
        -------
        query
            A :class:`peewee:Select` or :class:`peewee:ModelSelect` query.

        """

        pass

    def get_model(self):
        """Returns a Peewee model for the temporary table using reflection."""

        class Model(BaseModel):

            catalogid = peewee.BigIntegerField(primary_key=True)
            selected = peewee.BooleanField()
            cadence = peewee.TextField(null=True)
            priority = peewee.IntegerField()

            class Meta:
                database = self.database
                table_name = self.table_name
                schema = self.schema
                reflection_options = {'skip_foreign_keys': True}
                use_reflection = True

        if not Model.table_exists():
            raise TargetSelectionError(f'temporary table {self.path!r} does '
                                       'not exist.')

        return Model

    def run(self, tile=None, query_region=None, overwrite=False,
            **post_process_kawrgs):
        """Executes the query and post-process steps, and stores the results.

        This method calls `.build_query` and runs the returned query. The
        output of the query is stored in a temporary table whose schema and
        table name are defined when the object is instantiated.

        After the query has run, the `.post_process` routine is called if the
        method has been overridden for the given carton.

        Parameters
        ----------
        query_region : tuple
            A tuple defining the region over which the query should be
            performed, with the format ``(ra, dec, radius)`` in degrees. This
            will append a ``q3c_radial_query`` condition to the query.
        overwrite : bool
            Whether to overwrite the intermediary table if already exists.
        post_process_args : dict
            Keyword arguments to be passed to `.post_process`.

        Returns
        -------
        model : :class:`peewee:Model`
            The model for the intermediate table.

        """

        query_region = query_region or self.query_region

        if self.database.table_exists(self.table_name, schema=self.schema):
            if overwrite:
                log.info(f'Dropping table {self.path!r}.')
                self.drop_table()
            else:
                raise RuntimeError(f'Temporary table {self.path!r} '
                                   'already exists.')

        log.info('Running query ...')
        version_id = cdb.Version.get(plan=self.xmatch_plan).id

        # If build_query accepts a query_region parameter, call with the query
        # region. Otherwise will add the radial query condition later.
        if 'query_region' in self._build_query_signature.parameters:
            query = self.build_query(version_id, query_region=query_region)
        else:
            query = self.build_query(version_id)

        # Make sure the catalogid column is selected.
        if cdb.Catalog.catalogid not in query._returning:
            raise RuntimeError('catalogid is not being returned in query.')

        if query_region:
            if 'query_region' in self._build_query_signature.parameters:
                pass
            else:
                # This may be quite inefficient depending on the query.
                subq = query.alias('subq')
                query = (peewee.Select(columns=[peewee.SQL('subq.*')])
                         .from_(subq)
                         .join(cdb.Catalog,
                               on=(cdb.Catalog.catalogid == subq.c.catalogid))
                         .where(peewee.fn.q3c_radial_query(cdb.Catalog.ra,
                                                           cdb.Catalog.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2])))

        query_sql, params = query.sql()
        cursor = self.database.cursor()
        query_str = cursor.mogrify(query_sql, params).decode()

        log.debug(color_text(f'CREATE TABLE IF NOT EXISTS {self.path} AS ' +
                             query_str, 'darkgrey'))

        with self.database.atomic():
            with Timer() as timer:
                self.setup_transaction()
                self.database.execute_sql(f'CREATE TABLE IF NOT EXISTS '
                                          f'{self.path} AS ' + query_sql,
                                          params)

        log.info(f'Created table {self.path!r} in {timer.interval:.3f} s.')

        log.debug('Adding columns and indexes.')

        columns = [col.name
                   for col in self.database.get_columns(self.table_name,
                                                        self.schema)]

        if 'selected' not in columns:
            self.database.execute_sql(f'ALTER TABLE {self.path} '
                                      'ADD COLUMN selected BOOL DEFAULT TRUE;')

        if 'cadence' not in columns:
            self.database.execute_sql(f'ALTER TABLE {self.path} '
                                      'ADD COLUMN cadence VARCHAR DEFAULT NULL;')
        if 'priority' not in columns:
            self.database.execute_sql(f'ALTER TABLE {self.path} '
                                      'ADD COLUMN priority INTEGER '
                                      'DEFAULT NULL;')

        self.database.execute_sql(f'ALTER TABLE {self.path} '
                                  'ADD PRIMARY KEY (catalogid);')
        self.database.execute_sql(f'CREATE INDEX ON {self.path} (selected);')
        self.database.execute_sql(f'ANALYZE {self.path};')

        ResultsModel = self.get_model()

        n_rows = ResultsModel.select().count()
        log.debug(f'Table {self.path!r} contains {n_rows:,} rows.')

        log.debug('Running post-process.')
        with self.database.atomic():
            self.setup_transaction()
            self.post_process(ResultsModel, **post_process_kawrgs)

        self.has_run = True

        return ResultsModel

    def post_process(self, model, **kwargs):
        """Post-processes the temporary table.

        This method provides a framework for applying non-SQL operations on
        carton query. It receives the model for the temporary table and can
        perform any operation on it, including modifying the ``selected``
        column with a mask of targets to be used.

        This method can also be used to set the ``cadence`` column in the
        temporary table. This column will be used to set the target cadence if
        the carton `.cadence` attribute is not set.

        `.post_process` runs inside a database transaction so it's not
        necessary to create a new one, but savepoints can be added.

        Parameters
        ----------
        model : peewee:Model
            The model of the intermediate table.

        Returns
        -------
        mask : `tuple`
            The list of catalogids from the temporary table that should be
            selected as part of this carton. If `True` (the default), selects
            all the records.

        """

        return True

    def setup_transaction(self):
        """Setups the transaction locally modifying the datbase parameters.

        This method runs inside a transaction and can be overridden to set the
        parameters of the transaction manually. It applies to both `.run` and
        `.load`.

        """

        if 'database_options' not in self.config:
            return

        for param, value in self.config['database_options'].items():
            self.database.execute_sql(f'SET LOCAL {param} = {value!r};')

    def drop_table(self):
        """Drops the intermediate table if it exists."""

        self.database.execute_sql(f'DROP TABLE IF EXISTS {self.path};')

    def write_table(self, filename=None, mode='results', write=True):
        """Writes the selection to a FITS file.

        Parameters
        ----------
        filename : str
            The file to which to write the table. Defaults to
            ``<name>_<plan>.fits``.
        mode : str
            Defines what data to write. If ``'results'``, writes the
            intermediate table (usually just the ``catalogid`` column). If
            ``'targetdb'``, writes all the relevant columns for the targets
            loaded to ``targetdb`` for this carton and plan (must be used after
            `.load` has been called).
        write : bool
            Whether to write the table to disk. If `False`, just returns the
            table object.

        Returns
        -------
        table : `~astropy.table.Table`
            A table object with the selected results.

        """

        if not filename:
            if mode == 'results':
                filename = f'{self.name}_{self.plan}.fits.gz'
            else:
                filename = f'{self.name}_{self.plan}_targetdb.fits.gz'

        log.debug(f'Writing table to {filename}.')

        if mode == 'results':

            results_model = self.get_model()
            assert results_model.table_exists(), \
                'temporary table does not exist.'

            write_query = results_model.select()

            colnames = [field.name for field in write_query._returning]

        elif mode == 'targetdb':

            mag_fields = [field
                          for field in tdb.Magnitude._meta.fields.values()
                          if field.name not in ['pk', 'target_pk', 'target']]

            write_query = (tdb.Target
                           .select(tdb.Target,
                                   *mag_fields,
                                   tdb.CartonToTarget.priority,
                                   tdb.Cadence.label.alias('cadence'))
                           .join(tdb.CartonToTarget)
                           .join(tdb.Magnitude)
                           .join_from(tdb.CartonToTarget, tdb.Cadence,
                                      peewee.JOIN.LEFT_OUTER)
                           .join_from(tdb.CartonToTarget, tdb.Carton)
                           .join(tdb.Version)
                           .where(tdb.Carton.carton == self.name,
                                  tdb.Version.plan == self.plan,
                                  tdb.Version.target_selection >> True))

            colnames = []
            for col in write_query._returning:
                if isinstance(col, peewee.ForeignKeyField):
                    colnames.append(col.column_name)
                elif isinstance(col, peewee.Alias):
                    colnames.append(col._alias)
                else:
                    colnames.append(col.name)

        else:
            raise ValueError('invalud mode. Available modes are '
                             '"results" and "targetdb".')

        if not write_query.exists():
            raise TargetSelectionError('no records found.')

        results = write_query.tuples()
        results = ((col if col is not None else numpy.nan for col in row)
                   for row in tuple(results))

        warnings.filterwarnings(
            'ignore', message='.*converting a masked element to nan.*')

        carton_table = table.Table(rows=results, names=colnames, masked=True)

        if write:
            carton_table.write(filename, overwrite=True)

        return carton_table

    def load(self, overwrite=False):
        """Loads the output of the intermediate table into targetdb."""

        if self.check_targets():
            if overwrite:
                warnings.warn(f'Carton {self.name!r} with plan {self.plan!r} '
                              f'and tag {self.tag!r} already has targets '
                              'loaded. Dropping carton-to-target entries.',
                              TargetSelectionUserWarning)
                self.drop_carton()
            else:
                raise TargetSelectionError(f'Found existing targets for '
                                           f'carton {self.name!r} with plan '
                                           f'{self.plan!r} and tag '
                                           f'{self.tag!r}.')

        RModel = self.get_model()
        if not RModel.table_exists():
            raise TargetSelectionError(f'No temporary table found '
                                       f'{self.full}. Did you call run()?')

        has_targets = (RModel.select()
                       .where(RModel.selected >> True)
                       .exists())

        if not has_targets:
            raise TargetSelectionError('No targets found in '
                                       'intermediate table.')

        with self.database.atomic():
            self.setup_transaction()
            self._create_carton_metadata()
            self._load_targets(RModel)
            self._load_carton_to_target(RModel)
            if self.load_magnitudes:
                self._load_magnitudes(RModel)
            else:
                warnings.warn('Skipping magnitude load.',
                              TargetSelectionUserWarning)

            self.log.debug('Committing records and checking constraints.')

    def check_targets(self):
        """Check if data has been loaded for this carton and targeting plan."""

        has_targets = (tdb.CartonToTarget
                       .select()
                       .join(tdb.Carton)
                       .join(tdb.Version)
                       .where(tdb.Carton.carton == self.name,
                              tdb.Version.plan == self.plan,
                              tdb.Version.target_selection >> True)
                       .exists())

        return has_targets

    def _create_carton_metadata(self):

        mapper_pk = None
        category_pk = None

        # Create targeting plan in tdb.
        version, created = tdb.Version.get_or_create(plan=self.plan,
                                                     tag=self.tag,
                                                     target_selection=True)
        version_pk = version.pk

        if created:
            log.info(f'Created record in targetdb.version for '
                     f'{self.plan!r} with tag {self.tag!r}.')

        if (tdb.Carton.select()
                      .where(tdb.Carton.carton == self.name,
                             tdb.Carton.version_pk == version_pk)
                      .exists()):
            return

        # Create carton and associated values.
        if self.mapper:
            mapper, created_pk = tdb.Mapper.get_or_create(label=self.mapper)
            mapper_pk = mapper.pk
            if created:
                log.debug(f'Created mapper {self.mapper!r}')

        if self.category:
            category, created = tdb.Category.get_or_create(label=self.category)
            category_pk = category.pk
            if created:
                log.debug(f'Created category {self.category!r}')

        tdb.Carton.create(carton=self.name, category_pk=category_pk,
                          program=self.program, mapper_pk=mapper_pk,
                          version_pk=version_pk).save()

        log.debug(f'Created carton {self.name!r}')

    def _load_targets(self, RModel):
        """Load data from the intermediate table tp targetdb.target."""

        log.debug('loading data into targetdb.target.')

        n_inserted = tdb.Target.insert_from(
            cdb.Catalog.select(cdb.Catalog.catalogid,
                               cdb.Catalog.ra,
                               cdb.Catalog.dec,
                               cdb.Catalog.pmra,
                               cdb.Catalog.pmdec,
                               cdb.Catalog.parallax,
                               peewee.Value(EPOCH))
            .join(RModel, on=(cdb.Catalog.catalogid == RModel.catalogid))
            .where(RModel.selected >> True)
            .where(~peewee.fn.EXISTS(
                tdb.Target
                .select(peewee.SQL('1'))
                .where(tdb.Target.catalogid == RModel.catalogid))),
            [tdb.Target.catalogid,
             tdb.Target.ra,
             tdb.Target.dec,
             tdb.Target.pmra,
             tdb.Target.pmdec,
             tdb.Target.parallax,
             tdb.Target.epoch]).returning().execute()

        log.info(f'Inserted {n_inserted:,} new rows into targetdb.target.')

        return

    def _load_magnitudes(self, RModel):
        """Load magnitudes into targetdb.magnitude."""

        log.debug('Loading data into targetdb.magnitude.')

        Magnitude = tdb.Magnitude

        magnitude_paths = self.config['magnitudes']
        fields = [Magnitude.carton_to_target_pk]

        select_from = (RModel
                       .select(tdb.CartonToTarget.pk)
                       .join(tdb.Target,
                             on=(RModel.catalogid == tdb.Target.catalogid))
                       .join(tdb.CartonToTarget)
                       .join(tdb.Carton)
                       .join(tdb.Version)
                       .where(RModel.selected >> True)
                       .where(tdb.Carton.carton == self.name,
                              tdb.Version.plan == self.plan,
                              tdb.Version.tag == self.tag,
                              tdb.Version.target_selection >> True))

        for mag, mpath in magnitude_paths.items():

            fields.append(getattr(Magnitude, mag))
            if hasattr(RModel, mag):
                select_from = select_from.select_extend(getattr(RModel, mag))
                continue

            select_from = select_from.switch(tdb.Target)

            # For each node in the join list we check if the node model has
            # already been join and if so, switch the pointer to that model.
            # Otherwise do a LEFT OUTER join because we want all the rows in
            # the temporary table even if they don't have associated
            # magnitudes. We make sure we only use the "best" match from the
            # cross-match. No need to apply filters on cross-match version_id
            # because catalogid is unique across versions.

            for node in mpath:
                column = None
                if node == mpath[-1]:
                    node, column = node.split('.')
                node_model = self.database.models.get('catalogdb.' + node)
                joins = [model[0] for join in select_from._joins.values()
                         for model in join]
                if node_model in joins:
                    select_from = select_from.switch(node_model)
                else:
                    if node.startswith('catalog_to_'):
                        select_from = select_from.join(
                            node_model,
                            peewee.JOIN.LEFT_OUTER,
                            on=(tdb.Target.catalogid == node_model.catalogid))
                        select_from = (select_from
                                       .where((node_model.best >> True) |
                                              (node_model.catalogid >> None)))
                    else:
                        select_from = select_from.join(node_model,
                                                       peewee.JOIN.LEFT_OUTER)
                if column:
                    select_from = (select_from
                                   .select_extend(getattr(node_model, column)))

        n_inserted = (Magnitude
                      .insert_from(select_from, fields)
                      .returning().execute())

        log.info(f'Inserted {n_inserted:,} new rows into targetdb.magnitude.')

    def _load_carton_to_target(self, RModel):
        """Populate targetdb.carton_to_target."""

        log.debug('Loading data into targetdb.carton_to_target.')

        version_pk = tdb.Version.get(plan=self.plan, tag=self.tag,
                                     target_selection=True)
        carton_pk = tdb.Carton.get(carton=self.name, version_pk=version_pk).pk

        Target = tdb.Target
        CartonToTarget = tdb.CartonToTarget

        select_from = (RModel
                       .select(Target.pk,
                               carton_pk)
                       .join(Target,
                             on=(Target.catalogid == RModel.catalogid))
                       .where(RModel.selected >> True)
                       .where(~peewee.fn.EXISTS(
                           CartonToTarget
                           .select(peewee.SQL('1'))
                           .join(tdb.Carton)
                           .where(CartonToTarget.target_pk == Target.pk,
                                  CartonToTarget.carton_pk == carton_pk,
                                  tdb.Carton.version_pk == version_pk))))

        if self.cadence is not None:

            # Check that not both the carton cadence and the cadence column
            # are not null.
            if RModel.select().where(~(RModel.cadence >> None)).exists():
                raise TargetSelectionError('both carton cadence and target '
                                           'cadence defined. This is not '
                                           'allowed.')

            cadence_pk = tdb.Cadence.get(label=self.cadence)
            select_from = select_from.select_extend(cadence_pk)

        else:

            # If all cadences are null we'll set that as a value and save us
            # a costly join.
            if not RModel.select().where(~(RModel.cadence >> None)).exists():
                select_from = select_from.select_extend(peewee.SQL('null'))
            else:
                select_from = (select_from
                               .select_extend(tdb.Cadence.pk)
                               .switch(RModel)
                               .join(tdb.Cadence, 'LEFT OUTER JOIN',
                                     on=(tdb.Cadence.label == RModel.cadence)))

        if self.priority is None:
            select_from = select_from.select_extend(RModel.priority)
        else:
            select_from = select_from.select_extend(self.priority)

        # Now do the insert
        n_inserted = CartonToTarget.insert_from(
            select_from,
            [CartonToTarget.target_pk,
             CartonToTarget.carton_pk,
             CartonToTarget.cadence_pk,
             CartonToTarget.priority]).returning().execute()

        log.info(f'Inserted {n_inserted:,} rows '
                 'into targetdb.carton_to_target.')

    def drop_carton(self):
        """Drops the entry in ``targetdb.carton``."""

        version = (tdb.Version
                   .select()
                   .where(tdb.Version.plan == self.plan,
                          tdb.Version.tag == self.tag,
                          tdb.Version.target_selection >> True))

        if version.count() == 0:
            return

        tdb.Carton.delete().where(tdb.Carton.carton == self.name,
                                  tdb.Carton.version == version).execute()
