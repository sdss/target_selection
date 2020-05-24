#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-24
# @Filename: base.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import abc

import numpy
import peewee
from astropy import table

from sdssdb.peewee import BaseModel
from sdssdb.peewee.sdss5db import catalogdb as cdb
from sdssdb.peewee.sdss5db import targetdb as tdb
from sdsstools.color_print import color_text

from target_selection import __version__, config, log, manager
from target_selection.exceptions import TargetSelectionError
from target_selection.utils import Timer


class BaseCarton(metaclass=abc.ABCMeta):
    """A base class for target cartons.

    This class is not intended for direct instantiation. Instead, it must be
    subclassed and the relevant class attributes overridden with the values
    corresponding to the carton.

    Parameters
    ----------
    targeting_plan : str
        The target selection plan version.
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
    survey : str
        The survey associated with this carton.
    orm : str
        The ORM library to be used, ``peewee`` or ``sqlalchemy``.
    tag : str
        The version of the ``target_selection`` code used.
    tile : bool
        Whether to tile the query instead of running it all at once.
    tile_region : tuple
        A tuple defining the region over which to tile with the format
        ``(ra0, dec0, ra1, dec1)``. This setting is only applied if
        ``tile=True``; for a non-tiled delimitation of the query region, use
        the `q3c <https://github.com/segasai/q3c>`__ function
        ``q3c_poly_query`` when building the query.
    tile_num : int
        The number of tile nodes in which to divide the RA and Dec axes
        when tiling.

    """

    name = None
    cadence = None
    category = None
    survey = None

    tile = False
    tile_region = None
    tile_num = None

    def __init__(self, targeting_plan, schema=None, table_name=None):

        assert self.name, 'carton subclass must override name'
        assert self.category, 'carton subclass must override category'

        self.plan = targeting_plan
        self.tag = __version__

        if self.plan not in config:
            raise TargetSelectionError(f'({self.name}): cannot find plan '
                                       f'{self.plan!r} in config.')

        self.config = config[self.plan]

        if 'parameters' in self.config:
            self.parameters = self.config['parameters'].get(self.name, None)
        else:
            self.parameters = None

        try:
            self.xmatch_plan = self.config['xmatch_plan']
        except KeyError:
            raise TargetSelectionError(f'({self.name}): xmatch_plan '
                                       'not found in config.')

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

        self.has_run = False

    @property
    def path(self):
        """The schema-qualified path to the output table."""

        return self.schema + '.' + self.table_name

    @abc.abstractmethod
    def build_query(self, version_id):
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

    def _resolve_field(self, field, column=None):
        """Returns the field associated with a field, column, or alias."""

        if isinstance(field, peewee.Field):
            return (field.name, field)
        elif isinstance(field, peewee.Alias):
            return (field._alias, self._resolve_field(field.alias())[1])
        elif isinstance(field, peewee.CTE):
            for returning in field._query._returning:
                if returning.name == column.name:
                    return self._resolve_field(returning, column=column)
        elif isinstance(field, peewee.Column):
            return self._resolve_field(field.source, column=field)

        raise ValueError(f'failed resolving column {field!r}')

    def get_model(self, query=None, use_reflection=False):
        """Returns a Peewee model with the columns returned by a query."""

        if query is None:
            version_id = cdb.Version.get(plan=self.xmatch_plan).id
            query = self.build_query(version_id)

        class Model(BaseModel):

            catalogid = peewee.BigIntegerField(primary_key=True)
            selected = peewee.BooleanField()
            cadence = peewee.TextField(null=True)

            class Meta:
                database = self.database
                table_name = self.table_name
                schema = self.schema
                reflection_options = {'skip_foreign_keys': True}

        if use_reflection:
            Model.reflect()
            return Model

        returning_fields = query._returning

        for field in returning_fields:

            field_name, resolved_field = self._resolve_field(field)
            is_primary_key = resolved_field.primary_key

            if field_name in Model._meta.fields:
                continue

            new_field = resolved_field.__class__(primary_key=is_primary_key,
                                                 null=not is_primary_key)

            Model._meta.add_field(field_name, new_field)

        return Model

    def run(self, tile=None, tile_num=None, progress_bar=True, overwrite=False,
            **post_process_kawrgs):
        """Executes the query and post-process steps, and stores the results.

        This method calls `.build_query` and runs the returned query. The
        output of the query is stored in a temporary table whose schema and
        table name are defined when the object is instantiated. The target
        selection query can be run as a single query or tiled.

        After the query has run, the `.post_process` routine is called if the
        method has been overridden for the given carton. The returned mask
        is stored as a new column, ``selected``, in the intermediary table.

        Parameters
        ----------
        tile : bool
            Whether to perform multiple queries by tiling the sky with a
            rectangular grid instead of running a single query. If `True`,
            a filter using ``q3c_poly_query`` is added to the query.
        tile_num : int
            The number of tile nodes in which to divide the RA and Dec axes
            when tiling.
        progress_bar : bool
            Use a progress bar when tiling.
        overwrite : bool
            Whether to overwrite the intermediary table if already exists.
        post_process_args : dict
            Keyword arguments to be passed to `.post_process`.

        Returns
        -------
        model : :class:`peewee:Model`
            The model for the intermediate table.

        """

        tile = tile or self.tile
        tile_num = tile_num or self.tile_num

        if self.database.table_exists(self.table_name, schema=self.schema):
            if overwrite:
                log.info(f'dropping table {self.path!r}.')
                self.drop_table()
            else:
                raise RuntimeError(f'temporary table {self.path!r} '
                                   'already exists.')

        log.debug('building query ...')
        version_id = cdb.Version.get(plan=self.xmatch_plan).id
        query = self.build_query(version_id)

        # Make sure the catalogid column is selected.
        if cdb.Catalog.catalogid not in query._returning:
            raise RuntimeError(f'catalogid is not being returned in query.')

        # Dynamically create a model for the results table.
        ResultsModel = self.get_model(query)

        # Create sandbox table.
        self.database.create_tables([ResultsModel])
        log.debug(f'created table {self.table_name!r}')

        # Make the "selected" column detault to true. We cannot do this in
        # Peewee because field defaults are implemented only on the Python
        # side.
        with self.database.atomic():
            self.database.execute_sql(
                f'ALTER TABLE {self.schema}.{self.table_name} '
                'ALTER COLUMN selected SET DEFAULT true;')

        log.info(f'running query with tile={tile}')

        insert_fields = [ResultsModel._meta.fields[field]
                         for field in ResultsModel._meta.fields
                         if field not in ['selected', 'cadence']]

        if tile is False:

            log.debug(f'INSERT INTO {ResultsModel._meta.table_name}: ' +
                      color_text(str(query), 'darkgrey'))

            with self.database.atomic():
                with Timer() as timer:
                    n_rows = (ResultsModel
                              .insert_from(query, insert_fields)
                              .returning()
                              .execute())

        else:

            assert tile_num > 1, 'tile_num must be >= 2'

            if self.tile_region is None:
                ra_space = numpy.linspace(0, 360, num=tile_num)
                dec_space = numpy.linspace(-90, 90, num=tile_num)
            else:
                ra_space = numpy.linspace(self.tile_region[0],
                                          self.tile_region[2],
                                          num=tile_num)
                dec_space = numpy.linspace(self.tile_region[1],
                                           self.tile_region[3],
                                           num=tile_num)

            n_tiles = (len(ra_space) - 1) * (len(dec_space) - 1)

            if progress_bar:
                counter = manager.counter(total=n_tiles,
                                          desc=self.name,
                                          unit='ticks')
                counter.update(0)

            with self.database.atomic():
                with Timer() as timer:

                    nn = 1
                    for ii in range(len(ra_space) - 1):
                        for jj in range(len(dec_space) - 1):

                            ra0 = ra_space[ii]
                            ra1 = ra_space[ii + 1]
                            dec0 = dec_space[jj]
                            dec1 = dec_space[jj + 1]

                            polygon = (f'(({ra0:.1f}, {dec0:.1f}), '
                                       f'({ra0:.1f}, {dec1:.1f}), '
                                       f'({ra1:.1f}, {dec1:.1f}), '
                                       f'({ra1:.1f}, {dec0:.1f}))')

                            if isinstance(query, peewee. ModelSelect):
                                query_model = query.model
                            else:
                                query_model = query._cte_list[0].c

                            tile_query = query.where(
                                peewee.fn.q3c_poly_query(
                                    query_model.ra,
                                    query_model.dec,
                                    peewee.SQL(f'\'{polygon}\'::polygon')))

                            log.debug(f'tile {nn}/{n_tiles}: {polygon}')

                            n_rows = (ResultsModel
                                      .insert_from(tile_query, insert_fields)
                                      .returning()
                                      .execute())

                            nn += 1

                            if progress_bar:
                                counter.update()

        log.info(f'inserted {n_rows:,} rows into {self.path!r} '
                 f'in {timer.interval:.3f} s.')

        mask = self.post_process(ResultsModel, **post_process_kawrgs)
        if mask is True:
            log.debug('post-processing returned True. Selecting all records.')
            pass  # No need to do anything
        else:
            assert len(mask) > 0, 'the post-process list is empty'
            log.debug(f'selecting {len(mask)} records.')
            (ResultsModel
             .update({ResultsModel.selected: False})
             .where(ResultsModel._meta.primary_key.not_in(mask))
             .execute())

        self.has_run = True

        return ResultsModel

    def post_process(self, model, **kwargs):
        """Post-processes the temporary table.

        This method provides a framework for applying non-SQL operations on
        carton query. It receives the model for the temporary table and can
        perform any operation as long as it returns a tuple with the list of
        catalogids that should be selected.

        This method can also be used to set the ``cadence`` column in the
        temporary table. This column will be used to set the target cadence if
        the carton `.cadence` attribute is not set.

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

    def drop_table(self):
        """Drops the intermediate table if it exists."""

        self.database.execute_sql(f'DROP TABLE IF EXISTS {self.path};')

    def write_table(self, filename=None, model=None):
        """Writes the intermediate table to a FITS file.

        Parameters
        ----------
        filename : str
            The file to which to write the table. Defaults to
            ``<name>_<plan>.fits``.
        model : peewee:Model
            The model of the intermediate table. Defaults to use the
            model matching the carton query.

        """

        filename = filename or f'{self.name}_{self.plan}.fits'

        log.debug(f'writing table to {filename}.')

        results_model = model or self.get_model()

        write_query = results_model.select()

        colnames = [field.name for field in write_query._returning]
        carton_table = table.Table(rows=write_query.tuples(), names=colnames)
        carton_table.write(filename, overwrite=True)

    def load(self):
        """Loads the output of the intermediate table into targetdb."""

        if self.check_targets():
            raise TargetSelectionError(f'found existing targets for carton '
                                       f'{self.name!r} with plan '
                                       f'{self.plan!r}.')

        if not self.has_run:
            raise TargetSelectionError('the query needs to be run before '
                                       'calling load().')

        RModel = self.get_model(use_reflection=True)

        has_targets = (RModel.select()
                       .where(RModel.selected == True)  # noqa
                       .exists())

        if not has_targets:
            raise TargetSelectionError('no targets found in '
                                       'intermediate table.')

        with self.database.atomic():

            self._create_program_metadata()
            self._load_data(RModel)
            self._load_magnitudes(RModel)
            self._load_program_to_target(RModel)

    def check_targets(self):
        """Check if data has been loaded for this carton and targeting plan."""

        has_targets = (tdb.ProgramToTarget
                       .select()
                       .join(tdb.Program)
                       .join(tdb.Plan)
                       .where(tdb.Program.label == self.name,
                              tdb.Plan.label == self.plan)
                       .exists())

        return has_targets

    def _create_program_metadata(self):

        survey_pk = None
        category_pk = None

        # Create targeting plan in tdb.
        plan_pk, created = tdb.Plan.get_or_create(label=self.plan,
                                                  tag=self.tag,
                                                  target_selection=True)
        if created:
            log.info(f'created record in targetdb.plan for {self.plan!r}.')

        plan_pk = tdb.Plan.get(label=self.plan, target_selection=True)

        if (tdb.Program.select()
                       .where(tdb.Program.label == self.name,
                              tdb.Program.plan_pk == plan_pk)
                       .exists()):
            return

        # Create program and associated values.
        if self.survey:
            survey_pk, created_pk = tdb.Survey.get_or_create(label=self.name)
            if created:
                log.debug(f'created survey {self.survey!r}')

        if self.category:
            category_pk, created = tdb.Category.get_or_create(
                label=self.category)
            if created:
                log.debug(f'created category {self.category!r}')

        tdb.Program.create(label=self.name, category_pk=category_pk,
                           survey_pk=survey_pk, plan_pk=plan_pk)
        log.debug(f'created program {self.name!r}')

    def _load_data(self, RModel):
        """Load data from the intermediate table tp targetdb.target."""

        log.debug('loading data into targetdb.target.')

        n_inserted = tdb.Target.insert_from(
            cdb.Catalog.select(cdb.Catalog.catalogid,
                               cdb.Catalog.ra,
                               cdb.Catalog.dec,
                               cdb.Catalog.pmra,
                               cdb.Catalog.pmdec,
                               cdb.Catalog.parallax)
            .join(RModel, on=(cdb.Catalog.catalogid == RModel.catalogid))
            .where(RModel.selected == True)  # noqa
            .where(~peewee.fn.EXISTS(
                tdb.Target
                .select(peewee.SQL('1'))
                .where(tdb.Target.catalogid == RModel.catalogid))),
            [tdb.Target.catalogid,
             tdb.Target.ra,
             tdb.Target.dec,
             tdb.Target.pmra,
             tdb.Target.pmdec,
             tdb.Target.parallax]).returning().execute()

        log.info(f'inserted {n_inserted:,} rows into targetdb.target.')

        return

    def _load_magnitudes(self, RModel):
        """Load magnitudes into targetdb.magnitude."""

        log.debug('loading data into targetdb.magnitude.')

        Magnitude = tdb.Magnitude

        magnitude_paths = self.config['magnitudes']
        fields = [Magnitude.target_pk]

        select_from = (tdb.Target
                       .select(tdb.Target.pk)
                       .join(RModel,
                             on=(RModel.catalogid == tdb.Target.catalogid))
                       .join(cdb.Catalog,
                             on=(RModel.catalogid == cdb.Catalog.catalogid)))

        for mag, mpath in magnitude_paths.items():

            fields.append(getattr(Magnitude, mag))

            select_from = select_from.switch(cdb.Catalog)

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
                    select_from = select_from.join(node_model,
                                                   peewee.JOIN.LEFT_OUTER)
                if node.startswith('catalog_to_'):
                    select_from = (select_from
                                   .where((node_model.best >> True) |
                                          (node_model.catalogid >> None)))
                if column:
                    select_from = (select_from
                                   .select_extend(getattr(node_model, column)))

        n_inserted = (Magnitude
                      .insert_from(select_from, fields)
                      .returning().execute())

        log.debug(f'inserted {n_inserted:,} rows into targetdb.magnitude.')

    def _load_program_to_target(self, RModel):
        """Populate targetdb.program_to_target."""

        log.debug('loading data into targetdb.program_to_target.')

        plan_pk = tdb.Plan.get(label=self.plan, target_selection=True)
        program_pk = tdb.Program.get(label=self.name, plan_pk=plan_pk).pk

        Target = tdb.Target
        ProgramToTarget = tdb.ProgramToTarget

        select_from = (RModel
                       .select(Target.pk,
                               program_pk)
                       .join(Target,
                             on=(Target.catalogid == RModel.catalogid))
                       .where(~peewee.fn.EXISTS(
                           ProgramToTarget
                           .select(peewee.SQL('1'))
                           .join(tdb.Program)
                           .where(ProgramToTarget.target_pk == Target.pk,
                                  ProgramToTarget.program_pk == program_pk,
                                  tdb.Program.plan_pk == plan_pk))))

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

        # Now do the insert
        n_inserted = ProgramToTarget.insert_from(
            select_from,
            [ProgramToTarget.target_pk,
             ProgramToTarget.program_pk,
             ProgramToTarget.cadence_pk]).returning().execute()

        log.debug(f'inserted {n_inserted:,} rows into '
                  'targetdb.program_to_target.')
