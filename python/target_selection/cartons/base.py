#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-24
# @Filename: base.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import logging

import numpy
import peewee
from astropy import table
from playhouse.migrate import PostgresqlMigrator, migrate

from sdssdb.peewee.sdss5db import SDSS5dbModel, catalogdb, database, targetdb

from .. import config, log, manager


# For now use Gaia but eventually change to Catalog and Catalog.catalogid.
catalog_model = catalogdb.GaiaDR2Source
catalogid_field = catalogdb.GaiaDR2Source.source_id

default_fields = {'ra': catalogdb.GaiaDR2Source.ra,
                  'dec': catalogdb.GaiaDR2Source.dec,
                  'pmra': catalogdb.GaiaDR2Source.pmra,
                  'pmdec': catalogdb.GaiaDR2Source.pmdec,
                  'epoch': catalogdb.GaiaDR2Source.ref_epoch.alias('epoch')}


class CartonMeta(type):

    __require_overload__ = ['build_query']

    def __new__(meta, name, bases, dct):

        if bases:
            assert dct.get('name', None) is not None, 'name is not defined.'

            # Make sure we have overloaded the necessary methods. We do this
            # manually because we are already using one metaclass and cannot
            # use abc.
            for method in meta.__require_overload__:
                assert method in dct, f'overloading of {method!r} is required.'

        return super(CartonMeta, meta).__new__(meta, name, bases, dct)


class BaseCarton(metaclass=CartonMeta):
    """A base class for target cartons.

    This class is not intended for direct instantiation. Instead, it must be
    subclassed and the relevant class attributes overridden with the values
    corresponding to the carton.

    Parameters
    ----------
    targeting_version : str
        The version of the target selection.
    schema : str
        Schema in which the temporary table with the results of the
        query will be created.
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

    orm = 'peewee'

    def __init__(self, targeting_version, schema='sandbox', table_name=None):

        self.targeting_version = targeting_version

        self.database = targetdb.database
        self.schema = 'sandbox'
        self.table_name = table_name or f'temp_{self.name}'

        self.has_run = False

        if self.targeting_version in config:
            self.config = config[self.targeting_version].get(self.name, None)
        else:
            self.config = None

        if not self.orm == 'peewee':
            raise NotImplementedError('not implemented for SQLAlchemy.')

        assert self.database.connected, 'database is not connected.'

        assert self.name, 'carton subclass must override name'
        assert self.category, 'carton subclass must override category'

    def log(self, message, level=logging.INFO):
        """Logs a message with a header of the current target class name."""

        message = f'({self.name}): {message}'
        log.log(level, message)

    def _check_targetdb(self):
        """Checks that target, cadence, and survey are present in targetdb."""

        try:
            program = targetdb.Program.get(label=self.name)
        except peewee.DoesNotExist:
            raise RuntimeError(f'program {self.name!r} does not exist in targetdb.program.')

        if self.survey:
            assert program.survey and program.survey.label == self.survey, \
                f'{self.survey!r} not present in targetdb.survey.'
        else:
            assert program.survey is None, 'targetdb.survey should be empty but is not.'

        # TODO: uncomment this when the cadences are set in the DB.
        # if self.candece:
        #     assert (targetdb.Cadence
        #             .select()
        #             .where(targetdb.Cadence.label == self.cadence)
        #             .count() == 1), f'{self.cadence!r} does not exist in targetdb.cadence.'

        assert program.category and program.category.label == self.category, \
            f'{self.category!r} not present in targetdb.category.'

    def build_query(self):
        """Builds and returns the query.

        The ORM query for the target class. Note that this must be the
        *un-executed* query, which will be executed in `.run`. The select
        statement must be run on ``targetdb.catalog`` and needs to include,
        at least, ``catalogid``. Additional columns such as
        ``ra, dec, pmra, pmdec`` can be specified but are otherwise completed
        when calling `.run`. Magnitude columns can be included and propagated
        to ``targetdb`` but must be aliased as ``magnitude_<band>`` where
        ``<band>`` must be one of the columns in ``targetdb.magnitude``.

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

    def get_model_from_query(self, query=None):
        """Returns a Peewee model with the columns returned by a query."""

        if query is None:
            query = self.build_query()

        class Model(SDSS5dbModel):
            class Meta:
                table_name = self.table_name
                schema = self.schema
                primary_key = False

        returning_fields = query._returning

        for field in returning_fields:

            field_name, resolved_field = self._resolve_field(field)
            is_primary_key = resolved_field.primary_key

            new_field = resolved_field.__class__(primary_key=is_primary_key,
                                                 null=not is_primary_key)

            if is_primary_key:
                Model._meta.set_primary_key(field_name, new_field)
            else:
                Model._meta.add_field(field_name, new_field)

        return Model

    def run(self, tile=None, tile_num=None, progress_bar=True,
            call_post_process=True, **post_process_kawrgs):
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
        call_post_process : bool
            Whether to call the post-process routine.
        post_process_args : dict
            Keyword arguments to be passed to `.post_process`.

        Returns
        -------
        model : :class:`peewee:Model`
            The model for the intermediate table.

        """

        tile = tile or self.tile
        tile_num = tile_num or self.tile_num

        # Check to make sure that the program entries exist in targetdb.
        self._check_targetdb()

        if database.table_exists(self.table_name, schema=self.schema):
            raise RuntimeError(f'temporary table {self.table_name} already exists.')

        self.log('building query ...')
        query = self.build_query()

        # Tweak the query. Start by making sure the catalogid column is selected.
        if catalogid_field not in query._returning:
            raise RuntimeError(f'{catalogid_field} not being returned in query.')

        # Dynamically create a model for the results table.
        ResultsModel = self.get_model_from_query(query)

        # Create sandbox table.
        database.create_tables([ResultsModel])
        self.log(f'created table {self.table_name}', logging.DEBUG)

        self.log(f'running query with tile={tile}', logging.DEBUG)

        if tile is False:

            ResultsModel.insert_from(query, ResultsModel._meta.fields).execute()

        else:

            assert tile_num > 1, 'tile_num must be >= 2'

            if self.tile_region is None:
                ra_space = numpy.linspace(0, 360, num=tile_num)
                dec_space = numpy.linspace(-90, 90, num=tile_num)
            else:
                ra_space = numpy.linspace(self.tile_region[0], self.tile_region[2],
                                          num=tile_num)
                dec_space = numpy.linspace(self.tile_region[1], self.tile_region[3],
                                           num=tile_num)

            n_tiles = (len(ra_space) - 1) * (len(dec_space) - 1)

            if progress_bar:
                counter = manager.counter(total=n_tiles, desc=self.name, unit='ticks')
                counter.update(0)

            nn = 1
            for ii in range(len(ra_space) - 1):
                for jj in range(len(dec_space) - 1):

                    ra0 = ra_space[ii]
                    ra1 = ra_space[ii + 1]
                    dec0 = dec_space[jj]
                    dec1 = dec_space[jj + 1]

                    polygon = (f'(({ra0:.1f}, {dec0:.1f}), ({ra0:.1f}, {dec1:.1f}), '
                               f'({ra1:.1f}, {dec1:.1f}), ({ra1:.1f}, {dec0:.1f}))')

                    if isinstance(query, peewee. ModelSelect):
                        query_model = query.model
                    else:
                        query_model = query._cte_list[0].c

                    tile_query = query.where(
                        peewee.fn.q3c_poly_query(query_model.ra, query_model.dec,
                                                 peewee.SQL(f'\'{polygon}\'::polygon')))

                    self.log(f'tile {nn}/{n_tiles}: {polygon}', logging.DEBUG)

                    ResultsModel.insert_from(tile_query, ResultsModel._meta.fields).execute()

                    nn += 1

                    if progress_bar:
                        counter.update()

        # Add the selected column.
        selected_field = peewee.BooleanField(default=True, null=False)
        migrator = PostgresqlMigrator(database)
        migrate(migrator.set_search_path(self.schema),
                migrator.add_column(self.table_name, 'selected', selected_field))
        ResultsModel._meta.add_field('selected', selected_field)

        if call_post_process:
            mask = self.post_process(ResultsModel, **post_process_kawrgs)
            if mask is True:
                self.log('post-processing returned True. Selecting all records.')
                pass  # No need to do anything
            else:
                assert len(mask) > 0, 'the post-process list is empty'
                self.log(f'selecting {len(mask)} records.')
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

        if self.schema:
            self.database.execute_sql(f'DROP TABLE IF EXISTS {self.schema}.{self.table_name};')
        else:
            self.database.execute_sql(f'DROP TABLE IF EXISTS {self.table_name};')

    def write_table(self, filename=None, model=None):
        """Writes the intermediate table to a FITS file.

        Parameters
        ----------
        filename : str
            The file to which to write the table. Defaults to
            ``<name>_<version>.fits``.
        model : peewee:Model
            The model of the intermediate table. Defaults to use the
            model matching the carton query.

        """

        filename = filename or f'{self.name}_{self.targeting_version}.fits'

        log.debug(f'({self.name}): writing table to {filename}.')

        results_model = model or self.get_model_from_query()

        write_query = results_model.select()

        colnames = [field.name for field in write_query._returning]
        carton_table = table.Table(rows=write_query.tuples(), names=colnames)
        carton_table.write(filename, overwrite=True)
