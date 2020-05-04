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

from sdssdb.peewee.sdss5db import catalogdb, database, targetdb
from sdsstools.color_print import color_text

from target_selection import config, log, manager
from target_selection.exceptions import TargetSelectionError
from target_selection.utils import Timer


class BaseCarton(metaclass=abc.ABCMeta):
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
        query will be created. If `None`, tries to use the ``schema`` parameter
        from the configuration file for this version of target selection. If
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

    def __init__(self, targeting_version, schema=None, table_name=None):

        assert self.name, 'carton subclass must override name'
        assert self.category, 'carton subclass must override category'

        self.version = targeting_version

        if self.version not in config:
            raise TargetSelectionError(f'({self.name}): cannot find version '
                                       f'{self.version!r} in config.')

        self.config = config[self.version].get(self.name, None)

        try:
            self.xmatch_version = config[self.version]['xmatch_version']
        except KeyError:
            raise TargetSelectionError(f'({self.name}): xmatch_version '
                                       'not found in config.')

        self.database = targetdb.database
        assert self.database.connected, 'database is not connected.'

        self.schema = schema or self.config.get('schema', None) or 'sandbox'
        self.table_name = table_name or f'temp_{self.name}'

        if self.cadence:
            assert (targetdb.Cadence
                    .select()
                    .where(targetdb.Cadence.label == self.cadence)
                    .count() == 1), f'{self.cadence!r} does not exist in targetdb.cadence.'

        self.has_run = False
        log.header = f'({self.name}): '

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

    def get_model_from_query(self, query=None):
        """Returns a Peewee model with the columns returned by a query."""

        if query is None:
            version_id = catalogdb.Version.get(version=self.xmatch_version).id
            query = self.build_query(version_id)

        class Model(peewee.Model):

            catalogid = peewee.BigIntegerField(primary_key=True)
            selected = peewee.BooleanField()
            cadence = peewee.TextField(null=True)

            class Meta:
                database = self.database
                table_name = self.table_name
                schema = self.schema
                primary_key = False

        returning_fields = query._returning

        for field in returning_fields:

            field_name, resolved_field = self._resolve_field(field)
            is_primary_key = resolved_field.primary_key

            if field_name in Model._meta.fields:
                continue

            new_field = resolved_field.__class__(primary_key=is_primary_key,
                                                 null=not is_primary_key)

            if is_primary_key:
                Model._meta.set_primary_key(field_name, new_field)
            else:
                Model._meta.add_field(field_name, new_field)

        return Model

    def run(self, tile=None, tile_num=None, progress_bar=True,
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
        post_process_args : dict
            Keyword arguments to be passed to `.post_process`.

        Returns
        -------
        model : :class:`peewee:Model`
            The model for the intermediate table.

        """

        tile = tile or self.tile
        tile_num = tile_num or self.tile_num

        if database.table_exists(self.table_name, schema=self.schema):
            raise RuntimeError(f'temporary table {self.table_name} already exists.')

        log.info('building query ...')
        version_id = catalogdb.Version.get(version=self.xmatch_version).id
        query = self.build_query(version_id)

        # Make sure the catalogid column is selected.
        if catalogdb.Catalog.catalogid not in query._returning:
            raise RuntimeError(f'catalogid is not being returned in query.')

        # Dynamically create a model for the results table.
        ResultsModel = self.get_model_from_query(query)

        # Create sandbox table.
        database.create_tables([ResultsModel])
        log.debug(f'created table {self.table_name!r}')

        # Make the "selected" column detault to true. We cannot do this in
        # Peewee because field defaults are implemented only on the Python
        # side.
        with self.database.atomic():
            self.database.execute_sql(
                f'ALTER TABLE {self.schema}.{self.table_name} '
                'ALTER COLUMN selected SET DEFAULT true;')

        log.debug(f'running query with tile={tile}')

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
                counter = manager.counter(total=n_tiles, desc=self.name, unit='ticks')
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

        log.info(f'inserted {n_rows:,} rows into {self.table_name} '
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

        filename = filename or f'{self.name}_{self.version}.fits'

        log.debug(f'({self.name}): writing table to {filename}.')

        results_model = model or self.get_model_from_query()

        write_query = results_model.select()

        colnames = [field.name for field in write_query._returning]
        carton_table = table.Table(rows=write_query.tuples(), names=colnames)
        carton_table.write(filename, overwrite=True)
