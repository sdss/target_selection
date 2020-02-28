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

from sdssdb.peewee.sdss5db import SDSS5dbModel, catalogdb, database, targetdb

from .. import config, log


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

    Parameters
    ----------
    targeting_version : str
        The version of the target selection.

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

    """

    name = None
    cadence = None
    category = None
    survey = None

    orm = 'peewee'

    def __init__(self, targeting_version):

        self.targeting_version = targeting_version
        self.database = targetdb.database

        if self.targeting_version in config:
            self.config = config[self.targeting_version].get(self.name, None)
        else:
            self.config = None

        if not self.orm == 'peewee':
            raise NotImplementedError('not implemented for SQLAlchemy.')

        assert self.database.connected, 'database is not connected.'

    def log(self, message, level=logging.INFO):
        """Logs a message with a header of the current target class name."""

        message = f'({self.name}): {message}'
        log.log(level, message)

    def _check_targetdb(self):
        """Checks that target, cadence, and survey are present in targetdb."""

        try:
            program = targetdb.Program.get(label=self.name)
        except peewee.DoesNotExist:
            raise f'program {self.name} does not exist in targetdb.program.'

        if self.survey:
            assert program.survey and program.survey.label == self.survey, \
                f'{self.survey!r} not present in targetdb.survey.'
        else:
            assert program.survey is None, 'targetdb.survey should be empty but is not.'

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
            A `~peewee.Select` or `~peewee.ModelSelect` query.

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

    def _get_model_from_query(self, table_name, query):
        """Returns a Peewee model with the columns returned by a query."""

        class Model(SDSS5dbModel):
            class Meta:
                schema = 'catalogdb'
                primary_key = False

        Model._meta.table_name = table_name
        returning_fields = query._returning

        for field in returning_fields:

            field_name, resolved_field = self._resolve_field(field)
            is_primary_key = resolved_field.primary_key

            new_field = resolved_field.__class__(primary_key=is_primary_key,
                                                 null=not is_primary_key)

            Model._meta.add_field(field_name, new_field)

        return Model

    def run(self, schema='sandbox', table_name=None, tile=False, tile_num=21):
        """Executes the query and stores the results.

        Parameters
        ----------
        schema : str
            Schema in which the temporary table with the results of the
            query will be created.
        table_name : str
            The name of the temporary table. Defaults to ``temp_<name>`` where
            ``<name>`` is the target class name.
        tile : bool
            Whether to perform multiple queries by tiling the sky with a
            rectangular grid instead of running a single query. If `True`,
            a filter using ``q3c_poly_query`` is added to the query.
        tile_num : int
            The number of tile nodes in which to divide the RA and Dec axes
            when tiling.

        """

        self.log(f'running target selection for target class {self.name!r}.')

        # Check to make sure that the program entries exist in targetdb.
        self._check_targetdb()

        table_name = table_name or f'temp_{self.name}'

        if database.table_exists(table_name, schema=schema):
            raise RuntimeError(f'temporary table {table_name} already exists.')

        query = self.build_query()

        # Tweak the query. Start by making sure the catalogid column is selected.
        if catalogid_field not in query._returning:
            raise RuntimeError(f'{catalogid_field} not being returned in query.')

        # Dynamically create a model for the results table.
        ResultsModel = self._get_model_from_query(table_name, query)
        ResultsModel._meta.schema = schema

        # Create sandbox table.
        database.create_tables([ResultsModel])
        self.log(f'created table {table_name}', logging.DEBUG)

        self.log(f'running query with tile={tile}', logging.DEBUG)

        if tile is False:

            ResultsModel.insert_from(query, ResultsModel._meta.fields).execute()

        else:

            assert tile_num > 1, 'tile_num must be >= 2'

            ra_space = numpy.linspace(0, 20, num=tile_num)
            dec_space = numpy.linspace(0, 20, num=tile_num)
            n_tiles = (len(ra_space) - 1) * (len(dec_space) - 1)

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

                    self.log(f'executing tile {nn}/{n_tiles}: {polygon}', logging.DEBUG)

                    ResultsModel.insert_from(tile_query, ResultsModel._meta.fields).execute()

                    nn += 1

        return ResultsModel
