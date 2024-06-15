#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2021-04-29
# @Filename: tools.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os

import numpy
import peewee
from astropy.table import Table

from sdssdb.utils.ingest import copy_data, create_model_from_table

from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError
from target_selection.utils import vacuum_table


def get_file_carton(filename):
    """Returns a carton class that creates a carton based on a FITS file.
    The FITS file is located in the open_fiber_path which is specified in
    python/config/target_selection.yml.
    The list of FITS files to be loaded is specified in the
    file open_fiber_file_list.txt which is in the directory open_fiber_path.
    """

    # Import this here to prevent this module not being importable if the database
    # connection is not ready.
    from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                                 CatalogToLegacy_Survey_DR8,
                                                 CatalogToLegacy_Survey_DR10,
                                                 CatalogToPanstarrs1,
                                                 CatalogToTIC_v8,
                                                 CatalogToTwoMassPSC, Gaia_DR2,
                                                 Gaia_DR3, Legacy_Survey_DR8,
                                                 Legacy_Survey_DR10,
                                                 Panstarrs1, TIC_v8,
                                                 TwoMassPSC)

    class FileCarton(BaseCarton):

        can_offset = None  # Will be set in query.

        def __init__(self, targeting_plan, config_file=None, schema=None, table_name=None):

            self._file_path = filename

            self._table = Table.read(self._file_path)
            if self._table.masked:
                self._table = self._table.filled()

            unique_cartonname = numpy.unique(self._table['cartonname'])
            if len(unique_cartonname) == 1:
                self.name = unique_cartonname[0].lower()
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' contains more than one cartonname')

            unique_can_offset = numpy.unique(self._table['can_offset'])
            if len(unique_can_offset) > 1:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' contains more than one' +
                                           ' value of can_offset:' +
                                           ' can_offset values must be ' +
                                           ' all 0 or all 1')

            if (unique_can_offset[0] != 1) and (unique_can_offset[0] != 0):
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' can_offset can only be 0 or 1.' +
                                           ' can_offset is ' +
                                           str(unique_can_offset[0]))

            unique_inertial = numpy.unique(self._table['inertial'])
            if len(unique_inertial) > 2:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' contains more than two' +
                                           ' values of inertial:' +
                                           ' inertial values must be ' +
                                           ' 0 or 1')

            if (unique_inertial[0] != 1) and (unique_inertial[0] != 0):
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' inertial can only be 0 or 1.' +
                                           ' inertial is ' +
                                           str(unique_inertial[0]))

            # If there is only one inertial value then the above statement
            # is enough. Otherwise, we need to run the below check.
            if len(unique_inertial) == 2:
                if (unique_inertial[1] != 1) and (unique_inertial[1] != 0):
                    raise TargetSelectionError('error in get_file_carton(): ' +
                                               filename +
                                               ' inertial can only be 0 or 1.' +
                                               ' inertial is ' +
                                               str(unique_inertial[1]))

            # The valid_program list is from the output of the below command.
            # select distinct(program) from targetdb.carton order by program;
            #
            # mwm_bin is for future mwm binary star cartons
            valid_program = [
                'bhm_aqmes',
                'bhm_csc',
                'bhm_filler',
                'bhm_rm',
                'bhm_spiders',
                'commissioning',
                'mwm_bin',
                'mwm_cb',
                'mwm_dust',
                'mwm_erosita',
                'mwm_filler',
                'mwm_galactic',
                'mwm_gg',
                'mwm_halo',
                'mwm_legacy',
                'mwm_magcloud',
                'mwm_ob',
                'mwm_planet',
                'mwm_rv',
                'mwm_snc',
                'mwm_tessob',
                'mwm_tessrgb',
                'mwm_validation',
                'mwm_wd',
                'mwm_yso',
                'open_fiber',
                'ops',
                'ops_sky',
                'ops_std',
                'SKY']

            # The valid_category list is from CartonImportTable.pdf
            valid_category = ['science', 'standard_apogee',
                              'standard_boss', 'guide',
                              'sky_boss', 'sky_apogee', 'standard',
                              'sky', 'veto location boss',
                              'veto_location_apogee']

            # The valid_category list is from CartonImportTable.pdf
            valid_mapper = ['', 'MWM', 'BHM']

            unique_category = numpy.unique(self._table['category'])
            if len(unique_category) == 1:
                self.category = unique_category[0].lower()
                if (self.category not in valid_category):
                    raise TargetSelectionError('error in get_file_carton(): ' +
                                               filename +
                                               ' contains invalid category = ' +
                                               self.category)
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' contains more than one category')

            unique_program = numpy.unique(self._table['program'])
            if len(unique_program) == 1:
                self.program = unique_program[0].lower()
                if (self.program not in valid_program):
                    raise TargetSelectionError('error in get_file_carton(): ' +
                                               filename +
                                               ' contains invalid program = ' +
                                               self.program)
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' contains more than one program')

            unique_mapper = numpy.unique(self._table['mapper'])
            if len(unique_mapper) == 1:
                # We do not use lower() for mapper because
                # allowed values for mapper are '' or 'MWM' or 'BHM'.
                self.mapper = unique_mapper[0]
                if (self.mapper not in valid_mapper):
                    raise TargetSelectionError('error in get_file_carton(): ' +
                                               filename +
                                               ' contains invalid mapper = ' +
                                               self.mapper)
                if (self.mapper == ''):
                    self.mapper = None
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           ' contains more than one mapper')

            basename_fits = os.path.basename(filename)
            basename_parts = os.path.splitext(basename_fits)
            basename = basename_parts[0]
            carton_name_from_filename = basename.lower()

            if (self.name != carton_name_from_filename):
                raise TargetSelectionError('filename parameter of get_file_carton() and ' +
                                           'cartonname in FITS file do not match.' + '\n' +
                                           'carton_name_from_filename = ' +
                                           carton_name_from_filename +
                                           ' cartonname = ' + self.name)

            super().__init__(
                targeting_plan,
                config_file=config_file,
                schema=schema,
                table_name=table_name,
            )

            self._disable_query_log = True

        def build_query(self, version_id, query_region=None):

            self.log.debug(f'Processing file {self._file_path}.')

            # We need to copy the data to a temporary table so that we can
            # join on it. We could use a Peewee ValueList but for large tables
            # that will hit the limit of 1GB in PSQL.

            # Create model for temporary table from FITS table columns.
            # This works fine because we know there are no arrays.
            temp_table = self.name.lower() + '_temp'
            temp = create_model_from_table(temp_table, self._table)
            temp._meta.database = self.database
            temp.create_table(temporary=True)

            # Copy data.
            copy_data(self._table, self.database, temp_table)

            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("Gaia_DR3_Source_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("Gaia_DR2_Source_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("LegacySurvey_DR8_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("LegacySurvey_DR10_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("PanSTARRS_DR2_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("TwoMASS_ID")')
            vacuum_table(self.database, temp_table, vacuum=False, analyze=True)

            inertial_case = peewee.Case(
                None,
                ((temp.inertial.cast('boolean').is_null(), False),),
                temp.inertial.cast('boolean'))

            # Below we make aliases for the temp table column names.
            # For example, column temp.Gaia_DR3_Source_ID has alias 'gaia_source_id'.
            # However, we do not use the column name aliases later.
            # i.e. later we use the full column name temp.Gaia_DR3_Source_ID.
            query_common = (Catalog
                            .select(Catalog.catalogid,
                                    temp.Gaia_DR3_Source_ID.alias('gaia_dr3_source_id'),
                                    temp.Gaia_DR2_Source_ID.alias('gaia_source_id'),
                                    temp.LegacySurvey_DR8_ID.alias('ls_id8'),
                                    temp.LegacySurvey_DR10_ID.alias('ls_id10'),
                                    temp.PanSTARRS_DR2_ID.alias('catid_objid'),
                                    temp.TwoMASS_ID.alias('designation'),
                                    Catalog.ra,
                                    Catalog.dec,
                                    temp.delta_ra.cast('double precision'),
                                    temp.delta_dec.cast('double precision'),
                                    inertial_case.alias('inertial'),
                                    temp.cadence,
                                    temp.priority,
                                    temp.instrument,
                                    temp.can_offset.cast('boolean').alias('can_offset'),
                                    peewee.Value(0).alias('value'))
                            .distinct(Catalog.catalogid))

            query_gaia_dr3 = \
                (query_common
                 .join(CatalogToGaia_DR3)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(temp,
                       on=(temp.Gaia_DR3_Source_ID == Gaia_DR3.source_id))
                 .switch(Catalog)
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        (CatalogToGaia_DR3.best >> True) |
                        CatalogToGaia_DR3.best.is_null(),
                        Catalog.version_id == version_id))

            query_gaia_dr2 = \
                (query_common
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .join(temp,
                       on=(temp.Gaia_DR2_Source_ID == Gaia_DR2.source_id))
                 .switch(Catalog)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        (CatalogToTIC_v8.best >> True) |
                        CatalogToTIC_v8.best.is_null(),
                        Catalog.version_id == version_id))

            query_legacysurvey_dr8 = \
                (query_common
                 .join(CatalogToLegacy_Survey_DR8)
                 .join(Legacy_Survey_DR8)
                 .join(temp,
                       on=(temp.LegacySurvey_DR8_ID == Legacy_Survey_DR8.ls_id))
                 .switch(Catalog)
                 .where(CatalogToLegacy_Survey_DR8.version_id == version_id,
                        (CatalogToLegacy_Survey_DR8.best >> True) |
                        CatalogToLegacy_Survey_DR8.best.is_null(),
                        Catalog.version_id == version_id))

            query_legacysurvey_dr10 = \
                (query_common
                 .join(CatalogToLegacy_Survey_DR10)
                 .join(Legacy_Survey_DR10)
                 .join(temp,
                       on=(temp.LegacySurvey_DR10_ID == Legacy_Survey_DR10.ls_id))
                 .switch(Catalog)
                 .where(CatalogToLegacy_Survey_DR10.version_id == version_id,
                        (CatalogToLegacy_Survey_DR10.best >> True) |
                        CatalogToLegacy_Survey_DR10.best.is_null(),
                        Catalog.version_id == version_id))

            query_panstarrs_dr2 = \
                (query_common
                 .join(CatalogToPanstarrs1)
                 .join(Panstarrs1)
                 .join(temp,
                       on=(temp.PanSTARRS_DR2_ID == Panstarrs1.catid_objid))
                 .switch(Catalog)
                 .where(CatalogToPanstarrs1.version_id == version_id,
                        (CatalogToPanstarrs1.best >> True) |
                        CatalogToPanstarrs1.best.is_null(),
                        Catalog.version_id == version_id))

            query_twomass_psc = \
                (query_common
                 .join(CatalogToTwoMassPSC)
                 .join(TwoMassPSC)
                 .join(temp,
                       on=(temp.TwoMASS_ID == TwoMassPSC.designation))
                 .switch(Catalog)
                 .where(CatalogToTwoMassPSC.version_id == version_id,
                        (CatalogToTwoMassPSC.best >> True) |
                        CatalogToTwoMassPSC.best.is_null(),
                        Catalog.version_id == version_id))

            len_table = len(self._table)

            len_gaia_dr3 =\
                len(self._table[self._table['Gaia_DR3_Source_ID'] > 0])

            len_gaia_dr2 =\
                len(self._table[self._table['Gaia_DR2_Source_ID'] > 0])

            len_legacysurvey_dr8 =\
                len(self._table[self._table['LegacySurvey_DR8_ID'] > 0])

            len_legacysurvey_dr10 =\
                len(self._table[self._table['LegacySurvey_DR10_ID'] > 0])

            len_panstarrs_dr2 =\
                len(self._table[self._table['PanSTARRS_DR2_ID'] > 0])

            # TwoMass_ID corresponds to the designation column of
            # the table catalogdb.twomass_psc.
            # Since the designation column is a text column, below
            # we are comparing it to the string 'NA' and not the integer 0.
            #
            len_twomass_psc =\
                len(self._table[self._table['TwoMASS_ID'] != 'NA'])

            # There must be exactly one non-zero id per row else raise an exception.
            if ((len_gaia_dr3 + len_gaia_dr2 +
                 len_legacysurvey_dr8 + len_legacysurvey_dr10 +
                 len_panstarrs_dr2 + len_twomass_psc) != len_table):
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           '(len_gaia_dr3 + len_gaia_dr2 + ' +
                                           'len_legacysurvey_dr8 + len_legacysurvey_dr10 +' +
                                           'len_panstarrs_dr2 + len_twomass_psc) != ' +
                                           'len_table')

            if (len_gaia_dr3 > 0):
                is_gaia_dr3 = True
            else:
                is_gaia_dr3 = False

            if (len_gaia_dr2 > 0):
                is_gaia_dr2 = True
            else:
                is_gaia_dr2 = False

            if (len_legacysurvey_dr8 > 0):
                is_legacysurvey_dr8 = True
            else:
                is_legacysurvey_dr8 = False

            if (len_legacysurvey_dr10 > 0):
                is_legacysurvey_dr10 = True
            else:
                is_legacysurvey_dr10 = False

            if (len_panstarrs_dr2 > 0):
                is_panstarrs_dr2 = True
            else:
                is_panstarrs_dr2 = False

            if (len_twomass_psc > 0):
                is_twomass_psc = True
            else:
                is_twomass_psc = False

            query = None

            if (is_gaia_dr3 is True):
                if (query is None):
                    query = query_gaia_dr3
                else:
                    query = query | query_gaia_dr3

            if (is_gaia_dr2 is True):
                if (query is None):
                    query = query_gaia_dr2
                else:
                    query = query | query_gaia_dr2

            if (is_legacysurvey_dr8 is True):
                if (query is None):
                    query = query_legacysurvey_dr8
                else:
                    query = query | query_legacysurvey_dr8

            if (is_legacysurvey_dr10 is True):
                if (query is None):
                    query = query_legacysurvey_dr10
                else:
                    query = query | query_legacysurvey_dr10

            if (is_panstarrs_dr2 is True):
                if (query is None):
                    query = query_panstarrs_dr2
                else:
                    query = query | query_panstarrs_dr2

            if (is_twomass_psc is True):
                if (query is None):
                    query = query_twomass_psc
                else:
                    query = query | query_twomass_psc

            if (query is None):
                # At least one of the four boolean variables above
                # must be True, so we should not get here.
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           '(is_gaia_dr3 is False) and ' +
                                           '(is_gaia_dr2 is False) and ' +
                                           '(is_legacysurvey_dr8 is False) and ' +
                                           '(is_legacysurvey_dr10 is False) and ' +
                                           '(is_panstarrs_dr2 is False) and ' +
                                           '(is_twomass_psc is False)')

            return query

    return FileCarton


def create_table_as(query, table_name, schema=None, temporary=False,
                    database=None, execute=True, overwrite=False,
                    indices=[], analyze=True):
    """Creates a table from a query.

    Parameters
    ----------
    query
        A Peewee ``ModelSelect`` or a string with the query to create a table from.
    table_name
        The name of the table to create.
    schema
        The schema in which to create the table. If ``table_name`` is in the
        form ``schema.table``, the schema parameter is overridden by ``table_name``.
    temporary
        Whether to create a temporary table instead of a persistent one.
    database
        The database connection to use to execute the query. If not passed and the
        query is a ``ModelSelect``, the database will be inherited from the query
        model.
    execute
        Whether to actually execute the query. Requires ``database`` to be passed.
    overwrite
        If `True`, the table will be create even if a table with the same name
        already exists. Requires ``database`` or will be ignored.
    analyze
        Whether to ``VACUUM ANALIZE`` the new table.
        Only relevant if ``execute=True``.
    indices
        List of columns to create indices on. Only relevant if ``execute=True``.

    Returns
    -------
    create_query
        A tuple in which the first element is a Peewee ``Table`` for the created table
        (the table is bound to ``database`` if passed), and the ``CREATE TABLE AS``
        query as a string.

    """

    if '.' in table_name:
        schema, table_name = table_name.split('.')

    if schema is None and temporary is False:
        schema = 'public'
    elif temporary is True:
        schema = None

    path = f"{schema}.{table_name}" if schema else table_name
    create_sql = f'CREATE {"TEMPORARY " if temporary else ""}TABLE {path} AS '

    if database is None and isinstance(query, peewee.ModelSelect):
        database = query.model._meta.database

    if overwrite and database:
        database.execute_sql(f'DROP TABLE IF EXISTS {path};')

    query_sql, params = database.get_sql_context().sql(query).query()
    cursor = database.cursor()
    query_str = cursor.mogrify(query_sql, params).decode()

    if execute:
        if database is None:
            raise RuntimeError('Cannot execute query without a database.')

        with database.atomic():
            database.execute_sql(create_sql + query_sql, params)

        for index in indices:
            if isinstance(index, (list, tuple)):
                index = ','.join(index)
            database.execute_sql(f'CREATE INDEX ON {path} ({index})')

        if analyze:
            database.execute_sql(f'VACUUM ANALYZE {path}')

    table = peewee.Table(table_name, schema=schema).bind(database)

    create_str = create_sql + query_str
    return table, create_str
