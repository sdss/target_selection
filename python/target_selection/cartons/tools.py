#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2021-04-29
# @Filename: tools.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os
import warnings

import numpy
import peewee
from astropy.table import Table

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToLegacy_Survey_DR8,
                                             CatalogToPanstarrs1,
                                             CatalogToTIC_v8, Gaia_DR2,
                                             Legacy_Survey_DR8, Panstarrs1,
                                             TIC_v8, TwoMassPSC)
from sdssdb.utils.ingest import copy_data, create_model_from_table

from target_selection.cartons import BaseCarton
from target_selection.exceptions import (TargetSelectionError,
                                         TargetSelectionUserWarning)
from target_selection.utils import vacuum_table


def get_file_carton(filename):
    """Returns a carton class that creates a carton based on a FITS file.
    The FITS file is located in the open_fiber_path which is specified in
    python/config/target_selection.yml.
    The list of FITS files to be loaded is specified in the
    file open_fiber_file_list.txt which is in the directory open_fiber_path.
    """

    class FileCarton(BaseCarton):

        def __init__(self, targeting_plan, config_file=None, schema=None, table_name=None):

            self._file_path = filename

            self._table = Table.read(self._file_path)
            # historical
            # self._table.convert_bytestring_to_unicode()
            if self._table.masked:
                self._table = self._table.filled()

            unique_cartonname = numpy.unique(self._table['cartonname'])
            if len(unique_cartonname) == 1:
                self.name = unique_cartonname[0].lower()
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           'contains more than one cartonname')

            # The valid_program list is from the output of the below command.
            # select distinct(program) from targetdb.carton order by program;
            valid_program = [
                'bhm_aqmes',
                'bhm_csc',
                'bhm_filler',
                'bhm_rm',
                'bhm_spiders',
                'commissioning',
                'mwm_cb',
                'mwm_dust',
                'mwm_erosita',
                'mwm_filler',
                'mwm_galactic',
                'mwm_gg',
                'mwm_halo',
                'mwm_legacy',
                'mwm_ob',
                'mwm_planet',
                'mwm_rv',
                'mwm_snc',
                'mwm_tess_ob',
                'mwm_tessrgb',
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
                                               'contains invalid category = ' +
                                               self.category)
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           'contains more than one category')

            unique_program = numpy.unique(self._table['program'])
            if len(unique_program) == 1:
                self.program = unique_program[0].lower()
                if (self.program not in valid_program):
                    raise TargetSelectionError('error in get_file_carton(): ' +
                                               filename +
                                               'contains invalid program = ' +
                                               self.program)
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           'contains more than one program')

            unique_mapper = numpy.unique(self._table['mapper'])
            if len(unique_mapper) == 1:
                # We do not use lower() for mapper because
                # allowed values for mapper are '' or 'MWM' or 'BHM'.
                self.mapper = unique_mapper[0]
                if (self.mapper not in valid_mapper):
                    raise TargetSelectionError('error in get_file_carton(): ' +
                                               filename +
                                               'contains invalid mapper = ' +
                                               self.mapper)
                if(self.mapper == ''):
                    self.mapper = None
            else:
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           filename +
                                           'contains more than one mapper')

            basename_fits = os.path.basename(filename)
            basename_parts = os.path.splitext(basename_fits)
            basename = basename_parts[0]
            carton_name_from_filename = basename.lower()

            if (self.name != carton_name_from_filename):
                warnings.warn('filename parameter of get_file_carton() and ' +
                              'cartonname in FITS file do not match.',
                              TargetSelectionUserWarning)
                warnings.warn('carton_name_from_filename = ' + carton_name_from_filename +
                              ' cartonname = ' + self.name,
                              TargetSelectionUserWarning)

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

            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("Gaia_DR2_Source_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("LegacySurvey_DR8_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("PanSTARRS_DR2_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("TwoMASS_ID")')
            vacuum_table(self.database, temp_table, vacuum=False, analyze=True)

            inertial_case = peewee.Case(
                None,
                ((temp.inertial.cast('boolean').is_null(), False),),
                temp.inertial.cast('boolean'))

            query_common = (Catalog
                            .select(Catalog.catalogid,
                                    temp.Gaia_DR2_Source_ID.alias('gaia_source_id'),
                                    temp.LegacySurvey_DR8_ID.alias('ls_id'),
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
                                    peewee.Value(0).alias('value'))
                            .distinct(Catalog.catalogid))

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
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join(temp,
                       on=(temp.TwoMASS_ID == TwoMassPSC.designation))
                 .switch(Catalog)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        (CatalogToTIC_v8.best >> True) |
                        CatalogToTIC_v8.best.is_null(),
                        Catalog.version_id == version_id))

            len_table = len(self._table)

            len_gaia_dr2 =\
                len(self._table[self._table['Gaia_DR2_Source_ID'] > 0])

            len_legacysurvey_dr8 =\
                len(self._table[self._table['LegacySurvey_DR8_ID'] > 0])

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
            if ((len_gaia_dr2 + len_legacysurvey_dr8 +
                 len_panstarrs_dr2 + len_twomass_psc) != len_table):
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           '(len_gaia_dr2 + len_legacysurvey_dr8 + ' +
                                           'len_panstarrs_dr2 + len_twomass_psc) != ' +
                                           'len_table')

            if (len_gaia_dr2 > 0):
                is_gaia_dr2 = True
            else:
                is_gaia_dr2 = False

            if (len_legacysurvey_dr8 > 0):
                is_legacysurvey_dr8 = True
            else:
                is_legacysurvey_dr8 = False

            if (len_panstarrs_dr2 > 0):
                is_panstarrs_dr2 = True
            else:
                is_panstarrs_dr2 = False

            if (len_twomass_psc > 0):
                is_twomass_psc = True
            else:
                is_twomass_psc = False

            query = None

            if(is_gaia_dr2 is True):
                if(query is None):
                    query = query_gaia_dr2
                else:
                    query = query | query_gaia_dr2

            if(is_legacysurvey_dr8 is True):
                if(query is None):
                    query = query_legacysurvey_dr8
                else:
                    query = query | query_legacysurvey_dr8

            if(is_panstarrs_dr2 is True):
                if(query is None):
                    query = query_panstarrs_dr2
                else:
                    query = query | query_panstarrs_dr2

            if(is_twomass_psc is True):
                if(query is None):
                    query = query_twomass_psc
                else:
                    query = query | query_twomass_psc

            if(query is None):
                # At least one of the four boolean variables above
                # must be True, so we should not get here.
                raise TargetSelectionError('error in get_file_carton(): ' +
                                           '(is_gaia_dr2 is False) and ' +
                                           '(is_legacysurvey_dr8 is False) and ' +
                                           '(is_panstarrs_dr2 is False) and ' +
                                           '(is_twomass_psc is False)')

            if 'lambda_eff' in self._table.colnames:
                query = query.select_extend(temp.lambda_eff.alias('lambda_eff'))

            return query

    return FileCarton
