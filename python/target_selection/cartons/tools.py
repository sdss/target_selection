#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2021-04-29
# @Filename: tools.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import peewee
from astropy.table import Table

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToLegacy_Survey_DR8,
                                             CatalogToPanstarrs1,
                                             CatalogToTIC_v8, Gaia_DR2,
                                             Legacy_Survey_DR8, Panstarrs1,
                                             TIC_v8)
from sdssdb.utils.ingest import copy_data, create_model_from_table

from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError
from target_selection.utils import vacuum_table


def get_file_carton(
        filename,
        carton_name,
        carton_category,
        carton_program,
        replace_carton_name=True):
    """Returns a carton class that creates a carton based on a FITS file.
    The FITS file is located in the below location which is specified in
    python/config/target_selection.yml.
    open_fiber_path: $CATALOGDB_DIR/../open_fiber/0.5.0/
    """

    class FileCarton(BaseCarton):
        name = carton_name
        category = carton_category
        program = carton_program

        def __init__(self, targeting_plan, config_file=None, schema=None, table_name=None):

            self._file_path = filename

            self._table = Table.read(self._file_path)
            # self._table.convert_bytestring_to_unicode()
            if self._table.masked:
                self._table = self._table.filled()

            # Replace name of carton from file column.
            if 'cartonname' in self._table.columns and replace_carton_name:
                uniq_cname = numpy.unique(self._table['cartonname'])
                if len(uniq_cname) == 1:
                    self.name = uniq_cname[0].lower()

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

            # Create model for sandbox table from FITS table columns.
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
                 .join(CatalogToTIC_v8, peewee.JOIN.LEFT_OUTER)
                 .join(TIC_v8, peewee.JOIN.LEFT_OUTER)
                 .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER)
                 .join(temp,
                       on=(temp.Gaia_DR2_Source_ID == Gaia_DR2.source_id))
                 .switch(Catalog)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        (CatalogToTIC_v8.best >> True) |
                        CatalogToTIC_v8.best.is_null(),
                        Catalog.version_id == version_id))

            query_legacysurvey_dr8 = \
                (query_common
                 .join(CatalogToLegacy_Survey_DR8, peewee.JOIN.LEFT_OUTER)
                 .join(Legacy_Survey_DR8, peewee.JOIN.LEFT_OUTER)
                 .join(temp,
                       on=(temp.LegacySurvey_DR8_ID == Legacy_Survey_DR8.ls_id))
                 .switch(Catalog)
                 .where(CatalogToLegacy_Survey_DR8.version_id == version_id,
                        (CatalogToLegacy_Survey_DR8.best >> True) |
                        CatalogToLegacy_Survey_DR8.best.is_null(),
                        Catalog.version_id == version_id))

            query_panstarrs_dr2 = \
                (query_common
                 .join(CatalogToPanstarrs1, peewee.JOIN.LEFT_OUTER)
                 .join(Panstarrs1, peewee.JOIN.LEFT_OUTER)
                 .join(temp,
                       on=(temp.PanSTARRS_DR2_ID == Panstarrs1.catid_objid))
                 .switch(Catalog)
                 .where(CatalogToPanstarrs1.version_id == version_id,
                        (CatalogToPanstarrs1.best >> True) |
                        CatalogToPanstarrs1.best.is_null(),
                        Catalog.version_id == version_id))

            if len(self._table[self._table['Gaia_DR2_Source_ID'] > 0]):
                is_gaia_dr2 = True
            else:
                is_gaia_dr2 = False

            if len(self._table[self._table['LegacySurvey_DR8_ID'] > 0]):
                is_legacysurvey_dr8 = True
            else:
                is_legacysurvey_dr8 = False

            if len(self._table[self._table['PanSTARRS_DR2_ID'] > 0]):
                is_panstarrs_dr2 = True
            else:
                is_panstarrs_dr2 = False

            # We consider all 8 cases.
            if((is_gaia_dr2 is True) and
               (is_legacysurvey_dr8 is True) and
               (is_panstarrs_dr2 is True)):
                # query is a SQL union of the three queries on the RHS
                query = query_gaia_dr2 | query_legacysurvey_dr8 | query_panstarrs_dr2

            elif((is_gaia_dr2 is True) and
                 (is_legacysurvey_dr8 is True) and
                 (is_panstarrs_dr2 is False)):
                query = query_gaia_dr2 | query_legacysurvey_dr8

            elif((is_gaia_dr2 is True) and
                 (is_legacysurvey_dr8 is False) and
                 (is_panstarrs_dr2 is True)):
                query = query_gaia_dr2 | query_panstarrs_dr2

            elif((is_gaia_dr2 is True) and
                 (is_legacysurvey_dr8 is False) and
                 (is_panstarrs_dr2 is False)):
                query = query_gaia_dr2

            elif((is_gaia_dr2 is False) and
                 (is_legacysurvey_dr8 is True) and
                 (is_panstarrs_dr2 is True)):
                query = query_legacysurvey_dr8 | query_panstarrs_dr2

            elif((is_gaia_dr2 is False) and
                 (is_legacysurvey_dr8 is True) and
                 (is_panstarrs_dr2 is False)):
                query = query_legacysurvey_dr8

            elif((is_gaia_dr2 is False) and
                 (is_legacysurvey_dr8 is False) and
                 (is_panstarrs_dr2 is True)):
                query = query_panstarrs_dr2

            elif((is_gaia_dr2 is False) and
                 (is_legacysurvey_dr8 is False) and
                 (is_panstarrs_dr2 is False)):
                raise TargetSelectionError('error in get_file_carton():' +
                                           '(is_gaia_dr2 is False) and ' +
                                           '(is_legacysurvey_dr8 is False) and ' +
                                           '(is_panstarrs_dr2 is False)')
                query = None
            else:
                # we will not get here since we have
                # considered all 8 cases above
                query = None

            if 'lambda_eff' in self._table.colnames:
                query = query.select_extend(temp.lambda_eff.alias('lambda_eff'))

            return query

    return FileCarton
