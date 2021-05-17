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
                                             CatalogToTIC_v8, Gaia_DR2,
                                             Legacy_Survey_DR8, TIC_v8)

from .base import BaseCarton


def get_file_carton(
        filename,
        carton_name,
        carton_category,
        carton_program,
        replace_carton_name=True):
    """Returns a carton class that creates a carton based on a FITS file."""

    class FileCarton(BaseCarton):
        name = carton_name
        category = carton_category
        program = carton_program

        def __init__(self, targeting_plan, config_file=None, schema=None, table_name=None):

            self._file_path = filename

            self._table = Table.read(self._file_path)
            self._table.convert_bytestring_to_unicode()
            if self._table.masked:
                self._table = self._table.filled()

            # Replace name of carton from file column.
            if 'cartonname' in self._table.columns and replace_carton_name:
                uniq_cname = numpy.unique(self._table['cartonname'])
                if len(uniq_cname) == 1:
                    self.name = uniq_cname[0]

            super().__init__(
                targeting_plan,
                config_file=config_file,
                schema=schema,
                table_name=table_name,
            )

            self._disable_query_log = True

        def build_query(self, version_id, query_region=None):

            self.log.debug(f'Processing file {self._file_path}.')

            gaia_ids = (self._table[self._table['Gaia_DR2_Source_ID'] > 0]
                        ['Gaia_DR2_Source_ID'].tolist())
            ls8_ids = (self._table[self._table['LegacySurvey_DR8_ID'] > 0]
                       ['LegacySurvey_DR8_ID'].tolist())

            vl = peewee.ValuesList(self._table.as_array().tolist(),
                                   columns=self._table.colnames,
                                   alias='fits')

            gid_case = peewee.Case(
                None,
                ((vl.c.Gaia_DR2_Source_ID > 0, vl.c.Gaia_DR2_Source_ID),))
            ls_id_case = peewee.Case(
                None,
                ((vl.c.LegacySurvey_DR8_ID > 0, vl.c.LegacySurvey_DR8_ID),))
            inertial_case = peewee.Case(
                None,
                ((vl.c.inertial.cast('boolean').is_null(), False),),
                vl.c.inertial.cast('boolean'))

            query = (Catalog
                     .select(Catalog.catalogid,
                             gid_case.alias('gaia_source_id'),
                             ls_id_case.alias('ls_id'),
                             vl.c.ra.cast('double precision'),
                             vl.c.dec.cast('double precision'),
                             vl.c.delta_ra.cast('double precision'),
                             vl.c.delta_dec.cast('double precision'),
                             inertial_case.alias('inertial'),
                             vl.c.cadence,
                             vl.c.priority,
                             vl.c.instrument,
                             peewee.Value(0).alias('value'))
                     .join(CatalogToTIC_v8, peewee.JOIN.LEFT_OUTER)
                     .join(TIC_v8, peewee.JOIN.LEFT_OUTER)
                     .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER)
                     .switch(Catalog)
                     .join(CatalogToLegacy_Survey_DR8, peewee.JOIN.LEFT_OUTER)
                     .join(Legacy_Survey_DR8, peewee.JOIN.LEFT_OUTER)
                     .join(vl, on=((vl.c.Gaia_DR2_Source_ID == Gaia_DR2.source_id) |
                                   (vl.c.LegacySurvey_DR8_ID == Legacy_Survey_DR8.ls_id)))
                     .where(Gaia_DR2.source_id.in_(gaia_ids) |
                            Legacy_Survey_DR8.ls_id.in_(ls8_ids))
                     .where(Catalog.version_id == version_id,
                            ((CatalogToTIC_v8.best >> True) | (CatalogToTIC_v8.best.is_null())),
                            ((CatalogToLegacy_Survey_DR8.best >> True) |
                             (CatalogToLegacy_Survey_DR8.best.is_null()))))

            if 'lambda_eff' in self._table.colnames:
                query = query.select_extend(vl.c.lambda_eff.alias('lambda_eff'))

            return query

    return FileCarton
