#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_cb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGUVCat,
                                             CatalogToTIC_v8,
                                             GeometricDistances_Gaia_DR2,
                                             GUVCat, TIC_v8)

from . import BaseCarton


class MWM_CB_300_Carton(BaseCarton):
    """MWM Compact Binaries 300pc.

    Definition: Cross-match Gaia & Bailer-Jones distances by source_id,
    cross-match with GALEX including proper motion corrections. Define a
    single linear relation in absolute FUV magnitude vs FUV - NUV,
    select all objects with distances less than 300pc.

    SQL:
        (FUVmag - 5 * log10(r_est/10)) < 14 * (FUVmag-NUVmag) - 46 & r_est <300

    """

    name = 'mwm_cb_300pc'
    mapper = 'MWM'
    category = 'science'
    program = 'CB'

    def build_query(self, version_id, query_region=None):

        GD = GeometricDistances_Gaia_DR2

        FUV = GUVCat.fuv_mag
        NUV = GUVCat.nuv_mag

        FUV_abs = FUV - 5 * peewee.fn.log(GD.r_est / 10)

        query = (Catalog
                 .select(Catalog.catalogid,
                         GD.source_id,
                         GD.r_est,
                         FUV, NUV)
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(GD, on=(TIC_v8.gaia_int == GD.source_id))
                 .join_from(Catalog, CatalogToGUVCat)
                 .join(GUVCat)
                 .where(GD.r_est < 300,
                        FUV_abs < 14 * (FUV - NUV) - 46,
                        CatalogToGUVCat.version_id == version_id,
                        CatalogToGUVCat.best >> True,
                        CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                           Catalog.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
