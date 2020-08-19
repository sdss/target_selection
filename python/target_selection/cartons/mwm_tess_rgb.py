#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-07-01
# @Filename: mwm_cb_uvex.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import Catalog, CatalogToTIC_v8, TIC_v8

from . import BaseCarton


class MWM_TESS_RGB_Carton(BaseCarton):
    """MWM TESS RGB Carton.

    Definition:

        - Jmag - Kmag > 0.5 (get red stars)
        - Hmag < 12 (get needed SNR w/APOGEE)
        - Tmag < 13 (get detections of oscillations) [Tmag := TESS magnitude]
        - M_H (absolute H band magnitude) < 1 (get giants)
            where: MH = Hmag - 10 + 5.0 * log10(parallax)
            where the parallax is in mas.

    For v0 we'll also add |b|>20

    """

    name = 'mwm_tessrgb_core'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_tessrgb'
    cadence = None
    priority = 2800

    def build_query(self, version_id, query_region=None):

        MH = TIC_v8.hmag - 10 + 5 * fn.log(TIC_v8.plx)

        query = (TIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         TIC_v8.hmag,
                         TIC_v8.jmag,
                         TIC_v8.kmag,
                         TIC_v8.tmag,
                         TIC_v8.plx)
                 .join(CatalogToTIC_v8)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where((TIC_v8.jmag - TIC_v8.kmag) > 0.5)
                 .where(((TIC_v8.plx > 0) & (MH < 1)) | (TIC_v8.plx < 0),
                        fn.abs(TIC_v8.gallat) > 20)
                 .where(TIC_v8.hmag < 12)
                 .where(TIC_v8.tmag < 13))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
