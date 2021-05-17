#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-07-01
# @Filename: mwm_cb_uvex.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import Catalog, CatalogToTIC_v8, TIC_v8

from target_selection.cartons import BaseCarton


def h2exp(hmag, sn=100, exptime=15.0):
    """ This function takes in a hmag and given signal to noise and spits back
    the required time. Based on Hmag = 11 at S/N 100 in an hour.
    """

    # Scale the hmag based on t = (1 hour)*10^(0.4*(H-11))
    # Then I cut it up into 15 minute exposures.
    time = 60 * (sn**2 / 100.0**2) * 10**(0.4 * (hmag - 11))
    nexp = numpy.array(numpy.round(time / exptime))

    # Min value is 1
    nexp[(nexp == 0)] = 1

    # Set Nan's to nan
    nexp[numpy.isnan(hmag)] = numpy.nan

    return(nexp)


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
    instrument = 'APOGEE'
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

    def post_process(self, model, **kwargs):

        data = numpy.array(model.select(model.catalogid, model.hmag).tuples(),
                           dtype=[('catalogid', numpy.int64),
                                  ('hmag', numpy.float32)])
        n_exp = h2exp(data['hmag'], sn=80)

        values = ((int(data['catalogid'][ii]), 'bright_1x' + str(int(n_exp[ii]))
                   if not numpy.isnan(n_exp[ii]) else None)
                  for ii in range(len(data)))
        vl = peewee.ValuesList(values, columns=('catalogid', 'cadence'), alias='vl')

        (model
         .update(cadence=vl.c.cadence)
         .from_(vl)
         .where(model.catalogid == vl.c.catalogid)).execute()
