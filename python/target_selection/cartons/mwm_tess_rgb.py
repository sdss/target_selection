#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-07-01
# @Filename: mwm_tess_rgb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTIC_v8,
    Gaia_DR3,
    TIC_v8,
)

from target_selection.cartons import BaseCarton


def h2exp(hmag, sn=100, exptime=15.0):
    """This function takes in a hmag and given signal to noise and spits back
    the required time. Based on Hmag = 11 at S/N 100 in an hour.
    """

    # Scale the hmag based on t = (1 hour)*10^(0.4*(H-11))
    # Then I cut it up into 15 minute exposures.
    time = 60 * (sn**2 / 100.0**2) * 10 ** (0.4 * (hmag - 11))
    nexp = numpy.array(numpy.round(time / exptime))

    # Min value is 1
    nexp[(nexp == 0)] = 1

    # Set Nan's to nan
    nexp[numpy.isnan(hmag)] = numpy.nan

    return nexp


class MWM_TESS_RGB_apogee_Carton(BaseCarton):
    """MWM TESS RGB apogee Carton.

    Definition:

    For v0.5
        - Jmag - Kmag > 0.5 (get red stars)
        - Hmag < 12 (get needed SNR w/APOGEE)
        - Tmag < 13 (get detections of oscillations) [Tmag := TESS magnitude]
        - M_H (absolute H band magnitude) < 1 (get giants)
            where: MH = Hmag - 10 + 5.0 * log10(parallax)
            where the parallax is in mas.

    For v1
    remove the |b| > 20 criterion
    change the J-K limit to 0.3
    change the MH, Tmag criteria as described below
    h2exp remains the same, yes!

    Simplified Description of selection criteria:

    Jmag - Kmag > 0.3,
    Hmag < 12,
    MH <= 3 where MH=Hmag-10+5.0*log10(parallax)
    where the parallax is in mas.

    if MH < 1 then Tmag must be  < 13 else Tmag must be < -5*MH+18

    Gaia DR2 parameters to be converted to Gaia DR3:
    use Gaia DR3 parallaxes instead of Gaia DR2

    use joins tic_v8 → catalog_to_tic_v8 → catalog → catalog_to_gaia_dr3 → gaia_dr3

    Return columns: Unchanged

    Metadata:
    priority=2820
    cadence as before
    can_offset=True
    instrument=APOGEE

    """

    name = "mwm_tess_rgb_apogee"  # there is an underscore between tess and rgb
    mapper = "MWM"
    category = "science"
    program = "mwm_tessrgb"  # there is no underscore between tess and rgb
    instrument = "APOGEE"
    cadence = None
    priority = 2820
    can_offset = True

    def build_query(self, version_id, query_region=None):
        MH = TIC_v8.hmag - 10 + 5 * fn.log(Gaia_DR3.parallax)

        # v0.5
        #         query = (TIC_v8
        #                  .select(CatalogToTIC_v8.catalogid,
        #                          TIC_v8.hmag,
        #                          TIC_v8.jmag,
        #                          TIC_v8.kmag,
        #                          TIC_v8.tmag,
        #                          TIC_v8.plx)
        #                  .join(CatalogToTIC_v8)
        #                  .where(CatalogToTIC_v8.version_id == version_id,
        #                         CatalogToTIC_v8.best >> True)
        #                  .where((TIC_v8.jmag - TIC_v8.kmag) > 0.5)
        #                  .where(((TIC_v8.plx > 0) & (MH < 1)) | (TIC_v8.plx < 0),
        #                         fn.abs(TIC_v8.gallat) > 20)
        #                  .where(TIC_v8.hmag < 12)
        #                  .where(TIC_v8.tmag < 13))

        query = (
            Catalog.select(
                CatalogToTIC_v8.catalogid,
                TIC_v8.id,
                TIC_v8.hmag,
                TIC_v8.jmag,
                TIC_v8.kmag,
                TIC_v8.tmag,
                Gaia_DR3.parallax,
            )
            .join(CatalogToTIC_v8, on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .switch(Catalog)
            .join(CatalogToGaia_DR3, on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                CatalogToGaia_DR3.best >> True,
                (TIC_v8.jmag - TIC_v8.kmag) > 0.3,
                TIC_v8.hmag < 12,
                Gaia_DR3.parallax > 0,
                ((MH <= 1) & (TIC_v8.tmag < 13))
                | ((MH > 1) & (MH <= 3) & (TIC_v8.tmag < (-5 * MH + 18))),
            )
        )

        if query_region:
            query = query.join_from(CatalogToTIC_v8, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query

    def post_process(self, model, **kwargs):
        data = numpy.array(
            model.select(model.catalogid, model.hmag).tuples(),
            dtype=[("catalogid", numpy.int64), ("hmag", numpy.float32)],
        )
        n_exp = h2exp(data["hmag"], sn=80)

        # Below old code is for historical reference.
        # Below values1 is a generator.
        # values1 = ((int(data['catalogid'][ii]),
        #          'bright_flexible_' + str(int(n_exp[ii])) + 'x1'
        #           if not numpy.isnan(n_exp[ii]) else None)
        #          for ii in range(len(data)))

        # We use name 'values1' instead of 'values' since values()
        # is a Python built-in function.
        #
        # For int(n_exp[ii]) == 1, we use bright_1x1 since
        # there is no cadence bright_flexible_1x1
        #
        # Below values1 is a list and then later we convert it to a tuple.
        values1 = [None] * len(data)
        for ii in range(len(data)):
            if not numpy.isnan(n_exp[ii]):
                if int(n_exp[ii]) == 1:
                    current_cadence = "bright_1x1"
                else:
                    current_cadence = "bright_flexible_" + str(int(n_exp[ii])) + "x1"
                values1[ii] = (int(data["catalogid"][ii]), current_cadence)
            else:
                values1[ii] = None

        values1 = tuple(values1)
        vl = peewee.ValuesList(values1, columns=("catalogid", "cadence"), alias="vl")

        (
            model.update(cadence=vl.c.cadence).from_(vl).where(model.catalogid == vl.c.catalogid)
        ).execute()
