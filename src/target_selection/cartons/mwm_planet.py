#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2022-12-11
# @Filename: mwm_tess_planets.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToTIC_v8,
    TESS_TOI_v1,
    TIC_v8,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


def h2exp(hmag, sn=100, exptime=15.0):
    """
    This function takes in a hmag and given signal to noise and spits back
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


class MWM_TESS_2min_apogee_Carton(BaseCarton):
    """TESS Planets.
    Definition:
        All planet-candidate TOIs + Community TOIs (CTIOs), TESS
        Candidate-Target List (CTL) selected by a TIC Priority cut.
        All cut between 7 < H < 12.

    All the candidates are compiled in catalogdb.tess_toi_v1. Note that this list
    contains duplicates on ticid.

    We changed the carton class name, the variable 'name'
    and the temp table name in post_process().

    The old name of this cartoon is 'mwm_tess_planet'.
    """

    name = "mwm_tess_2min_apogee"
    program = "mwm_planet"
    category = "science"
    mapper = "MWM"
    instrument = "APOGEE"
    cadence = None  # cadence is set in post_process()
    priority = 2810
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            TESS_TOI_v1.select(
                CatalogToTIC_v8.catalogid,
                TESS_TOI_v1.ticid,
                TESS_TOI_v1.tess_disposition,
                TwoMassPSC.h_m.alias("hmag"),
                TESS_TOI_v1.target_type.alias("tess_target_type"),
            )
            .join(TIC_v8, on=(TESS_TOI_v1.ticid == TIC_v8.id))
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .switch(TIC_v8)
            .join(CatalogToTIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                TwoMassPSC.h_m > 7,
                TwoMassPSC.h_m < 12,
            )
            .distinct([TESS_TOI_v1.ticid])
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

    def post_process(self, model):
        cursor = self.database.execute_sql(
            "select catalogid, hmag from " + " sandbox.temp_mwm_tess_2min_apogee ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_hmag = float(output[i][1])

            numexp = h2exp(current_hmag, sn=80)

            # for numexp == 1, we use bright_1x1 since
            # there is no cadence bright_flexible_1x1
            if numexp == 1:
                current_cadence = "bright_1x1"
            else:
                current_cadence = "bright_flexible_" + str(int(numexp)) + "x1"

            self.database.execute_sql(
                " update sandbox.temp_mwm_tess_2min_apogee "
                + " set cadence = '"
                + current_cadence
                + "'"
                " where catalogid = " + str(current_catalogid) + ";"
            )
