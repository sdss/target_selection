#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2021-01-17
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import MWM_TESS_OB, Catalog, CatalogToTIC_v8, Gaia_DR2, TIC_v8

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# For example
# peewee Model name ---> postgres table name
# Gaia_DR2 --->'catalogdb.gaia_dr2_source'
# MWM_TESS_OB --->'catalogdb.mwm_tess_ob'


class MWM_TESS_OB_Carton(BaseCarton):
    """MWM TESS OB
    Owner: Nathan De Lee

    Shorthand name: mwm_tess_ob

    Mock Target catalogs with GAIA DR2 ids for the TESS
    CVZ OBAF stars

    Simplified Description of selection criteria:
    a random subsample of TIC targets in both TESS CVZs,
    selected according to type of their variability
    (i.e., if (var_type = eclipsing_binary .and. var_type = pulsator
     .and. coord = TESS_CVZ) select )

    Wiki page:

    Additional source catalogs needed:
    There necessary variability information is not currently available,
    so we are going to use mock targets from the following
    two lists of Gaia DR2 IDs:

    mwm_tess_ob_8x1_GaiaID.csv
    mwm_tess_ob_8x3_GaiaID.csv

    Additional cross-matching needed: None

    return columns:  Gaia ID, Hmag

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    There are three cadences that will be used by these targets.

    For all targets in mwm_tess_ob_8x1_GaiaID.csv use
        instrument =  "APOGEE"
        cadence = "bright_8x1"

    For targets in mwm_tess_ob_8x3_GaiaID.csv
    if H < 11 then use
        instrument = "APOGEE"
        cadence = "bright_8x2"

    if H > =11 then use
        instrument = "APOGEE"
        cadence = "bright_8x4"

    Pseudo SQL (optional):

    Implementation:

    Non-SQL implementation:

    Lead Contact: Andrew Tkachenko
    """

    name = "mwm_tess_ob"
    category = "science"
    instrument = None  # assigned in query
    cadence = None  # assigned in query
    program = "mwm_tessob"
    mapper = "MWM"
    priority = 2200  # from the table in the V0.5 carton wiki page

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                Gaia_DR2.phot_g_mean_mag.alias("gaia_dr2_g"),
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag,
                Gaia_DR2.parallax,
                MWM_TESS_OB.ra,
                MWM_TESS_OB.dec,
                MWM_TESS_OB.h_mag,
                MWM_TESS_OB.instrument,
                MWM_TESS_OB.cadence,
            )
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .join(MWM_TESS_OB, on=(Gaia_DR2.source_id == MWM_TESS_OB.gaia_dr2_id))
            .where(CatalogToTIC_v8.version_id == version_id, CatalogToTIC_v8.best >> True)
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
