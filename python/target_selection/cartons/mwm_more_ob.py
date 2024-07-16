#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2024-07-03
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


# from target_selection.exceptions import TargetSelectionError


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class Openfibertargets_mwm_more_ob_boss_Carton(BaseCarton):
    """openfibertargets_mwm_more_ob_boss
    Shorthand name: openfibertargets_mwm_more_ob_boss
    Link: More B stars
    Simplified Description of selection criteria:

    SELECT g.source_id
    FROM gaiadr3.gaia_source as g
    USING (source_id)
    INNER JOIN gaiadr3.tmass_best_neighbour AS xmatch
    ON g.source_id = xmatch.source_id
    INNER JOIN gaiadr1.tmass_original_valid AS tm
    ON tm.tmass_oid = xmatch.tmass_oid

    WHERE parallax < power (10, ((10. - tm.k_m + 0.) / 5.))

    AND parallax > power(10, ((10. - tm.k_m - 0.61) / 5.))

    AND g.phot_g_mean_mag < 16.
    AND tm.j_m - tm.k_m - 0.25 * (g.phot_g_mean_mag - tm.k_m) < 0.10
    AND tm.j_m - tm.k_m - 0.25 * (g.phot_g_mean_mag - tm.k_m) > -0.30
    AND tm.j_m - tm.h_m < 0.15 * (g.phot_g_mean_mag - tm.k_m) + 0.05
    AND tm.j_m - tm.h_m > 0.15 * (g.phot_g_mean_mag - tm.k_m) - 0.15
    AND g.phot_g_mean_mag > 2 * (g.phot_g_mean_mag - tm.k_m) + 3.0

    Return columns: Gaia DR3: source_id, ra, dec, parallax, pmra, pmdec,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag; 2MASS: j_m, h_m, k_m

    Metadata: can_offset=True
    Priority: 6085
    Cadence: bright_1x1
    Instrument: BOSS
    Program: open_fiber

    Lead contact: Eleonora Zari, Jaime Villaseñor

    """

    name = "openfibertargets_mwm_more_ob_boss"
    category = "science"
    instrument = "BOSS"
    cadence = "bright_1x1"
    program = "open_fiber"
    mapper = "MWM"
    priority = 6085
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.parallax,
                Gaia_DR3.pmra,
                Gaia_DR3.pmdec,
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag,
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                Gaia_DR3.phot_g_mean_mag < 16,
                Gaia_DR3.parallax < peewee.fn.power(10, ((10.0 - TwoMassPSC.k_m + 0.00) / 5.0)),
                Gaia_DR3.parallax > peewee.fn.power(10, ((10.0 - TwoMassPSC.k_m - 0.61) / 5.0)),
                TwoMassPSC.j_m
                - TwoMassPSC.k_m
                - 0.25 * (Gaia_DR3.phot_g_mean_mag - TwoMassPSC.k_m)
                < 0.10,
                TwoMassPSC.j_m
                - TwoMassPSC.k_m
                - 0.25 * (Gaia_DR3.phot_g_mean_mag - TwoMassPSC.k_m)
                > -0.30,
                TwoMassPSC.j_m - TwoMassPSC.h_m
                < 0.15 * (Gaia_DR3.phot_g_mean_mag - TwoMassPSC.k_m) + 0.05,
                TwoMassPSC.j_m - TwoMassPSC.h_m
                > 0.15 * (Gaia_DR3.phot_g_mean_mag - TwoMassPSC.k_m) - 0.15,
                Gaia_DR3.phot_g_mean_mag > 2 * (Gaia_DR3.phot_g_mean_mag - TwoMassPSC.k_m) + 3.0,
            )
        )

        # There can be cases in which the same catalogid has multiple entries
        # in a catalog_to_x table since the same physical object
        # may match with multiple catalogids.
        # Hence, we have the below code in the above query.
        #                CatalogToGaia_DR3.best >> True,
        #                CatalogToTwoMassPSC.best >> True,

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query
