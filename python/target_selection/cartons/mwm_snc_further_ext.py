#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2024-07-08
# @Filename: mwm_snc_further_ext.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


class MWM_snc_further_ext_Base_Carton(BaseCarton):
    """
    MWM_snc_further_ext_Base_Carton is a base carton.

    Actual cartons are implemented as subclasses of MWM_SNC_further_ext_Base_Carton
    for different conditions in the WHERE clause.

    https://sdss-wiki.atlassian.net/wiki/spaces/MWM/pages/132579391/A+Further+Extension+of+the+SNC
    """

    def build_query(self, version_id, query_region=None):
        ll = Gaia_DR3.l
        bb = Gaia_DR3.b

        # Dense regions (Galactic plane, SMC, LMC).
        gal_cut = (
            ((ll <= 180) & (bb < (-0.139 * ll + 25)) & (bb > (0.139 * ll - 25)))
            | ((ll > 180) & (bb > (-0.139 * ll + 25)) & (bb < (0.139 * ll - 25)))
            | (fn.sqrt(fn.pow(ll - 303.2, 2) + 2 * fn.pow(bb + 44.4, 2)) < 5)
            | (fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8)
        )

        cte = (
            Gaia_DR3.select(Gaia_DR3.source_id)
            .where(
                Gaia_DR3.parallax < 10,
                Gaia_DR3.parallax - Gaia_DR3.parallax_error > 2,
                Gaia_DR3.phot_g_mean_mag + 5 * fn.log10(Gaia_DR3.parallax / 1000) + 5
                > 19.125 - 0.03225 * (1000 / Gaia_DR3.parallax),
            )
            .cte("plx_mag_cte", materialized=True)
        )

        query = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.l,
                Gaia_DR3.b,
                Gaia_DR3.parallax,
                Gaia_DR3.parallax_error,
                Gaia_DR3.phot_g_mean_mag.alias("gaia_dr3_phot_g_mean_mag"),
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag,
                Gaia_DR3.ruwe,
                Gaia_DR3.phot_bp_rp_excess_factor,
                Gaia_DR3.astrometric_excess_noise,
            )
            .join(CatalogToGaia_DR3)
            .join_from(Gaia_DR3, cte, on=(Gaia_DR3.source_id == cte.c.source_id))
            .where(
                (
                    (Gaia_DR3.astrometric_excess_noise < 2)
                    & (Gaia_DR3.ruwe < 1.2)
                    & (
                        Gaia_DR3.phot_bp_rp_excess_factor
                        > 1
                        + 0.015 * fn.pow(Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag, 2)
                    )
                    & (
                        Gaia_DR3.phot_bp_rp_excess_factor
                        < 1.3
                        + 0.06 * fn.pow(Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag, 2)
                    )
                    & gal_cut
                )
                | ~gal_cut
            )
            .where(CatalogToGaia_DR3.version_id == version_id, CatalogToGaia_DR3.best >> True)
            .with_cte(cte)
        )
        return query


class Openfibertargets_mwm_snc_further_ext_apogee_bright_Carton(MWM_snc_further_ext_Base_Carton):
    """Shorthand name: openfibertargets_mwm_snc_further_ext_apogee_bright
    REPLACES openfibertargets_nov2020_24
    Return columns: Gaia Galactic longitude and latitude, parallax, parallax_error,
        G, BP, RP, ruwe, phot_bp_rp_excess_factor, astrometric_excess_noise, and 2MASS H
    Metadata: can_offset=True
    Priority: 6085
    Cadence: bright_1x1
    Instrument: APOGEE
    Program: open_fiber
    Lead contact: Ilija Medan
    """

    name = "openfibertargets_mwm_snc_further_ext_apogee_bright"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    program = "open_fiber"
    mapper = "MWM"
    priority = 6085
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = (
            query.switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC)
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                TwoMassPSC.h_m < 11,
            )
        )

        return query


class Openfibertargets_mwm_snc_further_ext_apogee_dark_Carton(MWM_snc_further_ext_Base_Carton):
    """Shorthand name: openfibertargets_mwm_snc_further_ext_apogee_dark
    REPLACES openfibertargets_nov2020_24
    Return columns: Gaia Galactic longitude and latitude, parallax, parallax_error,
        G, BP, RP, ruwe, phot_bp_rp_excess_factor, astrometric_excess_noise, and 2MASS H
    Metadata: can_offset=True
    Priority: 6085
    Cadence: dark_1x1
    Instrument: APOGEE
    Program: open_fiber
    Lead contact: Ilija Medan
    """

    name = "openfibertargets_mwm_snc_further_ext_apogee_dark"
    category = "science"
    instrument = "APOGEE"
    cadence = "dark_1x1"
    program = "open_fiber"
    mapper = "MWM"
    priority = 6085
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = (
            query.switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC)
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                TwoMassPSC.h_m < 11,
            )
        )

        return query


class Openfibertargets_mwm_snc_further_ext_boss_Carton(MWM_snc_further_ext_Base_Carton):
    """Shorthand name: openfibertargets_mwm_snc_further_ext_boss
    REPLACES openfibertargets_nov2020_24
    Return columns: Gaia Galactic longitude and latitude, parallax, parallax_error,
    G, BP, RP, ruwe, phot_bp_rp_excess_factor, astrometric_excess_noise, and 2MASS H
    Metadata: can_offset=True
    Priority: 6086
    Cadence: bright_1x1 for stars with G < 16 and dark_1x1 for stars with G > 16
    Instrument: BOSS
    Program: open_fiber
    Lead contact: Ilija Medan
    """

    name = "openfibertargets_mwm_snc_further_ext_boss"
    category = "science"
    instrument = "BOSS"
    cadence = None  # cadence is set in post_process()
    program = "open_fiber"
    mapper = "MWM"
    priority = 6086
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)

        query = query.where(
            Gaia_DR3.phot_g_mean_mag > 10,
            Gaia_DR3.phot_g_mean_mag < 20,
            Gaia_DR3.radial_velocity.is_null(),
        )

        return query

    def post_process(self, model):
        # Note that self is a parameter of post_process() above but
        # self is not a parameter of post_process() below
        # since we are using super().
        super().post_process(model)

        # set cadence bright_1x1 for stars with G < 16 and dark_1x1 for stars with G > 16

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr3_phot_g_mean_mag from "
            + " sandbox.temp_openfibertargets_mwm_snc_further_ext_boss ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_g = output[i][1]

            if current_g < 16:
                current_cadence = "bright_1x1"
            else:
                current_cadence = "dark_1x1"

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_openfibertargets_mwm_snc_further_ext_boss "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )
