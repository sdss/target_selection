#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-02
# @Filename: mwm_snc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import random

from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


class MWM_SNC_100pc_Carton(BaseCarton):
    """MWM Solar Neighbourhood Census 100pc.

    Definition:

        - Divide the sky in crowded/not crowded regions via simple geometric
          cuts in (l,b).
        - Select all stars that fall within their one-sigma uncertainty within
          100pc in the not crowded regions.
        - Select all stars that fall within their one-sigma uncertainty
          within 100pc and have an astrometric_excess_noise < 2 in the
          crowded regions.

    Pseudo SQL:

        parallax-parallax_error>10 &&
        ( (astrometric_excess_noise < 2 &&
           ( (l <= 180 && b < -0.139 * l + 25 && b > 0.139 * l - 25) ||
             (l > 180 && b > -0.139 * l + 25 && b < 0.139 * l - 25) ||
             (sqrt(pow(l - 303.2, 2) + 2 * pow(b + 44.4, 2)) < 5) ||
             (sqrt(pow(l - 280.3, 2) + 2 * pow(b + 33.0, 2)) < 8) )
        ) ||
        !( (l <= 180 && b < -0.139 * l + 25 && b > 0.139 * l - 25) ||
           (l > 180 && b > -0.139 * l + 25 && b < 0.139 * l - 25) ||
           (sqrt(pow(l - 303.2, 2) + 2 * pow(b + 44.4, 2)) < 5) ||
           (sqrt(pow(l - 280.3, 2) + 2 * pow(b + 33.0, 2)) < 8) )
        )

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

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
            .where(Gaia_DR3.parallax > 0, (Gaia_DR3.parallax - Gaia_DR3.parallax_error) > 10)
            .cte("plx_cte", materialized=True)
        )

        query = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra,
                Gaia_DR3.dec,
                Gaia_DR3.phot_g_mean_mag,
                ll,
                bb,
            )
            .join(CatalogToGaia_DR3)
            .join_from(Gaia_DR3, cte, on=(Gaia_DR3.source_id == cte.c.source_id))
            .where(((Gaia_DR3.astrometric_excess_noise < 2) & gal_cut) | ~gal_cut)
            .where(CatalogToGaia_DR3.version_id == version_id, CatalogToGaia_DR3.best >> True)
            .with_cte(cte)
        )

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


class MWM_SNC_100pc_APOGEE_Carton(MWM_SNC_100pc_Carton):
    """MWM SNC 100pc targets to be observed with APOGEE."""

    name = "mwm_snc_100pc_apogee"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2705
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region=query_region)
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


class MWM_SNC_100pc_BOSS_Carton(MWM_SNC_100pc_Carton):
    """MWM SNC 100pc targets to be observed with BOSS.

    Shorthand name:  mwm_snc_100pc_boss
    https://sdss-wiki.atlassian.net/wiki/spaces/OPS/pages/13707660/Cartons+for+post+v1

    Existing carton code:
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_snc.py

    Simplified Description of selection criteria :
    Select everything within 100pc (i.e. parallax > 10mas) and Gaia-G
    fainter than 10th mag. Add some constraints on astrometric excess noise
    in crowded regions (Galactic plane, LMC, SMC, towards Galactic center).
    For details, see here. Pseudo-code: For the pseduo-code, I will shorten
    it and use the definition of "gal_cut" from the existing carton code.
    Also, the red text below is what is being added to the target selection
    compared to the original v1 carton.
    (parallax > 0) & (parallax - parallax_error > 10) & (phot_g_mean_mag > 10) &
    ((astrometric_excess_noise < 2) & (ll >= 20) & (ll <= 340) & gal_cut) |
    ((astrometric_excess_noise < 2) & (ruwe < 1.2) &
    ((1 + 0.015 * (phot_bp_mean_mag - phot_rp_mean_mag)^2) < phot_bp_rp_excess_factor) &
    (phot_bp_rp_excess_factor <
    (1.3 + 0.06 * (phot_bp_mean_mag - phot_rp_mean_mag)^2)) &
    ((ll < 20) | (ll > 340)) & gal_cut) | ~gal_cut)

    Return columns: Same columns as before + astrometric_excess_noise,
    phot_bp_rp_excess_factor, ruwe, phot_bp_mean_mag , phot_rp_mean_mag

    Metadata:
    Priority: 1800
    Cadence:  if G < 16, bright_2x1 ; if G > 16, dark_2x1
    Instrument: BOSS
    can_offset = True
    Lead contact: Ilija Medan

    """

    name = "mwm_snc_100pc_boss"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "BOSS"
    cadence = None  # cadence is set in post_process()
    priority = 1800
    can_offset = True

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

        query = super().build_query(version_id, query_region)
        query = query.where(
            Gaia_DR3.parallax > 0,
            Gaia_DR3.parallax - Gaia_DR3.parallax_error > 10,
            Gaia_DR3.phot_g_mean_mag > 10,
            ((Gaia_DR3.astrometric_excess_noise < 2) & (ll >= 20) & (ll <= 340) & gal_cut)
            | (
                (Gaia_DR3.astrometric_excess_noise < 2)
                & (Gaia_DR3.ruwe < 1.2)
                & (
                    (1 + 0.015 * fn.pow(Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag, 2))
                    < Gaia_DR3.phot_bp_rp_excess_factor
                )
                & (
                    (
                        Gaia_DR3.phot_bp_rp_excess_factor
                        < (
                            1.3
                            + 0.06
                            * fn.pow(Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag, 2)
                        )
                    )
                    & ((ll < 20) | (ll > 340))
                    & gal_cut
                )
                | ~gal_cut
            ),
        )

        return query

    def post_process(self, model, **kwargs):
        # G < 16 => cadence = bright_2x1
        model.update(cadence="bright_2x1").where(model.phot_g_mean_mag < 16).execute()

        # G > 16 => cadence = dark_2x1
        model.update(cadence="dark_2x1").where(model.phot_g_mean_mag >= 16).execute()

        return model


class MWM_SNC_100pc_BOSS_single_Carton(MWM_SNC_100pc_BOSS_Carton):
    """MWM SNC 100pc targets to be observed with BOSS.

    Shorthand name:  mwm_snc_100pc_boss_single
    https://sdss-wiki.atlassian.net/wiki/spaces/OPS/pages/13707660/Cartons+for+post+v1
    Simplified Description of selection criteria :
    NOTE this is a subcarton of mwm_snc_100pc_boss (see above) where the only change
    is the cadence and priority.

    Metadata:
    Priority: 1801
    Cadence:  if G < 16, bright_1x1 ; if G > 16, dark_1x1
    Instrument: BOSS
    can_offset = True
    Lead contact: Ilija Medan
    """

    name = "mwm_snc_100pc_boss_single"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "BOSS"
    cadence = None  # cadence is set in post_process()
    priority = 1801
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        return query

    def post_process(self, model, **kwargs):
        # G < 16 => cadence = bright_1x1
        model.update(cadence="bright_1x1").where(model.phot_g_mean_mag < 16).execute()

        # G > 16 => cadence = dark_1x1
        model.update(cadence="dark_1x1").where(model.phot_g_mean_mag >= 16).execute()

        return model


class MWM_SNC_Ext_Carton(BaseCarton):
    """MWM Solar Neighbourhood Census extension.

    Shorthand name: mwm_snc_ext_apogee
    Simplified Description of selection criteria:
        Above 100pc, apply a linear cut in absolute G-mag and distance to select increasingly
        more massive stars at larger distances. 2MASS H-band brighter than 11mag.
        Add some constraints on astrometric excess noise in crowded regions
        (Galactic plane, LMC, SMC).
    Pseudo-code:
        ((parallax<10) && (phot_g_mean_mag+5*log10(parallax/1000)+5) <
          19.125-0.03225*(1000/parallax))  &&
          ( (astrometric_excess_noise < 2 &&
            ( (l<=180 && b<-0.139*l + 25 && b>0.139*l - 25) ||
              (l>180 && b>-0.139*l + 25 && b <0.139*l - 25) ||
              (sqrt( pow(l-303.2,2) + 2*pow(b+44.4,2) ) < 5) ||
              (sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8))) ||
            !( (l<=180 && b<-0.139*l + 25 && b>0.139*l - 25) ||
               (l>180 && b>-0.139*l + 25 && b <0.139*l - 25) ||
               (sqrt( pow(l-303.2,2) + 2*pow(b+44.4,2) ) < 5) ||
               (sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8) ) ) && 2MASS_H < 11


    Shorthand name: mwm_snc_ext_boss
    Simplified Description of selection criteria:
        Above 100pc, apply a linear cut in absolute G-mag and distance to select increasingly
        more massive stars at larger distances. Gaia G-band fainter than 10mag.
        Add some constraints on astrometric excess noise in crowded regions
        (Galactic plane, LMC, SMC).
    Pseudo-code:
        ((parallax<10) && (phot_g_mean_mag+5*log10(parallax/1000)+5) <
          19.125-0.03225*(1000/parallax))  &&
          ( (astrometric_excess_noise < 2 &&
            ( (l<=180 && b<-0.139*l + 25 && b>0.139*l - 25) ||
              (l>180 && b>-0.139*l + 25 && b <0.139*l - 25) ||
              (sqrt( pow(l-303.2,2) + 2*pow(b+44.4,2) ) < 5) ||
              (sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8))) ||
            !( (l<=180 && b<-0.139*l + 25 && b>0.139*l - 25) ||
               (l>180 && b>-0.139*l + 25 && b <0.139*l - 25) ||
               (sqrt( pow(l-303.2,2) + 2*pow(b+44.4,2) ) < 5) ||
               (sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8) ) ) && phot_g_mag > 10

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

        mag_cut = (
            Gaia_DR3.phot_g_mean_mag + 5 * fn.log(Gaia_DR3.parallax / 1000) + 5
        ) < 19.125 - 0.03225 * (1000 / Gaia_DR3.parallax)

        query = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra,
                Gaia_DR3.dec,
                Gaia_DR3.phot_g_mean_mag,
                ll,
                bb,
            )
            .join(CatalogToGaia_DR3)
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                Gaia_DR3.parallax > 0,
                Gaia_DR3.parallax < 10,
                ((Gaia_DR3.astrometric_excess_noise < 2) & gal_cut) | ~gal_cut,
                mag_cut,
            )
        )

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


class MWM_SNC_Ext_APOGEE_Base_Carton(MWM_SNC_Ext_Carton):
    """SNC extension for APOGEE. See base carton for details.
    This is a base carton. The derived cartons are
    mwm_snc_ext_main_apogee: priority = 2705
    mwm_snc_ext_filler_apogee: priority = 7200
    """

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region=query_region)
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


class MWM_SNC_Ext_main_APOGEE_Carton(MWM_SNC_Ext_APOGEE_Base_Carton):
    """See base carton for details."""

    name = "mwm_snc_ext_main_apogee"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2705
    can_offset = True

    def post_process(self, model):
        """
        a set random seed, select 1/16 of the targets
        """
        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_snc_ext_main_apogee " + "set selected = false;"
        )

        # The below "order by catalogid" ensures that the random selection
        # further below gives the same result every time we run this carton.
        cursor = self.database.execute_sql(
            "select catalogid from "
            + " sandbox.temp_mwm_snc_ext_main_apogee "
            + " order by catalogid;"
        )

        output = cursor.fetchall()

        # This random seed must be the same as in mwm_snc_ext_filler_apogee
        random.seed(1234)
        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_random = random.randrange(16)
            if current_random == 0:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_main_apogee "
                    + " set selected = true "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )


class MWM_SNC_Ext_filler_APOGEE_Carton(MWM_SNC_Ext_APOGEE_Base_Carton):
    """See base carton for details."""

    name = "mwm_snc_ext_filler_apogee"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 7200
    can_offset = True

    def post_process(self, model):
        """
        a set random seed, select 15/16 of the targets
        """
        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_snc_ext_filler_apogee " + "set selected = false;"
        )

        # The below "order by catalogid" ensures that the random selection
        # further below gives the same result every time we run this carton.
        cursor = self.database.execute_sql(
            "select catalogid from "
            + " sandbox.temp_mwm_snc_ext_filler_apogee "
            + " order by catalogid;"
        )

        output = cursor.fetchall()

        # This random seed must be the same as in mwm_snc_ext_main_apogee
        random.seed(1234)
        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_random = random.randrange(16)
            if current_random != 0:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_filler_apogee "
                    + " set selected = true "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )


class MWM_SNC_Ext_BOSS_Base_Carton(MWM_SNC_Ext_Carton):
    """SNC extension for BOSS. See base carton for details.
    This is a base carton. The derived cartons are
    mwm_snc_ext_main_boss: priority = 2705
    mwm_snc_ext_filler_boss: priority = 7200"""

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag > 10)

        return query


class MWM_SNC_Ext_main_BOSS_Carton(MWM_SNC_Ext_BOSS_Base_Carton):
    """See base carton for details."""

    # cadence must be None here so that
    # it can be set in post_process().
    # If cadence is not None here then
    # it cannot be set in post_process().
    name = "mwm_snc_ext_main_boss"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "BOSS"
    cadence = None  # cadence is set in post_process()
    priority = 2705
    can_offset = True

    def post_process(self, model):
        """
        a set random seed, select 1/20 of the targets
        """

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_snc_ext_main_boss " + "set selected = false;"
        )

        # The below "order by catalogid" ensures that the random selection
        # further below gives the same result every time we run this carton.
        cursor = self.database.execute_sql(
            "select catalogid from "
            + " sandbox.temp_mwm_snc_ext_main_boss "
            + " order by catalogid;"
        )

        output = cursor.fetchall()

        # This random seed must be the same as in mwm_snc_ext_filler_boss
        random.seed(5678)
        for i in range(len(output)):
            current_catalogid = output[i][0]

            current_random = random.randrange(20)
            if current_random == 0:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_main_boss "
                    + " set selected = true "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )

            if Gaia_DR3.phot_g_mean_mag < 16:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_main_boss "
                    + " set cadence = 'bright_flexible_2x1' "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )
            else:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_main_boss "
                    + " set cadence = 'dark_flexible_2x1' "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )


class MWM_SNC_Ext_filler_BOSS_Carton(MWM_SNC_Ext_BOSS_Base_Carton):
    """See base carton for details."""

    # cadence must be None here so that
    # it can be set in post_process().
    # If cadence is not None here then
    # it cannot be set in post_process().
    name = "mwm_snc_ext_filler_boss"
    category = "science"
    program = "mwm_snc"
    mapper = "MWM"
    instrument = "BOSS"
    cadence = None  # cadence is set in post_process()
    priority = 7200
    can_offset = True

    def post_process(self, model):
        """
        a set random seed, select 19/20 of the targets
        """

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_snc_ext_filler_boss " + "set selected = false;"
        )

        # The below "order by catalogid" ensures that the random selection
        # further below gives the same result every time we run this carton.
        cursor = self.database.execute_sql(
            "select catalogid from "
            + " sandbox.temp_mwm_snc_ext_filler_boss "
            + " order by catalogid;"
        )

        output = cursor.fetchall()

        # This random seed must be the same as in mwm_snc_ext_main_boss
        random.seed(5678)
        for i in range(len(output)):
            current_catalogid = output[i][0]

            current_random = random.randrange(20)
            if current_random != 0:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_filler_boss "
                    + " set selected = true "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )

            if Gaia_DR3.phot_g_mean_mag < 16:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_filler_boss "
                    + " set cadence = 'bright_flexible_2x1' "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )
            else:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_filler_boss "
                    + " set cadence = 'dark_flexible_2x1' "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )
