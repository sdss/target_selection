#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-02
# @Filename: mwm_snc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import random

from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             CatalogToTwoMassPSC,
                                             Gaia_DR3, TwoMassPSC)

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
        gal_cut = (((ll <= 180) & (bb < (-0.139 * ll + 25)) & (bb > (0.139 * ll - 25)))
                   | ((ll > 180) & (bb > (-0.139 * ll + 25)) & (bb < (0.139 * ll - 25)))
                   | (fn.sqrt(fn.pow(ll - 303.2, 2) + 2 * fn.pow(bb + 44.4, 2)) < 5)
                   | (fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8))

        cte = Gaia_DR3.select(Gaia_DR3.source_id).where(
            Gaia_DR3.parallax > 0, (Gaia_DR3.parallax - Gaia_DR3.parallax_error)
            > 10).cte('plx_cte', materialized=True)

        query = (Gaia_DR3.select(
            CatalogToGaia_DR3.catalogid,
            Gaia_DR3.source_id,
            Gaia_DR3.phot_g_mean_mag,
            ll,
            bb,
        ).join(CatalogToGaia_DR3).join_from(
            Gaia_DR3, cte, on=(Gaia_DR3.source_id == cte.c.source_id)).where((
                (Gaia_DR3.astrometric_excess_noise < 2) & gal_cut) | ~gal_cut).where(
                    CatalogToGaia_DR3.version_id == version_id,
                    CatalogToGaia_DR3.best >> True).with_cte(cte))

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                ))

        return query


class MWM_SNC_100pc_APOGEE_Carton(MWM_SNC_100pc_Carton):
    """MWM SNC 100pc targets to be observed with APOGEE."""

    name = 'mwm_snc_100pc_apogee'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 1805
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region=query_region)
        query = (query.switch(CatalogToGaia_DR3).join(
            CatalogToTwoMassPSC, on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid
                                     )).join(TwoMassPSC).where(
                                         CatalogToTwoMassPSC.version_id == version_id,
                                         CatalogToTwoMassPSC.best >> True,
                                         TwoMassPSC.h_m < 11))

        return query


class MWM_SNC_100pc_BOSS_Carton(MWM_SNC_100pc_Carton):
    """MWM SNC 100pc targets to be observed with BOSS."""

    name = 'mwm_snc_100pc_boss'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'BOSS'
    priority = 1800
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag > 10)

        return query

    def post_process(self, model, **kwargs):

        # G > 16 => cadence = dark_2x1
        model.update(cadence='dark_2x1').where(model.phot_g_mean_mag > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence='bright_2x1').where(model.phot_g_mean_mag < 16).execute()

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
        gal_cut = (((ll <= 180) & (bb < (-0.139 * ll + 25)) & (bb > (0.139 * ll - 25)))
                   | ((ll > 180) & (bb > (-0.139 * ll + 25)) & (bb < (0.139 * ll - 25)))
                   | (fn.sqrt(fn.pow(ll - 303.2, 2) + 2 * fn.pow(bb + 44.4, 2)) < 5)
                   | (fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8))

        mag_cut = ((Gaia_DR3.phot_g_mean_mag + 5 * fn.log(Gaia_DR3.parallax / 1000) + 5)
                   < 19.125 - 0.03225 * (1000 / Gaia_DR3.parallax))

        query = (Gaia_DR3.select(
            CatalogToGaia_DR3.catalogid,
            Gaia_DR3.source_id,
            Gaia_DR3.phot_g_mean_mag,
            ll,
            bb,
        ).join(CatalogToGaia_DR3).where(
            CatalogToGaia_DR3.version_id == version_id, CatalogToGaia_DR3.best >> True,
            Gaia_DR3.parallax > 0, Gaia_DR3.parallax < 10,
            ((Gaia_DR3.astrometric_excess_noise < 2) & gal_cut) | ~gal_cut, mag_cut))

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                ))

        return query


class MWM_SNC_Ext_APOGEE_Carton(MWM_SNC_Ext_Carton):
    """SNC extension for APOGEE. See base carton for details."""

    # priority must be None here so that
    # it can be set in post_process().
    # If priority is not None here then it cannot be set in post_process().
    name = 'mwm_snc_ext_apogee'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = None    # priority is set in post_process()
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region=query_region)
        query = (query.switch(CatalogToGaia_DR3).join(
            CatalogToTwoMassPSC,
            on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid
                )).join(TwoMassPSC).where(CatalogToTwoMassPSC.version_id == version_id,
                                          CatalogToTwoMassPSC.best >> True, TwoMassPSC.h_m < 11))

        return query

    def post_process(self, model):
        """
        a set random seed, put 1/4 of them at priority 2705 and
         the rest at priority 4900.
        """

        cursor = self.database.execute_sql(
            "select catalogid from " +
            " sandbox.temp_mwm_snc_ext_apogee " +
            " order by catalogid;")

        output = cursor.fetchall()

        random.seed(1234)
        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_random = random.randrange(4)
            if (current_random == 0):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_apogee set priority = 2705 " +
                    " where catalogid = " + str(current_catalogid) + ";")
            else:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_apogee set priority = 4900 " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_SNC_Ext_BOSS_Carton(MWM_SNC_Ext_Carton):
    """SNC extension for BOSS. See base carton for details."""

    # cadence and priority must be None here so that
    # they can be set in post_process().
    # If cadence and priority are not None here then
    # they cannot be set in post_process().
    name = 'mwm_snc_ext_boss'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'BOSS'
    cadence = None    # cadence is set in post_process()
    priority = None    # priority is set in post_process()
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag > 10)

        return query

    def post_process(self, model):
        """
        a set random seed, put 1/4 of them at priority 2705 and
         the rest at priority 4900.
        """

        cursor = self.database.execute_sql(
            "select catalogid from " +
            " sandbox.temp_mwm_snc_ext_boss " +
            " order by catalogid;")

        output = cursor.fetchall()

        random.seed(5678)
        for i in range(len(output)):
            current_catalogid = output[i][0]

            current_random = random.randrange(4)
            if (current_random == 0):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_boss set priority = 2705 " +
                    " where catalogid = " + str(current_catalogid) + ";")
            else:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_boss set priority = 4900 " +
                    " where catalogid = " + str(current_catalogid) + ";")

            if (Gaia_DR3.phot_g_mean_mag < 16):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_boss set cadence = 'bright_2x1' " +
                    " where catalogid = " + str(current_catalogid) + ";")
            else:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_snc_ext_boss set cadence = 'dark_2x1' " +
                    " where catalogid = " + str(current_catalogid) + ";")
