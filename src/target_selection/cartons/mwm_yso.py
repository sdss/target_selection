#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    MIPSGAL,
    AllWise,
    Catalog,
    CatalogToAllWise,
    CatalogToGaia_DR3,
    CatalogToTIC_v8,
    CatalogToTwoMassPSC,
    Gaia_DR2,
    Gaia_DR3,
    Sagitta_EDR3,
    TIC_v8,
    TwoMassPSC,
    YSO_Clustering,
    Zari18pms,
)

from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class MWM_YSO_Disk_APOGEE_Carton(BaseCarton):
    """YSOs - Disk APOGEE (IR excess).

    Shorthand name: mwm_yso_disk_apogee

    old class name: MWM_YSO_S1_Carton
    old shorthand name: mwm_yso_s1

    Simplified Description of selection criteria:
    selection of YSOs based on IR excess,
    with WISE colors W1-W2>0.25, W2-W3>0.5, W3-W4>1.5,
    closer than parallax>0.3, and brighter than H<13
    (should have ~21.5K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: Gaia, 2mass, allwise
    Additional cross-matching needed:
    Note: Using the Gaia xmatch somehow misses half the sources.
    Selection was done on the allwise catalog that
    had 2mass photometry,
    and then the resulting selection was crossmatched against against
    Gaia with 1" search radius.
    Return columns: Gaia id, 2mass id, allwise id, G, BP, RP,
    J, H, K, W1, W2, W3, W4,parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: h_m<13 and w1mpro-w2mpro>0.25 and
    w2mpro-w3mpro>0.5 and w3mpro-w4mpro>1.5 and parallax>0.3

    """

    name = "mwm_yso_disk_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                AllWise.designation.alias("allwise_designation"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag.alias("gaia_dr3_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR3.parallax,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .switch(CatalogToGaia_DR3)
            .join(CatalogToAllWise, on=(CatalogToGaia_DR3.catalogid == CatalogToAllWise.catalogid))
            .join(AllWise, on=(CatalogToAllWise.target_id == AllWise.cntr))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.best >> True,
                CatalogToAllWise.best >> True,
                TwoMassPSC.h_m < 13,
                (AllWise.w1mpro - AllWise.w2mpro) > 0.25,
                (AllWise.w2mpro - AllWise.w3mpro) > 0.50,
                (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                Gaia_DR3.parallax > 0.3,
            )
        )

        # There can be cases in which the same catalogid has multiple entries
        # in a catalog_to_x table since the same physical object
        # may match with multiple catalogids.
        # Hence, we have the below code in the above query.
        #                CatalogToGaia_DR3.best >> True,
        #                CatalogToTwoMassPSC.best >> True,
        #                CatalogToAllWise.best >> True,

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


class MWM_YSO_Disk_APOGEE_Single_Carton(MWM_YSO_Disk_APOGEE_Carton):
    name = "mwm_yso_disk_apogee_single"
    cadence = "bright_1x1"
    priority = 2706


class MWM_YSO_Disk_BOSS_Carton(BaseCarton):
    """YSOs - Disk BOSS (IR excess).

    Shorthand name: mwm_yso_disk_boss

    old class name: MWM_YSO_S1_Carton
    old shorthand name: mwm_yso_s1

    Simplified Description of selection criteria:
    selection of YSOs based on IR excess,
    with WISE colors W1-W2>0.25, W2-W3>0.5, W3-W4>1.5,
    closer than parallax>0.3, and brighter than H<13
    (should have ~21.5K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: Gaia, 2mass, allwise
    Additional cross-matching needed:
    Note: Using the Gaia xmatch somehow misses half the sources.
    Selection was done on the allwise catalog that
    had 2mass photometry,
    and then the resulting selection was crossmatched against against
    Gaia with 1" search radius.
    Return columns: Gaia id, 2mass id, allwise id, G, BP, RP,
    J, H, K, W1, W2, W3, W4,parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    boss_bright_3x1 if RP<14.76 |
    boss_bright_4x1 if RP<15.075 |
    boss_bright_5x1 if RP<15.29 |
    boss_bright_6x1 if RP<15.5
    Pseudo SQL (optional):
    Implementation: phot_rp_mean_mag<15.5 and w1mpro-w2mpro>0.25 and
    w2mpro-w3mpro>0.5 and w3mpro-w4mpro>1.5 and parallax>0.3

    Comments: Split from mwm_yso_s1 to request BOSS observations,
    same color selection but assigning cadence and faint limit for carton based
    on RP instead of H
    """

    name = "mwm_yso_disk_boss"
    category = "science"
    instrument = None  # instrument is set in post_process()
    cadence = None  # cadence is set in post_process()
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                AllWise.designation.alias("allwise_designation"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag.alias("gaia_dr3_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR3.parallax,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .switch(CatalogToGaia_DR3)
            .join(CatalogToAllWise, on=(CatalogToGaia_DR3.catalogid == CatalogToAllWise.catalogid))
            .join(AllWise, on=(CatalogToAllWise.target_id == AllWise.cntr))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.best >> True,
                CatalogToAllWise.best >> True,
                Gaia_DR3.phot_rp_mean_mag < 15.5,
                (AllWise.w1mpro - AllWise.w2mpro) > 0.25,
                (AllWise.w2mpro - AllWise.w3mpro) > 0.50,
                (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                Gaia_DR3.parallax > 0.3,
            )
        )

        # Gaia_DR3 pweewee model class corresponds to
        # table catalogdb.gaia_dr3_source.
        #

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

    def post_process(self, model):
        """
        cadence options for these targets:
        boss_bright_3x1 if RP<14.76 |
        boss_bright_4x1 if RP<15.075 |
        boss_bright_5x1 if RP<15.29 |
        boss_bright_6x1 if RP<15.5
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr3_rp from " + " sandbox.temp_mwm_yso_disk_boss ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_rp = output[i][1]

            if current_rp < 14.76:
                current_instrument = "BOSS"
                current_cadence = "bright_3x1"
            elif current_rp < 15.075:
                current_instrument = "BOSS"
                current_cadence = "bright_4x1"
            elif current_rp < 15.29:
                current_instrument = "BOSS"
                current_cadence = "bright_5x1"
            elif current_rp < 15.5:
                current_instrument = "BOSS"
                current_cadence = "bright_6x1"
            else:
                # All cases should be covered above so we should not get here.
                current_instrument = None
                current_cadence = None
                raise TargetSelectionError(
                    "error in mwm_yso_disk_boss "
                    + "post_process(): "
                    + "instrument = None, cadence= None"
                )

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_disk_boss "
                    + " set instrument = '"
                    + current_instrument
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_disk_boss "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )


class MWM_YSO_Disk_BOSS_Single_Carton(MWM_YSO_Disk_BOSS_Carton):
    name = "mwm_yso_disk_boss_single"
    cadence = "bright_1x1"
    priority = 2706
    instrument = "BOSS"

    def post_process(self, model):
        pass


# The below class is from commit of Feb 2, 2022.
# The only changes are
# (1) can_offset = True
# (2) addition of the
# last condition in the WHERE clause and
# (3) addition of TwoMassPSC.ra and TwoMassPSC.decl
# in the SELECT clause.
class MWM_YSO_Embedded_APOGEE_Carton(BaseCarton):
    """YSOs - Embedded APOGEE (optically invisible).

    Shorthand name: mwm_yso_embedded_apogee

    old class name: MWM_YSO_S2_Carton
    old shorthand name: mwm_yso_s2

    Wiki page: v0.5
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass+allwise, gaia
    (allow sources that lack gaia xmatch)
    Additional cross-matching needed:
    Note: Using the Gaia xmatch somehow misses half the sources.
    Selection was done on the allwise catalog
    that had 2mass photometry,
    and then the resulting selection was crossmatched
    against against Gaia with 1" search radius.
    Return columns: Gaia id, 2mass id, allwise id, G, BP, RP,
    J, H, K, W1, W2, W3, W4
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):

    Pseudo SQL (optional): v1.0
    Implementation: h_m<13 and
    (phot_g_mean_mag>18.5 or phot_g_mean_mag is null)
    and j_m-h_m>1
    and h_m-ks_m>0.5 << v0.5 and v1.0 condition on H-K
    and w1mpro-w2mpro>0.5
    and w2mpro-w3mpro>1
    and w3mpro-w4mpro>1.5
    and w3mpro-w4mpro>(w1mpro-w2mpro)*0.8+1.1
    and (H-K)>0.65*(J-H)-0.25 <<< This is additional condition on H-K
                                  compared to v0.5
    """

    name = "mwm_yso_embedded_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            AllWise.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                TwoMassPSC.ra.alias("twomass_psc_ra"),
                TwoMassPSC.decl.alias("twomass_psc_decl"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                AllWise.designation.alias("allwise_designation"),
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag.alias("gaia_dr2_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR2.parallax,
            )
            .join(TIC_v8, on=(TIC_v8.allwise == AllWise.designation))
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .switch(TIC_v8)
            .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .switch(TIC_v8)
            .join(CatalogToTIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                TwoMassPSC.h_m < 13,
                (Gaia_DR2.phot_g_mean_mag > 18.5) | (Gaia_DR2.phot_g_mean_mag >> None),
                (AllWise.j_m_2mass - AllWise.h_m_2mass) > 1.0,
                (AllWise.h_m_2mass - AllWise.k_m_2mass) > 0.5,
                (AllWise.w1mpro - AllWise.w2mpro) > 0.50,
                (AllWise.w2mpro - AllWise.w3mpro) > 1.00,
                (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                (AllWise.w3mpro - AllWise.w4mpro) > (AllWise.w1mpro - AllWise.w2mpro) * 0.8 + 1.1,
                (
                    (AllWise.h_m_2mass - AllWise.k_m_2mass)  # added for v1.0
                    > (0.65 * (AllWise.j_m_2mass - AllWise.h_m_2mass) - 0.25)
                ),
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


class MWM_YSO_Embedded_APOGEE_Single_Carton(MWM_YSO_Embedded_APOGEE_Carton):
    name = "mwm_yso_embedded_apogee_single"
    cadence = "bright_1x1"
    priority = 2706


class MWM_YSO_Nebula_APOGEE_Carton(BaseCarton):
    """YSOs - Nebula APOGEE(optically invisible, WISE saturated).

    Shorthand name: mwm_yso_nebula_apogee

    old class name: MWM_YSO_S2_5_Carton
    old shorthand name: mwm_yso_s2_5

    Simplified Description of selection criteria:
    selection of YSOs, brighter than H<15,
    saturated (blank) W4 with W2-W3>4,
    or saturated W3 and W2, with J-H>1.1.
    Some contaminants from scanning are
    filtered on the plane of the sky:
    all the targets should be within 5 deg of the plane+
    few sources that can be located
    further south of the plane if l>180
    (should have ~1.2K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass, allwise
    Additional cross-matching needed:
    Return columns: 2mass id, allwise id, J, H, K, W1, W2, W3, W4
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: h_m<13 and
    (w2mpro-w3mpro>4 and w4mpro is null) or
    (w3mpro is null and w4mpro is null and j_m-h_m>1.1)
    and (b>-5 or l>180) and b<-5

    """

    name = "mwm_yso_nebula_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    # Above implementation has below clause
    # and (b>-5 or l>180) and b<-5
    # Replace (b>-5 or l>180) and b<-5 as below based on the text.
    # In words:
    # all the targets should be within 5 deg of the plane+
    # few sources that can be
    # located further south of the plane if l>180
    # Hence:
    # ((b>-5) and (b<5)) or ((b<-5) and (l > 180))
    #  l, b in Gaia_DR2 are gallong and gallat in TIC_v8.
    # We are using the values from Gaia since
    # TIC propagates the coordinates back to epoch 2000.0
    # (b>-5 or l>180) and b<-5
    # S2_5 query below has the same part before where() as S2 query.

    def build_query(self, version_id, query_region=None):
        query = (
            AllWise.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                TwoMassPSC.ra.alias("twomass_psc_ra"),
                TwoMassPSC.decl.alias("twomass_psc_decl"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                AllWise.designation.alias("allwise_designation"),
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag.alias("gaia_dr2_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR2.parallax,
            )
            .join(TIC_v8, on=(TIC_v8.allwise == AllWise.designation))
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .switch(TIC_v8)
            .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .switch(TIC_v8)
            .join(CatalogToTIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                TwoMassPSC.h_m < 13,
                (((AllWise.w2mpro - AllWise.w3mpro) > 4) & (AllWise.w4mpro >> None))
                | (
                    (AllWise.w3mpro >> None)
                    & (AllWise.w4mpro >> None)
                    & ((AllWise.j_m_2mass - AllWise.h_m_2mass) > 1.1)
                ),
                ((Gaia_DR2.b > -5) & (Gaia_DR2.b < 5))
                | ((Gaia_DR2.b < -5) & (Gaia_DR2.l > 180))
                | ((Gaia_DR2.b >> None) & (Gaia_DR2.l >> None)),
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


class MWM_YSO_Nebula_APOGEE_Single_Carton(MWM_YSO_Nebula_APOGEE_Carton):
    name = "mwm_yso_nebula_apogee_single"
    cadence = "bright_1x1"
    priority = 2706


class MWM_YSO_Variable_APOGEE_Carton(BaseCarton):
    """YSOs - Variable APOGEE (pre-main sequence optical variables).

    Shorthand name: mwm_yso_variable_apogee

    old class name: MWM_YSO_S3_Carton
    old shorthand name: mwm_yso_s3

    Simplified Description of selection criteria:
    selection of YSOs brighter than H<13, closer than parallax>0.3.
    Filter on the position of the HR diagram to
    select cool pre-main sequence stars,
    with BP-RP>13, (BP-RP)*2.5+2.5>M_G, (BP-RP)*2.5-1<M_G,
    requiring variability in g,bp,rp>0.02
    (with var_x defined as
    sqrt(phot_x_n_obs)/phot_x_mean_flux_over_error),
    have relations in variability of
    var_g<var_bp<var_g^0.75, 0.75*var_g<var_rp<var_g^0.95,
    and log10(var_bp)*5+11<M_BP, in which M_x is the absolute mag
    (should have ~52.7K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass, gaia
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation:
    phot_g_mean_mag < 18.5 and h_m <13 and parallax >0.3 and
    bp_rp*2.5+2.5 > phot_g_mean_mag-5*(log10(1000/parallax)-1) and
    bp_rp*2.5-1 < phot_g_mean_mag-5*(log10(1000/parallax)-1) and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error>
    sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error>
    sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error*0.75 and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error<
    power(sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error,0.75) and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error<
    power(sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error,0.95) and
    log10(sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error)*5+11<
    phot_bp_mean_mag-5*(log10(1000/parallax)-1) and
    bp_rp>1.3 and sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error>0.02
    and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error>0.02 and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error>0.02

    """

    name = "mwm_yso_variable_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                AllWise.designation.alias("allwise_designation"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag.alias("gaia_dr3_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR3.parallax,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .switch(CatalogToGaia_DR3)
            .join(CatalogToAllWise, on=(CatalogToGaia_DR3.catalogid == CatalogToAllWise.catalogid))
            .join(AllWise, on=(CatalogToAllWise.target_id == AllWise.cntr))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.best >> True,
                CatalogToAllWise.best >> True,
                Gaia_DR3.phot_g_mean_mag < 18.5,
                TwoMassPSC.h_m < 13,
                Gaia_DR3.parallax > 0.3,
                Gaia_DR3.bp_rp * 2.5 + 2.5
                > Gaia_DR3.phot_g_mean_mag - 5 * (peewee.fn.log(1000 / Gaia_DR3.parallax) - 1),
                Gaia_DR3.bp_rp * 2.5 - 1
                < Gaia_DR3.phot_g_mean_mag - 5 * (peewee.fn.log(1000 / Gaia_DR3.parallax) - 1),
                peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                > peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error,
                peewee.fn.sqrt(Gaia_DR3.phot_rp_n_obs) / Gaia_DR3.phot_rp_mean_flux_over_error
                > peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs)
                / Gaia_DR3.phot_g_mean_flux_over_error
                * 0.75,
                peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                < peewee.fn.power(
                    peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error,
                    0.75,
                ),
                peewee.fn.sqrt(Gaia_DR3.phot_rp_n_obs) / Gaia_DR3.phot_rp_mean_flux_over_error
                < peewee.fn.power(
                    peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error,
                    0.95,
                ),
                peewee.fn.log(
                    peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                )
                * 5
                + 11
                < Gaia_DR3.phot_bp_mean_mag - 5 * (peewee.fn.log(1000 / Gaia_DR3.parallax) - 1),
                Gaia_DR3.bp_rp > 1.3,
                peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error
                > 0.02,
                peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                > 0.02,
                peewee.fn.sqrt(Gaia_DR3.phot_rp_n_obs) / Gaia_DR3.phot_rp_mean_flux_over_error
                > 0.02,
            )
        )

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


class MWM_YSO_Variable_APOGEE_Single_Carton(MWM_YSO_Variable_APOGEE_Carton):
    name = "mwm_yso_variable_apogee_single"
    cadence = "bright_1x1"
    priority = 2706


class MWM_YSO_Variable_BOSS_Carton(BaseCarton):
    """YSOs - Variable BOSS (pre-main sequence optical variables).

    Shorthand name: mwm_yso_variable_boss

    old class name: MWM_YSO_S3_Carton
    old shorthand name: mwm_yso_s3

    Simplified Description of selection criteria:
    selection of YSOs brighter than H<13, closer than parallax>0.3.
    Filter on the position of the HR diagram to
    select cool pre-main sequence stars,
    with BP-RP>13, (BP-RP)*2.5+2.5>M_G, (BP-RP)*2.5-1<M_G,
    requiring variability in g,bp,rp>0.02
    (with var_x defined as
    sqrt(phot_x_n_obs)/phot_x_mean_flux_over_error),
    have relations in variability of
    var_g<var_bp<var_g^0.75, 0.75*var_g<var_rp<var_g^0.95,
    and log10(var_bp)*5+11<M_BP, in which M_x is the absolute mag
    (should have ~52.7K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass, gaia
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    boss_bright_3x1 if RP<14.76 |
    boss_bright_4x1 if RP<15.075 |
    boss_bright_5x1 if RP<15.29 |
    boss_bright_6x1 if RP<15.5
    Pseudo SQL (optional):
    Implementation:
    phot_rp_mean_mag<15.5 and phot_g_mean_mag < 18.5 and h_m <13 and parallax >0.3 and
    bp_rp*2.5+2.5 > phot_g_mean_mag-5*(log10(1000/parallax)-1) and
    bp_rp*2.5-1 < phot_g_mean_mag-5*(log10(1000/parallax)-1) and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error>
    sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error>
    sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error*0.75 and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error<
    power(sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error,0.75) and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error<
    power(sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error,0.95) and
    log10(sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error)*5+11<
    phot_bp_mean_mag-5*(log10(1000/parallax)-1) and
    bp_rp>1.3 and sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error>0.02
    and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error>0.02 and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error>0.02

    Comments: Split from mwm_yso_s3 to request BOSS observations,
    RP magnitude check added to the previous selection
    """

    name = "mwm_yso_variable_boss"
    category = "science"
    instrument = None  # instrument is set in post_process()
    cadence = None  # cadence is set in post_process()
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                AllWise.designation.alias("allwise_designation"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag.alias("gaia_dr3_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR3.parallax,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .switch(CatalogToGaia_DR3)
            .join(CatalogToAllWise, on=(CatalogToGaia_DR3.catalogid == CatalogToAllWise.catalogid))
            .join(AllWise, on=(CatalogToAllWise.target_id == AllWise.cntr))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.best >> True,
                CatalogToAllWise.best >> True,
                Gaia_DR3.phot_rp_mean_mag < 15.5,
                Gaia_DR3.phot_g_mean_mag < 18.5,
                TwoMassPSC.h_m < 13,
                Gaia_DR3.parallax > 0.3,
                Gaia_DR3.bp_rp * 2.5 + 2.5
                > Gaia_DR3.phot_g_mean_mag - 5 * (peewee.fn.log(1000 / Gaia_DR3.parallax) - 1),
                Gaia_DR3.bp_rp * 2.5 - 1
                < Gaia_DR3.phot_g_mean_mag - 5 * (peewee.fn.log(1000 / Gaia_DR3.parallax) - 1),
                peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                > peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error,
                peewee.fn.sqrt(Gaia_DR3.phot_rp_n_obs) / Gaia_DR3.phot_rp_mean_flux_over_error
                > peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs)
                / Gaia_DR3.phot_g_mean_flux_over_error
                * 0.75,
                peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                < peewee.fn.power(
                    peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error,
                    0.75,
                ),
                peewee.fn.sqrt(Gaia_DR3.phot_rp_n_obs) / Gaia_DR3.phot_rp_mean_flux_over_error
                < peewee.fn.power(
                    peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error,
                    0.95,
                ),
                peewee.fn.log(
                    peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                )
                * 5
                + 11
                < Gaia_DR3.phot_bp_mean_mag - 5 * (peewee.fn.log(1000 / Gaia_DR3.parallax) - 1),
                Gaia_DR3.bp_rp > 1.3,
                peewee.fn.sqrt(Gaia_DR3.phot_g_n_obs) / Gaia_DR3.phot_g_mean_flux_over_error
                > 0.02,
                peewee.fn.sqrt(Gaia_DR3.phot_bp_n_obs) / Gaia_DR3.phot_bp_mean_flux_over_error
                > 0.02,
                peewee.fn.sqrt(Gaia_DR3.phot_rp_n_obs) / Gaia_DR3.phot_rp_mean_flux_over_error
                > 0.02,
            )
        )

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

    def post_process(self, model):
        """
        cadence options for these targets:
        boss_bright_3x1 if RP<14.76 |
        boss_bright_4x1 if RP<15.075 |
        boss_bright_5x1 if RP<15.29 |
        boss_bright_6x1 if RP<15.5
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr3_rp from " + " sandbox.temp_mwm_yso_variable_boss ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_rp = output[i][1]

            if current_rp < 14.76:
                current_instrument = "BOSS"
                current_cadence = "bright_3x1"
            elif current_rp < 15.075:
                current_instrument = "BOSS"
                current_cadence = "bright_4x1"
            elif current_rp < 15.29:
                current_instrument = "BOSS"
                current_cadence = "bright_5x1"
            elif current_rp < 15.5:
                current_instrument = "BOSS"
                current_cadence = "bright_6x1"
            else:
                # All cases should be covered above so we should not get here.
                current_instrument = None
                current_cadence = None
                raise TargetSelectionError(
                    "error in mwm_yso_variable_boss "
                    + "post_process(): "
                    + "instrument = None, cadence= None"
                )

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_variable_boss "
                    + " set instrument = '"
                    + current_instrument
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_variable_boss "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )


class MWM_YSO_Variable_BOSS_Single_Carton(MWM_YSO_Variable_BOSS_Carton):
    name = "mwm_yso_variable_boss_single"
    cadence = "bright_1x1"
    priority = 2706
    instrument = "BOSS"

    def post_process(self, model):
        pass


class MWM_YSO_CMZ_APOGEE_Carton(BaseCarton):
    """YSOs - Central Molecular Zone APOGEE.

    Shorthand name: mwm_yso_cmz_apogee

    old class name: MWM_YSO_CMZ_Carton
    old shorthand name: mwm_yso_cmz

    Simplified Description of selection criteria:
    selection of sources in the central molecular zone
    based on spitzer fluxes from mipsgal.
     brighter than H<13, have color 8.0-24>2.5, and
     have parallax<0.2 or lack a Gaia xmatch.
    (should have ~3.2K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: mipsgal
    Additional cross-matching needed: the table has xmatch included
    Return columns:
    mipsgal id, 2mass id, j, h, k, 3.6, 4.8, 8.0, 24 mag
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    'apogee_bright_3x1'
    Pseudo SQL (optional):
    Implementation: Hmag<13 and _8_0_-_24_>2.5 and
    (parallax<0.2 or parallax is null)

    DONE (email of 1/6/20232 with jsg) check left outer join
    TODO Below comments are for v0.5. modify for v1.0.

    For CMZ, the raw sql query would be:
    select ct.catalogid from mipsgal m
    join twomass_psc t on twomass_name = designation
    join tic_v8 tic on tic.twomass_psc = t.designation
    left outer join gaia_dr2_source g on g.source_id = tic.gaia_int
    join catalog_to_tic_v8 ct on ct.target_id = tic.id
    where m.hmag < 13 and
    (m.mag_8_0 - m.mag_24) > 2.5 and
    (g.parallax < 0.2 or g.parallax is null)
    and ct.version_id = 13 and ct.best is true;

    Note you only need one left outer join between TIC and Gaia
    (all MIPSGAL targets have a counterpart in 2MASS,
    and all 2MASS have an entry in TIC,
    but not all the TIC entries have a Gaia counterpart).

    Comments: Formerly mwm_yso_cmz, removed check on the position on the sky:

    Removed below condition.
    l is glon (galactic longitude)
    b is glat (galactic latitude)
    All four statements below are equivalent.
    (l> 358 or l< 2) and
    b between -1 and 1

    (m.glon > 358 or m.glon < 2) and
    (m.glat > -1 and m.glat < 1) and

    Sources are within 2 degrees in l and
     1 degree in b from the galactic center,

    (MIPSGAL.glon > 358) | (MIPSGAL.glon < 2),
    (MIPSGAL.glat > -1) & (MIPSGAL.glat < 1),

    """

    name = "mwm_yso_cmz_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    # mipsgal is a subset of 2MASS
    # mipsgal can be joined to twomass_psc via
    # mipsgal.twomass_name = TwoMassPSC.designation.
    #
    # mipsgal is a subset of 2MASS
    #
    # 2MASS is not a subset of Gaia_DR3
    # Gaia_DR3 is not a subset of 2MASS
    #
    # table catalogdb.mipsgal
    # Foreign-key constraints:
    #    "twomass_name_fk" FOREIGN KEY (twomass_name)
    # REFERENCES twomass_psc(designation)
    #
    # Since the below condition is true,
    # we do not need a between to Catalog and CatalogToTwoMassPSC
    #
    # Catalog.catalogid == CatalogToTwoMassPSC.catalogid
    #
    # We can remove the join with Catalog in all the cartons
    # since catalogid is completely unique (even across different version_id)
    # so the join with Catalog doesn't give us anything extra and it's a costly join.

    def build_query(self, version_id, query_region=None):
        # after the first left outer join, the later joins must
        # also be left outer joins
        query = (
            CatalogToTwoMassPSC.select(
                CatalogToTwoMassPSC.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                MIPSGAL.mag_3_6,
                MIPSGAL.mag_4_5,
                MIPSGAL.mag_5_8,
                MIPSGAL.mag_8_0,
                MIPSGAL.mag_24,
                MIPSGAL.hmag,
                Gaia_DR3.parallax,
                MIPSGAL.glon,
                MIPSGAL.glat,
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .join(MIPSGAL, on=(MIPSGAL.twomass_name == TwoMassPSC.designation))
            .switch(CatalogToTwoMassPSC)
            .join(
                CatalogToGaia_DR3,
                peewee.JOIN.LEFT_OUTER,
                on=(CatalogToTwoMassPSC.catalogid == CatalogToGaia_DR3.catalogid),
            )
            .join(
                Gaia_DR3,
                peewee.JOIN.LEFT_OUTER,
                on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id),
            )
            .distinct(CatalogToTwoMassPSC.catalogid)
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,  # See below for CatalogToGaia_DR3.best
                MIPSGAL.hmag < 13,
                (MIPSGAL.mag_8_0 - MIPSGAL.mag_24) > 2.5,
                ((Gaia_DR3.parallax < 0.2) & (CatalogToGaia_DR3.best >> True))
                | (Gaia_DR3.parallax >> None),
            )
        )
        # above condition (Gaia_DR3.parallax >> None) ensures that
        # that we get the rows from the left outer join

        if query_region:
            query = query.join_from(CatalogToTwoMassPSC, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


class MWM_YSO_CMZ_APOGEE_Single_Carton(MWM_YSO_CMZ_APOGEE_Carton):
    name = "mwm_yso_cmz_apogee_single"
    cadence = "bright_1x1"
    priority = 2706


class MWM_YSO_Cluster_APOGEE_Carton(BaseCarton):
    """YSOs - Cluster APOGEE Catalog
    Shorthand name: mwm_yso_cluster_apogee

    old class name: MWM_YSO_Cluster_Carton
    old shorthand name: mwm_yso_cluster

    Simplified Description of selection criteria:
    Selecting the clustered sources from
    the catalog of clustered structures,
    with age<7.5 dex and brighter than H<13 mag.
    (should have ~45.5K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: Kounkel+20 clustered catalog
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: age<7.5 and h<13

    """

    name = "mwm_yso_cluster_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    # yso_clustering is a subset of gaia and
    # can be joined to gaia_dr2_source via source_id.
    #
    # table catalogdb.yso_clustering
    # Foreign-key constraints:
    #    "yso_clustering_source_id_fkey" FOREIGN KEY (source_id)
    # REFERENCES gaia_dr2_source(source_id)

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                YSO_Clustering.twomass,
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag.alias("gaia_dr2_rp"),
                YSO_Clustering.j,
                YSO_Clustering.h,
                YSO_Clustering.k,
                Gaia_DR2.parallax,
            )
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .join(YSO_Clustering, on=(Gaia_DR2.source_id == YSO_Clustering.source_id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                YSO_Clustering.h < 13,
                YSO_Clustering.age < 7.5,
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


class MWM_YSO_Cluster_APOGEE_Single_Carton(MWM_YSO_Cluster_APOGEE_Carton):
    name = "mwm_yso_cluster_apogee_single"
    cadence = "bright_1x1"
    priority = 2706


class MWM_YSO_Cluster_BOSS_Carton(BaseCarton):
    """YSOs - Cluster BOSS Catalog
    Shorthand name: mwm_yso_cluster_boss

    old class name: MWM_YSO_Cluster_Carton
    old shorthand name: mwm_yso_cluster

    Simplified Description of selection criteria:
    Selecting the clustered sources from
    the catalog of clustered structures,
    with age<7.5 dex and brighter than rp<15.5 mag.

    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: Kounkel+20 clustered catalog
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    cadence options for these targets:
    boss_bright_3x1 if RP<14.76 |
    boss_bright_4x1 if RP<15.075 |
    boss_bright_5x1 if RP<15.29 |
    boss_bright_6x1 if RP<15.5
    Pseudo SQL (optional):
    Implementation: age<7.5 and rp<15.5

    Comments: Split from Cluster to request BOSS observations,
    assigning cadence and faint limit for carton based on RP instead of H
    """

    name = "mwm_yso_cluster_boss"
    category = "science"
    instrument = None  # instrument is set in post_process()
    cadence = None  # cadence is set in post_process()
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2705
    can_offset = True

    # yso_clustering is a subset of gaia and
    # can be joined to gaia_dr2_source via source_id.
    #
    # table catalogdb.yso_clustering
    # Foreign-key constraints:
    #    "yso_clustering_source_id_fkey" FOREIGN KEY (source_id)
    # REFERENCES gaia_dr2_source(source_id)

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                YSO_Clustering.twomass,
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag.alias("gaia_dr2_rp"),
                YSO_Clustering.j,
                YSO_Clustering.h,
                YSO_Clustering.k,
                Gaia_DR2.parallax,
            )
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .join(YSO_Clustering, on=(Gaia_DR2.source_id == YSO_Clustering.source_id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                Gaia_DR2.phot_rp_mean_mag < 15.5,
                YSO_Clustering.age < 7.5,
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

    def post_process(self, model):
        """
        cadence options for these targets:
        boss_bright_3x1 if RP<14.76 |
        boss_bright_4x1 if RP<15.075 |
        boss_bright_5x1 if RP<15.29 |
        boss_bright_6x1 if RP<15.5
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr2_rp from " + " sandbox.temp_mwm_yso_cluster_boss ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_rp = output[i][1]

            if current_rp < 14.76:
                current_instrument = "BOSS"
                current_cadence = "bright_3x1"
            elif current_rp < 15.075:
                current_instrument = "BOSS"
                current_cadence = "bright_4x1"
            elif current_rp < 15.29:
                current_instrument = "BOSS"
                current_cadence = "bright_5x1"
            elif current_rp < 15.5:
                current_instrument = "BOSS"
                current_cadence = "bright_6x1"
            else:
                # All cases should be covered above so we should not get here.
                current_instrument = None
                current_cadence = None
                raise TargetSelectionError(
                    "error in mwm_yso_cluster_boss "
                    + "post_process(): "
                    + "instrument = None, cadence= None"
                )

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_cluster_boss "
                    + " set instrument = '"
                    + current_instrument
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_cluster_boss "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )


class MWM_YSO_Cluster_BOSS_Single_Carton(MWM_YSO_Cluster_BOSS_Carton):
    name = "mwm_yso_cluster_boss_single"
    cadence = "bright_1x1"
    priority = 2706
    instrument = "BOSS"

    def post_process(self, model):
        pass


class MWM_YSO_PMS_APOGEE_Sagitta_EDR3_Carton(BaseCarton):
    """
    YSOs - Pre-main sequence, APOGEE
    Shorthand name: mwm_yso_pms_apogee_sagitta_edr3
    Comments: For v1.0, we use sagitta_edr3 instead of sagitta.
    Simplified Description of selection criteria:
    Selecting the clustered sources from the catalog of vetted
    pre-main sequence stars
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    https://wiki.sdss.org/pages/viewpage.action?spaceKey=OPS&title=Cartons+for+v1.0
    Additional source catalogs needed: catalogdb.sagitta_edr3
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    apogee_bright_3x1 (for 7 < H < 13)
    Implementation: (in sagitta_edr3) & h<13
    lead contact:Marina Kounkel
    """

    # peewee Model name ---> postgres table name
    # Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
    # Gaia_DR3(CatalogdbModel)--->'gaia_dr3_source'
    # Zari18pms(CatalogdbModel)--->'catalogdb.zari18pms'
    # Zari18ums(CatalogdbModel)--->'catalogdb.zari18ums'
    # Sagitta_EDR3(CatalogdbModel)--->'catalogdb.sagitta_edr3'
    # TwoMassPSC(CatalogdbModel)--->'catalogdb.twomass_psc'

    name = "mwm_yso_pms_apogee_sagitta_edr3"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2700
    can_offset = True

    def build_query(self, version_id, query_region=None):
        # join with Sagitta_EDR3 (we use Gaia_DR3 for query)
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag.alias("gaia_dr3_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR3.parallax,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .switch(Gaia_DR3)
            .join(Sagitta_EDR3, on=(Gaia_DR3.source_id == Sagitta_EDR3.source_id))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                TwoMassPSC.h_m < 13,
            )
        )

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


class MWM_YSO_PMS_APOGEE_Sagitta_EDR3_Single_Carton(MWM_YSO_PMS_APOGEE_Sagitta_EDR3_Carton):
    name = "mwm_yso_pms_apogee_sagitta_edr3_single"
    cadence = "bright_1x1"
    priority = 2701  # pms so different from other single cartons


class MWM_YSO_PMS_APOGEE_Zari18pms_Carton(BaseCarton):
    """
    YSOs - Pre-main sequence, APOGEE
    Shorthand name: mwm_yso_pms_apogee_zari18pms
    Comments: NA
    Simplified Description of selection criteria:
    Selecting the clustered sources from the catalog of vetted
    pre-main sequence stars
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    https://wiki.sdss.org/pages/viewpage.action?spaceKey=OPS&title=Cartons+for+v1.0
    Additional source catalogs needed: catalogdb.zari18pms
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    apogee_bright_3x1 (for 7 < H < 13)
    Implementation: (in zari18pms) & h<13
    lead contact:Marina Kounkel
    """

    # peewee Model name ---> postgres table name
    # Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
    # Gaia_DR3(CatalogdbModel)--->'gaia_dr3_source'
    # Zari18pms(CatalogdbModel)--->'catalogdb.zari18pms'
    # Zari18ums(CatalogdbModel)--->'catalogdb.zari18ums'
    # Sagitta_EDR3(CatalogdbModel)--->'catalogdb.sagitta_edr3'
    # TwoMassPSC(CatalogdbModel)--->'catalogdb.twomass_psc'

    name = "mwm_yso_pms_apogee_zari18pms"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_3x1"
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2700
    can_offset = True

    def build_query(self, version_id, query_region=None):
        # join with Zari18pms (we use Gaia_DR2 for query and not Gaia_DR3)
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                TwoMassPSC.pts_key,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag.alias("gaia_dr2_rp"),
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR2.parallax,
            )
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .switch(TIC_v8)
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .switch(Gaia_DR2)
            .join(Zari18pms, on=(Gaia_DR2.source_id == Zari18pms.source))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                TwoMassPSC.h_m < 13,
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


class MWM_YSO_PMS_APOGEE_Zari18pms_Single_Carton(MWM_YSO_PMS_APOGEE_Zari18pms_Carton):
    name = "mwm_yso_pms_apogee_zari18pms_single"
    cadence = "bright_1x1"
    priority = 2701  # pms so different from other single cartons


class MWM_YSO_PMS_BOSS_Sagitta_EDR3_Carton(BaseCarton):
    """
    YSOs - Pre-main sequence, BOSS
    Shorthand name: mwm_yso_pms_boss_sagitta_edr3
    Comments: For v1.0, we use sagitta_edr3 instead of sagitta.
    Simplified Description of selection criteria:
    Selecting the clustered sources from the catalog of vetted
    pre-main sequence stars
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    https://wiki.sdss.org/pages/viewpage.action?spaceKey=OPS&title=Cartons+for+v1.0
    Additional source catalogs needed: catalogdb.sagitta_edr3
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets:
    boss_bright_3x1 if RP<14.76 |
    boss_bright_4x1 if RP<15.075 |
    boss_bright_5x1 if RP<15.29 |
    boss_bright_6x1 if RP<15.5
    Implementation: (in sagitta_edr3) & rp<15.5
    lead contact:Marina Kounkel
    """

    # peewee Model name ---> postgres table name
    # Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
    # Gaia_DR3(CatalogdbModel)--->'gaia_dr3_source'
    # Zari18pms(CatalogdbModel)--->'catalogdb.zari18pms'
    # Zari18ums(CatalogdbModel)--->'catalogdb.zari18ums'
    # Sagitta_EDR3(CatalogdbModel)--->'catalogdb.sagitta_edr3'
    # TwoMassPSC(CatalogdbModel)--->'catalogdb.twomass_psc'

    name = "mwm_yso_pms_boss_sagitta_edr3"
    category = "science"
    instrument = None  # instrument is set in post_process()
    cadence = None  # cadence is set in post_process()
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2700
    can_offset = True

    def build_query(self, version_id, query_region=None):
        # join with Sagitta_EDR3 (we use Gaia_DR3 for query)
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag.alias("gaia_dr3_rp"),
                Gaia_DR3.parallax,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(Sagitta_EDR3, on=(Gaia_DR3.source_id == Sagitta_EDR3.source_id))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                Gaia_DR3.phot_rp_mean_mag < 15.5,
            )
        )

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

    def post_process(self, model):
        """
        cadence options for these targets:
        boss_bright_3x1 if RP<14.76 |
        boss_bright_4x1 if RP<15.075 |
        boss_bright_5x1 if RP<15.29 |
        boss_bright_6x1 if RP<15.5
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr3_rp from " + " sandbox.temp_mwm_yso_pms_boss_sagitta_edr3 ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_rp = output[i][1]

            if current_rp < 14.76:
                current_instrument = "BOSS"
                current_cadence = "bright_3x1"
            elif current_rp < 15.075:
                current_instrument = "BOSS"
                current_cadence = "bright_4x1"
            elif current_rp < 15.29:
                current_instrument = "BOSS"
                current_cadence = "bright_5x1"
            elif current_rp < 15.5:
                current_instrument = "BOSS"
                current_cadence = "bright_6x1"
            else:
                # All cases should be covered above so we should not get here.
                current_instrument = None
                current_cadence = None
                raise TargetSelectionError(
                    "error in mwm_yso_pms_boss_sagitta_edr3 "
                    + "post_process(): "
                    + "instrument = None, cadence= None"
                )

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_pms_boss_sagitta_edr3 "
                    + " set instrument = '"
                    + current_instrument
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_pms_boss_sagitta_edr3 "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )


class MWM_YSO_PMS_BOSS_Sagitta_EDR3_Single_Carton(MWM_YSO_PMS_BOSS_Sagitta_EDR3_Carton):
    name = "mwm_yso_pms_boss_sagitta_edr3_single"
    cadence = "bright_1x1"
    priority = 2701  # pms so different from other single cartons
    instrument = "BOSS"

    def post_process(self, model):
        pass


class MWM_YSO_PMS_BOSS_Zari18pms_Carton(BaseCarton):
    """
    YSOs - Pre-main sequence, BOSS
    Shorthand name: mwm_yso_pms_boss_zari18pms
    Comments: NA
    Simplified Description of selection criteria:
    Selecting the clustered sources from the catalog of vetted
    pre-main sequence stars
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    https://wiki.sdss.org/pages/viewpage.action?spaceKey=OPS&title=Cartons+for+v1.0
    Additional source catalogs needed: catalogdb.zari18pms
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets:
    boss_bright_3x1 if RP<14.76 |
    boss_bright_4x1 if RP<15.075 |
    boss_bright_5x1 if RP<15.29 |
    boss_bright_6x1 if RP<15.5
    Implementation: (in zari18pms) & rp<15.5
    lead contact:Marina Kounkel
    """

    # peewee Model name ---> postgres table name
    # Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
    # Gaia_DR3(CatalogdbModel)--->'gaia_dr3_source'
    # Zari18pms(CatalogdbModel)--->'catalogdb.zari18pms'
    # Zari18ums(CatalogdbModel)--->'catalogdb.zari18ums'
    # Sagitta_EDR3(CatalogdbModel)--->'catalogdb.sagitta_edr3'
    # TwoMassPSC(CatalogdbModel)--->'catalogdb.twomass_psc'

    name = "mwm_yso_pms_boss_zari18pms"
    category = "science"
    instrument = None  # instrument is set in post_process()
    cadence = None  # cadence is set in post_process()
    program = "mwm_yso"
    mapper = "MWM"
    priority = 2700
    can_offset = True

    def build_query(self, version_id, query_region=None):
        # join with Zari18pms (we use Gaia_DR2 for query and not Gaia_DR3)
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid,
                Gaia_DR2.source_id,
                Gaia_DR2.ra.alias("gaia_dr2_ra"),
                Gaia_DR2.dec.alias("gaia_dr2_dec"),
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag.alias("gaia_dr2_rp"),
                Gaia_DR2.parallax,
            )
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .switch(TIC_v8)
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .switch(Gaia_DR2)
            .join(Zari18pms, on=(Gaia_DR2.source_id == Zari18pms.source))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                Gaia_DR2.phot_rp_mean_mag < 15.5,
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

    def post_process(self, model):
        """
        cadence options for these targets:
        boss_bright_3x1 if RP<14.76 |
        boss_bright_4x1 if RP<15.075 |
        boss_bright_5x1 if RP<15.29 |
        boss_bright_6x1 if RP<15.5
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr2_rp from " + " sandbox.temp_mwm_yso_pms_boss_zari18pms ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_rp = output[i][1]

            if current_rp < 14.76:
                current_instrument = "BOSS"
                current_cadence = "bright_3x1"
            elif current_rp < 15.075:
                current_instrument = "BOSS"
                current_cadence = "bright_4x1"
            elif current_rp < 15.29:
                current_instrument = "BOSS"
                current_cadence = "bright_5x1"
            elif current_rp < 15.5:
                current_instrument = "BOSS"
                current_cadence = "bright_6x1"
            else:
                # All cases should be covered above so we should not get here.
                current_instrument = None
                current_cadence = None
                raise TargetSelectionError(
                    "error in mwm_yso_pms_boss_zari18pms "
                    + "post_process(): "
                    + "instrument = None, cadence= None"
                )

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_pms_boss_zari18pms "
                    + " set instrument = '"
                    + current_instrument
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_yso_pms_boss_zari18pms "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )


class MWM_YSO_PMS_BOSS_Zari18pms_Single_Carton(MWM_YSO_PMS_BOSS_Zari18pms_Carton):
    name = "mwm_yso_pms_boss_zari18pms_single"
    cadence = "bright_1x1"
    priority = 2701  # pms so different from other single cartons
    instrument = "BOSS"

    def post_process(self, model):
        pass
