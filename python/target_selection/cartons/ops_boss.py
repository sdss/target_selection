#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-07-19
# @Filename: ops_boss.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2, TIC_v8, TwoMassPSC,
                                             eBOSS_Target_v5)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py

# peewee Model name ---> postgres table name
# Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
# TwoMassPSC(CatalogdbModel) --->'twomass_psc'
# eBOSS_Target_v5(CatalogdbModel)--->'ebosstarget_v5'

class OPS_BOSS_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds
    Selection criteria:
    --- criteria for ops_BOSS_stds ---
    #calculate distance modulus (could convert to use BailerJones distance)
    distMod = 5.*np.log10(1000./gaia_DR2.parallax)-5.

    #calculate absolute g and k magnitudes
    abs_gmag = gaia_DR2.phot_g_mean_mag - distMod
    abs_kmag = twomass_psc.k_m - distMod

    meet_std_criteria =
    np.where( ( ( (gaia_DR2.phot_bp_mean_mag -
                   gaia_DR2.phot_rp_mean_mag) >= 0.65) &
                  ( (gaia_DR2.phot_bp_mean_mag -
                     gaia_DR2.phot_rp_mean_mag) <= 0.8) &
                  ( (abs_gmag) >= 3.5) & ( (abs_gmag) <= 5.5) )
             )

    Lead contact: Kevin Covey
    """

    name = 'ops_boss_stds'
    category = 'standard'
    cadence = None
    program = 'std'
    mapper = None

    def build_query(self, version_id, query_region=None):

        distMod = 5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0
        abs_gmag = Gaia_DR2.phot_g_mean_mag - distMod

        # We need a join with Catalog because
        # when using coordinates in the cartons we want to always
        # use Catalog.ra and Catalog.dec.
        # This is because those coordinates have
        # all been put in a common epoch 2015.5.
        # Below, we are using Catalog.ra and Catalog.dec
        # inside the "if query_region" block.

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.parallax > 0,
                        (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag) >= 0.65,
                        (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag) <= 0.8,
                        abs_gmag >= 3.5,
                        abs_gmag <= 5.5))
        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            # Here we do not need a join with Catalog since query already contains
            # a join with Catalog above.
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class OPS_BOSS_Red_Stds_No_Deredden_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_red_stds_no_deredden

    Selection Criteria:
    This carton OPS_BOSS_Red_Stds_No_Deredden_Carton is for the case bp_rp_excess < 0.
    bp_rp_excess is defined below.

    (Another carton OPS_BOSS_Red_Stds_Deredden_Carton is for the case bp_rp_excess >= 0)

    #calculate distance modulus (could convert to use BailerJones distance)
    distMod = 5.*np.log10(1000./gaia_DR2.parallax)-5.

    #calculate absolute g and k magnitudes
    abs_gmag = gaia_DR2.phot_g_mean_mag - distMod

    #calculate E(BP-RP) with (BP-RP)_0 = 0.725
    bp_rp_excess = bp - rp - 0.725
    #Here 0.725 = mode bp-rp color of eBOSS flux standards in SDSS footprint

    #only deredden if bp_rp_excess is positive
    worth_dereddening = np.where( bp_rp_excess >= 0)

    #calculate a_g and a_k using coefficients from Wang & Chen 2019
    #(use 0 otherwise; np.zero commands truncated for brevity)
    ag[worth_dereddening] = 1.890*bp_rp_excess[worth_dereddening]
    ak[worth_dereddening] = 0.186*bp_rp_excess[worth_dereddening]

    #This carton has worth_dereddening = False since the where clause has
    #Gaia_DR2.bp_rp-0.725 < 0
    #Hence,
    ag = 0
    ak = 0

    --- criteria for ops_BOSS_redstds ---

    meet_reddened_criteria =
    np.where( ( ((gaia_DR2.phot_g_mean_mag - ag) - (twomass_psc.k_m - ak)) >= 1.1) &
              ( ((gaia_DR2.phot_g_mean_mag - ag) - (twomass_psc.k_m - ak)) <= 1.6) &
                ( (abs_gmag - ag) >= 3) & ( (abs_gmag - ag) <= 5.5) )

    Lead contact: Kevin Covey
    """

    name = 'ops_boss_red_stds_no_deredden'
    category = 'standard'
    cadence = None
    program = 'std'
    mapper = None

    def build_query(self, version_id, query_region=None):

        distMod = 5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0
        abs_gmag = Gaia_DR2.phot_g_mean_mag - distMod
        bp_rp_excess = Gaia_DR2.bp_rp - 0.725
        ag = 0
        ak = 0

        # Even though ag=0 and ak=0, we keep ag and ak in the below expressions
        # so that it is easy to compare with the np.where() formula above.
        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.parallax > 0,
                        bp_rp_excess < 0,
                        ((Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak)) >= 1.1,
                        ((Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak)) <= 1.6,
                        (abs_gmag - ag) >= 3,
                        (abs_gmag - ag) <= 5.5))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class OPS_BOSS_Red_Stds_Deredden_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_red_stds_deredden
    
    Selection Criteria:
    This carton OPS_BOSS_Red_Stds_Deredden_Carton
    is for the case bp_rp_excess >= 0.
    bp_rp_excess is defined below.

    (Another carton OPS_BOSS_Red_Stds_No_Deredden_Carton
    is for the case bp_rp_excess < 0)

    #calculate distance modulus
    #(could convert to use BailerJones distance)
    distMod = 5.*np.log10(1000./gaia_DR2.parallax)-5.

    #calculate absolute g and k magnitudes
    abs_gmag = gaia_DR2.phot_g_mean_mag - distMod

    #calculate E(BP-RP) with (BP-RP)_0 = 0.725
    bp_rp_excess = bp - rp - 0.725
    #Here 0.725 = mode bp-rp color of eBOSS flux standards in SDSS footprint

    #only deredden if bp_rp_excess is positive
    worth_dereddening = np.where( bp_rp_excess >= 0)

    #calculate a_g and a_k using coefficients from Wang & Chen 2019
    #(use 0 otherwise; np.zero commands truncated for brevity)
    ag[worth_dereddening] = 1.890*bp_rp_excess[worth_dereddening]
    ak[worth_dereddening] = 0.186*bp_rp_excess[worth_dereddening]

    #This carton has worth_dereddening==True since the where clause has
    #Gaia_DR2.bp_rp-0.725 >= 0
    #Hence,
    ag = 1.890*bp_rp_excess = 1.890*(Gaia_DR2.bp_rp - 0.725)
    ak = 0.186*bp_rp_excess = 0.186*(Gaia_DR2.bp_rp - 0.725)

    --- criteria for ops_BOSS_redstds ---

    meet_reddened_criteria =
    np.where( ( ((gaia_DR2.phot_g_mean_mag - ag) - (twomass_psc.k_m - ak)) >= 1.1) &
              ( ((gaia_DR2.phot_g_mean_mag - ag) - (twomass_psc.k_m - ak)) <= 1.6) &
                ( (abs_gmag - ag) >= 3) & ( (abs_gmag - ag) <= 5.5) )

    Lead contact: Kevin Covey
    """

    name = 'ops_boss_red_stds_deredden'
    category = 'standard'
    cadence = None
    program = 'std'
    mapper = None

    def build_query(self, version_id, query_region=None):

        distMod = 5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0
        abs_gmag = Gaia_DR2.phot_g_mean_mag - distMod
        bp_rp_excess = Gaia_DR2.bp_rp - 0.725
        ag = 1.890 * (Gaia_DR2.bp_rp - 0.725)
        ak = 0.186 * (Gaia_DR2.bp_rp - 0.725)

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.parallax > 0,
                        bp_rp_excess >= 0,
                        ((Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak)) >= 1.1,
                        ((Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak)) <= 1.6,
                        (abs_gmag - ag) >= 3,
                        (abs_gmag - ag) <= 5.5))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class OPS_eBOSS_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_eboss_stds
    Selection Criteria:
    The code of this carton is based on the below SQL.
    This returns 298885 rows.
    SELECT DISTINCT ON (e.objid_targeting) e.objid_targeting, c2t.catalogid
           FROM ebosstarget_v5 e JOIN tic_v8 t ON t.sdss = e.objid_targeting
           JOIN catalog_to_tic_v8 c2t ON c2t.target_id = t.id
           WHERE (e.eboss_target1 & pow(2, 50)::bigint) > 0 OR
                (e.eboss_target1 & pow(2, 51)::bigint) > 0 OR
                (e.eboss_target1 & pow(2, 52)::bigint) > 0;

    As shown below, more than one c2t.catalogid may correspond
    to the same e.objid_targeting.
    SELECT count(c2t.catalogid ), e.objid_targeting
    FROM ebosstarget_v5 e JOIN tic_v8 t ON t.sdss = e.objid_targeting
    JOIN catalog_to_tic_v8 c2t ON c2t.target_id = t.id
    WHERE (e.eboss_target1 & pow(2, 50)::bigint) > 0 OR
          (e.eboss_target1 & pow(2, 51)::bigint) > 0 OR
          (e.eboss_target1 & pow(2, 52)::bigint) > 0  GROUP BY e.objid_targeting;

    count |   objid_targeting
    -------+---------------------
        2 | 1237662530065399951
        2 | 1237671956455489893
    etc.

    Hence,
    SELECT COUNT(DISTINCT e.objid_targeting) etc.
    returns 298885

    and
    SELECT count(DISTINCT c2t.catalogid ) etc.
    returns 642787
    lead contact: Kevin Covey
    """

    name = 'ops_eboss_stds'
    category = 'standard'
    cadence = None
    program = 'std'
    mapper = None

    def build_query(self, version_id, query_region=None):

        selection_condition = (
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 50).
             cast("bigint")) > 0) |
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 51).
             cast("bigint")) > 0) |
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 52).
             cast("bigint")) > 0))

        # We have distinct(eBOSS_Target_v5.objid_targeting) at the end
        # since the table ebosstarget_v5 has duplicate values.
        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(eBOSS_Target_v5,
                       on=(TIC_v8.sdss == eBOSS_Target_v5.objid_targeting))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        selection_condition)
                 .distinct(eBOSS_Target_v5.objid_targeting))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query
