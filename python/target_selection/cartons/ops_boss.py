#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-07-19
# @Filename: ops_boss.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2,
                                             TIC_v8, TwoMassPSC)

# from target_selection import log
from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


# --- variables invoked in the criteria below ---
#
#  #step not shown:
# using gaia+2mass crossmatch primary key to link gaia+2mass counterparts.
#
#  #calculate distance modulus (could convert to use BailerJones distance)
#  distMod = 5.*np.log10(1000./gaia_DR2.parallax)-5.
#
#  #calculate absolute g and k magnitudes
#  abs_gmag = gaia_DR2.phot_g_mean_mag - distMod
#  abs_kmag = twomass_psc.k_m - distMod
#
#  #calculate E(BP-RP) with (BP-RP)_0 = 0.725
#  bp_rp_excess = bp - rp - 0.725
# 0.725 = mode bp-rp color of eBOSS flux standards in SDSS footprint
#
#  #only deredden if bp_rp_excess is positive
#  worth_dereddening = np.where( bp_rp_excess >= 0)
#
#  #calculate a_g and a_k using coefficients from Wang & Chen 2019
#  (use 0 otherwise; np.zero commands truncated for brevity)
#  ag[worth_dereddening] = 1.890*bp_rp_excess[worth_dereddening]
#  ak[worth_dereddening] = 0.186*bp_rp_excess[worth_dereddening]
#  abp[worth_dereddening] = 2.429*bp_rp_excess[worth_dereddening]
#  arp[worth_dereddening] = 1.429*bp_rp_excess[worth_dereddening]
#
#
#


class OPS_BOSS_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds
    lead contact: Kevin Covey
    """

    # peewee Model name ---> postgres table name
    # Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
    # TwoMassPSC(CatalogdbModel) --->'twomass_psc'

    # --- criteria for ops_BOSS_stds ---
#
# meet_std_criteria =
# np.where( ( ( (gaia_DR2.phot_bp_mean_mag -
#                gaia_DR2.phot_rp_mean_mag) >= 0.65) &
#               ( (gaia_DR2.phot_bp_mean_mag -
#                  gaia_DR2.phot_rp_mean_mag) <= 0.8) &
#               ( (abs_gmag) >= 3.5) & ( (abs_gmag) <= 5.5) )
#          )
#
#
#
#
    name = 'ops_boss_stds'
    category = 'OPS'
    cadence = None
    program = 'OPS'
    mapper = 'OPS'

    def build_query(self, version_id, query_region=None):

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
                        (Gaia_DR2.phot_g_mean_mag -
                         (5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0)) >= 3.5,
                        (Gaia_DR2.phot_g_mean_mag -
                         (5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0)) <= 5.5))
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


class OPS_BOSS_Red_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_red_stds
    lead contact: Kevin Covey
    """


# --- criteria for ops_BOSS_redstds ---
#
# meet_reddened_criteria =
# np.where( ( ((gaia_DR2.phot_g_mean_mag - ag) - (twomass_psc.k_m - ak)) >= 1.1) &
#           ( ((gaia_DR2.phot_g_mean_mag - ag) - (twomass_psc.k_m - ak)) <= 1.6) &
#             ( (abs_gmag - ag) >= 3) & ( (abs_gmag - ag) <= 5.5) )
#

    name = 'mwm_boss_red_stds'
    category = 'OPS'
    cadence = None
    program = 'OPS'
    mapper = 'OPS'

    def build_query(self, version_id, query_region=None):

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
                        ((Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak)) >= 1.1,
                        ((Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak)) <= 1.6,
                        ((Gaia_DR2.phot_g_mean_mag -
                          (5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0)) - ag) >= 3,
                        ((Gaia_DR2.phot_g_mean_mag -
                          (5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0)) - ag) <= 5.5))

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
