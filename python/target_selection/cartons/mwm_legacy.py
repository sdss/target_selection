#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-10-16
# @Filename: mwm_legacy.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToTIC_v8, Gaia_DR2,
                                             SDSS_APOGEE_AllStarMerge_r13,
                                             TIC_v8
                                             )

from target_selection.cartons import BaseCarton

# TODO
# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee model : postgres table name
# Gaia_DR2 : catalogdb.gaia_dr2_source  
# SDSS_APOGEE_AllStarMerge_r13 : catalogdb.sdss_apogeeallstarmerge_r13 
# Catalog : catalogdb.catalog 
# CatalogToTIC_v8 : catalogdb.catalog_to_tic_v8 
# TIC_v8 : catalogdb.tic_v8
#         : catalogdb.catalog_to_sdss_dr16_apogeestar   
#
class MWM_Legacy_ir2opt_Carton(BaseCarton):
    """MWM APOGEE targets.

Shorthand name: mwm_legacy_ir2opt

What is it?: APOGEE targets to be observed with BOSS fibres.
To be used as a filler for MWM-led plates.

Simplified Description of selection criteria: Select all Gaia targets that
have an APOGEE counterpart (i.e., have an entry in sdss_apogeeallstarmerge_r13)
with 14 < G < 18, BP > 13, and RP > 13.

Wiki page: NA

Additional source catalogs needed: Gaia, sdss_apogeeallstarmerge_r13

Additional cross-matching needed: NA

Return columns: catalog_id, source_id, apogee_id, phot_g_mean_mag,
phot_rp_mean_mag, phot_bp_mean_mag

cadence options for these targets 
(list all options, even though no single target will receive more than one):  
NA
    """

    name = 'mwm_legacy_ir2opt'
    category = 'science'
    cadence = None
    program = 'mwm_legacy'
    mapper = 'MWM'
    priority = None

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(SDSS_APOGEE_AllStarMerge_r13, TODO on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(AllWise, on=(TIC_v8.allwise == AllWise.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 13,
                        (AllWise.w1mpro - AllWise.w2mpro) > 0.25,
                        (AllWise.w2mpro - AllWise.w3mpro) > 0.50,
                        (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                        Gaia_DR2.parallax > 0.3))

        # Gaia_DR2 peewee model class corresponds to
        # table catalogdb.gaia_dr2_source.

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


