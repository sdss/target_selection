#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-10-16
# @Filename: mwm_legacy.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2,
                                             SDSS_APOGEE_AllStarMerge_r13,
                                             TIC_v8, TwoMassPSC)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee model : postgres table name
# Catalog : catalogdb.catalog
# CatalogToTIC_v8 : catalogdb.catalog_to_tic_v8
# Gaia_DR2 : catalogdb.gaia_dr2_source
# SDSS_APOGEE_AllStarMerge_r13 : catalogdb.sdss_apogeeallstarmerge_r13
# TIC_v8 : catalogdb.tic_v8
# TwoMassPSC : catalogdb.twomass_psc
#

# We associate each APOGEE target with a 2MASS source since the
# vast majority of APOGEE targets were selected from the 2MASS PSC.
# APOGEE_ID is essentially the same as
# the “designation” column of the 2MASS PSC;
# We just have to strip off the “2M” in the APOGEE_ID.
#
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# sdss5db=# select apogee_id from  catalogdb.sdss_apogeeallstarmerge_r13 limit 2;
#      apogee_id
# --------------------
#  2M14044120-1550575
#  2M14033676-1554164
# (2 rows)
#
# sdss5db=# select ltrim(apogee_id,'2M') from
#  catalogdb.sdss_apogeeallstarmerge_r13 limit 2;
#       ltrim
# ------------------
#  14044120-1550575
#  14033676-1554164
# (2 rows)
#
#
# sdss5db=# select count(1) from
#  catalogdb.sdss_apogeeallstarmerge_r13 where dist_src='gaia';
#  count
# -------
#      0
# (1 row)
#
# Hence use trim(dist_src) like below:
#
# sdss5db=# select count(1) from
#  catalogdb.sdss_apogeeallstarmerge_r13 where trim(dist_src)='gaia';
#  count
# --------
#  487508
# (1 row)
#


class MWM_Legacy_ir2opt_Carton(BaseCarton):
    """MWM APOGEE targets.

Shorthand name: mwm_legacy_ir2opt

What is it?: APOGEE targets to be observed with BOSS fibres.
To be used as a filler for MWM-led plates.

Simplified Description of selection criteria: Select all Gaia targets that
have an APOGEE counterpart (i.e., have an entry in sdss_apogeeallstarmerge_r13)
with 13 < G < 18, BP > 13, and RP > 13.

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
    priority = 6100

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         TwoMassPSC.pts_key)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join(SDSS_APOGEE_AllStarMerge_r13,
                       on=(TwoMassPSC.designation ==
                           peewee.fn.ltrim(SDSS_APOGEE_AllStarMerge_r13.apogee_id,
                                           '2M')))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.phot_g_mean_mag.between(13, 18),
                        Gaia_DR2.phot_bp_mean_mag > 13,
                        Gaia_DR2.phot_rp_mean_mag > 13))

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
