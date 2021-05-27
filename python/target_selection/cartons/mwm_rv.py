#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_rv.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import math

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (AllWise, Catalog, CatalogToTIC_v8,
                                             Gaia_DR2,
                                             SDSS_APOGEE_AllStarMerge_r13,
                                             TIC_v8, TwoMassPSC)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py

# 2.2.1. Long Baseline (Legacy Targets)

# Shorthand name: mwm_rv_long (Not an actual target class.
# Defined as a parent sample to select from for 2.2.1.x sections below)

# Simplified Description of selection criteria:
# APOGEE-1&-2 "main survey" targets with at least 3 previous epochs
# brighter than H of 12.8 that have a Gaia parallax.

# Wiki page: Binaries and Substellar Companions

# Additional source catalogs needed: Just need sdss_apogeeallstarmerge_r13

# Additional cross-matching needed: (None)

# Return columns: apogee_id, nvisits, ra, dec, pmra, pmdec, h, baseline, fields

# cadence options for these targets (list all options,
# even though no single target will receive more than one):
# (See entries for this in 2.2.1.x below)

# Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
# pmra, pmdec, h, baseline, fields FROM sdss_apogeeallstarmerge_r13
# WHERE h<12.8 AND nvisits>=3 AND dist_src == 'gaia' AND
# [targflags contains 'APOGEE_SHORT' OR 'APOGEE_INTERMEDIATE' OR
# 'APOGEE_LONG' OR '*APOGEE2_*BIN_’] AS mwm_rv_long

# Implementation:

# Non-SQL implementation:

# lead contact: Nicholas Troup
# Target selection final for v0?: No

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

mwm_rv_long_condition = (SDSS_APOGEE_AllStarMerge_r13.h < 12.8,
                         SDSS_APOGEE_AllStarMerge_r13.nvisits >= 3,
                         peewee.fn.trim(SDSS_APOGEE_AllStarMerge_r13.dist_src) ==
                         'gaia',
                         (SDSS_APOGEE_AllStarMerge_r13.targflags %
                          '%APOGEE_SHORT%') |
                         (SDSS_APOGEE_AllStarMerge_r13.targflags %
                          '%APOGEE_INTERMEDIATE%') |
                         (SDSS_APOGEE_AllStarMerge_r13.targflags %
                          '%APOGEE_LONG%') |
                         (SDSS_APOGEE_AllStarMerge_r13.targflags %
                          '%APOGEE2_%BIN_%'))


class MWM_RV_Long_FPS_Carton(BaseCarton):
    """3.2.1.3. Long Baseline Targets for FPS

    Shorthand name: mwm_rv_long_fps

    Simplified Description of selection criteria:
     Select from long-baseline targets (above) with H brighter than 11.5

    Wiki page: Binaries and Substellar Companions (Page Under Construction)

    Additional source catalogs needed: Select from mwm_rv_long (2.2.1)

    Additional cross-matching needed: (None)

    Return columns: apogee_id, nvisits, ra, dec,
     pmra, pmdec, h, baseline, fields (Same as 2.2.1 above)

    cadence options for these targets
     (list all options, even though no single target will receive more than one):
      If H>10.8 then use mwm_rv_<nn>x2, otherwise use mwm_rv_<nn>x1,
       where <nn> = 3*ceiling((18-nvisits)/3)
      if <nn> is less than 6 then
         set <nn> = 6

    Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
     pmra, pmdec, h, baseline, fields FROM mwm_rv_long WHERE h<11.5

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
    """
    name = 'mwm_rv_long_fps'
    category = 'science'
    instrument = 'APOGEE'
    cadence = None  # cadence is set in post_process()
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2500

    # peewee Model name ---> postgres table name
    # SDSS_APOGEE_AllStarMerge_r13(CatalogdbModel)--->'sdss_apogeeallstarmerge_r13'

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid,
                         SDSS_APOGEE_AllStarMerge_r13.apogee_id,
                         SDSS_APOGEE_AllStarMerge_r13.nvisits,
                         SDSS_APOGEE_AllStarMerge_r13.ra.alias('allstarmerge_ra'),
                         SDSS_APOGEE_AllStarMerge_r13.dec.alias('allstarmerge_dec'),
                         SDSS_APOGEE_AllStarMerge_r13.pmra.alias('allstarmerge_pmra'),
                         SDSS_APOGEE_AllStarMerge_r13.pmdec.alias('allstarmerge_pmdec'),
                         SDSS_APOGEE_AllStarMerge_r13.h,
                         SDSS_APOGEE_AllStarMerge_r13.baseline,
                         SDSS_APOGEE_AllStarMerge_r13.fields)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join(SDSS_APOGEE_AllStarMerge_r13,
                       on=(TwoMassPSC.designation ==
                           peewee.fn.ltrim(SDSS_APOGEE_AllStarMerge_r13.apogee_id,
                                           '2M')))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        *mwm_rv_long_condition,
                        SDSS_APOGEE_AllStarMerge_r13.h < 11.5))
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

    def post_process(self, model):
        """
        If H>10.8 then use bright_<nn>x2, otherwise use bright_<nn>x1,
        where <nn> = 3*ceiling((18-nvisits)/3)
        if <nn> is less than 6 then
            set <nn> = 6
        """

        cursor = self.database.execute_sql(
            "select catalogid, nvisits, h from " +
            " sandbox.temp_mwm_rv_long_fps ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_nvisits = output[i][1]
            current_h = output[i][2]

            nn = 3 * math.ceil((18 - current_nvisits) / 3)
            if(nn < 6):
                nn = 6

            if(current_h > 10.8):
                current_cadence = 'bright_' + str(nn) + 'x2'
            else:
                current_cadence = 'bright_' + str(nn) + 'x1'

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_rv_long_fps " +
                    " set cadence = '" + current_cadence + "'"
                    " where catalogid = " + str(current_catalogid) + ";")


# 2.2.2. Short Baseline (Fresh Targets)
#
# Shorthand name: mwm_rv_short (Not an actual target class.
# Defined as a parent sample to select from)
#
# Simplified Description of selection criteria:
# 2MASS Targets between H = 7 and H=12.8 with the same
# selection criteria as APOGEE-2 main survey
# with the additional requirement of having a Gaia parallax.
#
# Wiki page: Binaries and Substellar Companions (Page Under Construction)
#
# Additional source catalogs needed: 2MASS, WISE, Gaia
#
# Additional cross-matching needed:
# The 3 catalogs above should be cross-matched if not already
#
# Return columns: ra, dec, pmra, pmdec, (from Gaia); h_m (from 2MASS)
#
# cadence options for these targets
# (list all options, even though no single target will
# receive more than one):  (See entries for this in 2.2.1.x below)
#
# Pseudo SQL (optional):  SELECT catalog.catalogid,
# catalog.ra, catalog.ra, catalog.dec, catalog.pmra, catalog.pmdec,
# catalog.epoch, twomass_psc.h_m, FROM catalog
# JOIN gaia_dr2_source, twomass_psc, allWise
# WHERE h_m < 12.8 AND ((j_m-k_m) - (1.5*0.918*(h_m-w2mpro-0.05))) >= 0.5
# AND (j_msigcom <= 0.1 AND h_msigcom<=0.1 AND k_msigcom <= 0.1) AND w2_sigmpro <=0.1
# AND (ph_qual= 'AAA' OR ph_qual= 'AAB' OR ph_qual= 'ABA'
# OR ph_qual= 'BAA' OR ph_qual= 'ABB' OR ph_qual= 'BAB'
# OR ph_qual= 'BBA' OR ph_qual = 'BBB')
# AND prox >= 6 AND cc_flag='000' AND gal_contam='000'
# AND (rd_flg='111' OR rd_flg='112'  OR rd_flg='121'  OR rd_flg='211'
#   OR rd_flg='122'  OR rd_flg='212'  OR rd_flg='221'  OR rd_flg='222')
# AND ext_key=Null AS mwm_rv_short
#
# Implementation:
#
# Non-SQL implementation:
#
# lead contact: Nicholas Troup
#
# Target selection final for v0?: No

mwm_rv_short_condition = (TwoMassPSC.h_m < 12.8,
                          ((TwoMassPSC.j_m - TwoMassPSC.k_m) -
                           (1.5 * 0.918 *
                            (TwoMassPSC.h_m - AllWise.w2mpro - 0.05))) >= 0.5,
                          TwoMassPSC.j_msigcom <= 0.1,
                          TwoMassPSC.h_msigcom <= 0.1,
                          TwoMassPSC.k_msigcom <= 0.1,
                          AllWise.w2sigmpro <= 0.1,
                          ((TwoMassPSC.ph_qual == 'AAA') |
                           (TwoMassPSC.ph_qual == 'AAB') |
                           (TwoMassPSC.ph_qual == 'ABA') |
                           (TwoMassPSC.ph_qual == 'BAA') |
                           (TwoMassPSC.ph_qual == 'ABB') |
                           (TwoMassPSC.ph_qual == 'BAB') |
                           (TwoMassPSC.ph_qual == 'BBA') |
                           (TwoMassPSC.ph_qual == 'BBB')),
                          TwoMassPSC.prox >= 6,
                          TwoMassPSC.cc_flg == '000',
                          TwoMassPSC.gal_contam == 0,
                          ((TwoMassPSC.rd_flg == '111') |
                           (TwoMassPSC.rd_flg == '112') |
                           (TwoMassPSC.rd_flg == '121') |
                           (TwoMassPSC.rd_flg == '211') |
                           (TwoMassPSC.rd_flg == '122') |
                           (TwoMassPSC.rd_flg == '212') |
                           (TwoMassPSC.rd_flg == '221') |
                           (TwoMassPSC.rd_flg == '222')),
                          TwoMassPSC.ext_key >> None,
                          Gaia_DR2.parallax.is_null(False))


class MWM_RV_Short_FPS_Carton(BaseCarton):
    """3.2.2.3. Short Baseline Targets for FPS

    Shorthand name: mwm_rv_short_fps

    Simplified Description of selection criteria:
     Select from short-baseline targets (above) with H brighter than 10.8

    Wiki page: Binaries and Substellar Companions (Page Under Construction)

    Additional source catalogs needed: Select from mwm_rv_short (2.2.2)

    Additional cross-matching needed: (None)

    cadence options for these targets (list all options,
     even though no single target will receive more than one): mwm_rv_18x1

    Pseudo SQL (optional): SELECT catalogid, ra, dec,
     pmra, pmdec, epoch, h_m FROM mwm_rv_short WHERE h_m<10.8

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
    """
    name = 'mwm_rv_short_fps'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'bright_18x1'
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2510

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.ra.alias('gaia_dr2_ra'),
                         Gaia_DR2.dec.alias('gaia_dr2_dec'),
                         Gaia_DR2.pmra.alias('gaia_dr2_pmra'),
                         Gaia_DR2.pmdec.alias('gaia_dr2_pmdec'),
                         TwoMassPSC.h_m.alias('twomass_h_m'))
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(AllWise, on=(TIC_v8.allwise == AllWise.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        *mwm_rv_short_condition,
                        TwoMassPSC.h_m < 10.8))
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
