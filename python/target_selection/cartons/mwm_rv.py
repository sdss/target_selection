#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_rv.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import math

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (AllWise, Catalog,
                                             CatalogToGaia_DR3,
                                             CatalogToTIC_v8,
                                             CatalogToTwoMassPSC, Gaia_DR3,
                                             SDSS_DR17_APOGEE_Allstarmerge,
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

# Additional source catalogs needed: Just need sdss_dr17_apogee_allstarmerge

# Additional cross-matching needed: (None)

# Return columns: apogee_id, nvisits, ra, dec, pmra, pmdec, h, baseline, fields

# cadence options for these targets (list all options,
# even though no single target will receive more than one):
# (See entries for this in 2.2.1.x below)

# Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
# pmra, pmdec, h, baseline, fields FROM sdss_dr17_apogee_allstarmerge
# WHERE h<12.8 AND nvisits>=3 AND dist_src == 'GaiaDR3' AND
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
# We for sdss_dr17_apogee_allstarmerge we do not
# have to strip off the “2M” in the APOGEE_ID.
# (old: for sdss_apogeeallstarmerge_r13 we had to strip off the
#  “2M” in the APOGEE_ID.)
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# sdss5db=# select apogee_id from catalogdb.sdss_dr17_apogee_allstarmerge limit 2;
#    apogee_id
# ------------------
#  19140272-1554055
#  19155129-1617591
# (2 rows)

# ######### start comment for old sdss_apogeeallstarmerge_r13 #######################
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
# ######### end comment for old sdss_apogeeallstarmerge_r13 ##########

# Below we use GaiaEDR3 and not GaiaDR3 since the dist_src column
# does not have GaiaDR3

mwm_rv_long_condition = (SDSS_DR17_APOGEE_Allstarmerge.h < 11.5,  # old 12.8,
                         SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 6,  # old 3,
                         peewee.fn.trim(SDSS_DR17_APOGEE_Allstarmerge.dist_src) ==
                         'GaiaEDR3',
                         (SDSS_DR17_APOGEE_Allstarmerge.targflags %
                          '%APOGEE_SHORT%') |
                         (SDSS_DR17_APOGEE_Allstarmerge.targflags %
                          '%APOGEE_INTERMEDIATE%') |
                         (SDSS_DR17_APOGEE_Allstarmerge.targflags %
                          '%APOGEE_LONG%') |
                         (SDSS_DR17_APOGEE_Allstarmerge.targflags %
                          '%APOGEE2_%BIN_%'))


# old class name MWM_RV_Long_FPS_Carton(BaseCarton):
class MWM_bin_rv_long_Carton(BaseCarton):
    """3.2.1.3. Long Baseline Targets for FPS

    Shorthand name: mwm_bin_rv_long (old name mwm_rv_long_fps)

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

    name = 'mwm_bin_rv_long'  # old name = 'mwm_rv_long_fps'
    category = 'science'
    instrument = 'APOGEE'
    cadence = None  # cadence is set in post_process()
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = None  # priority is set in post_process()

    # peewee Model name ---> postgres table name
    # SDSS_DR17_APOGEE_Allstarmerge(CatalogdbModel)--->'sdss_dr17_apogee_allstarmerge'
    # There is no gaia_dr3 in the below query so we do not have to do
    # major modification of the old query.
    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid,
                         SDSS_DR17_APOGEE_Allstarmerge.apogee_id,
                         SDSS_DR17_APOGEE_Allstarmerge.nvisits,
                         SDSS_DR17_APOGEE_Allstarmerge.ra.alias('allstarmerge_ra'),
                         SDSS_DR17_APOGEE_Allstarmerge.dec.alias('allstarmerge_dec'),
                         SDSS_DR17_APOGEE_Allstarmerge.pmra.alias('allstarmerge_pmra'),
                         SDSS_DR17_APOGEE_Allstarmerge.pmdec.alias('allstarmerge_pmdec'),
                         SDSS_DR17_APOGEE_Allstarmerge.h,
                         SDSS_DR17_APOGEE_Allstarmerge.baseline,
                         SDSS_DR17_APOGEE_Allstarmerge.fields,
                         SDSS_DR17_APOGEE_Allstarmerge.teff_avg,  # old teff
                         SDSS_DR17_APOGEE_Allstarmerge.logg_avg,  # old logg
                         SDSS_DR17_APOGEE_Allstarmerge.dist,
                         SDSS_DR17_APOGEE_Allstarmerge.dist_src
                         )
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join(SDSS_DR17_APOGEE_Allstarmerge,
                       on=(TwoMassPSC.designation ==
                           peewee.fn.ltrim(SDSS_DR17_APOGEE_Allstarmerge.apogee_id, '2M')))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        *mwm_rv_long_condition,
                        SDSS_DR17_APOGEE_Allstarmerge.h < 11.5))
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
        For cadence:
        If H>10.8 then use bright_<nn>x2, otherwise use bright_<nn>x1,
        where <nn> = 3*ceiling((18-nvisits)/3)
        if <nn> is less than 6 then
            set <nn> = 6

        Teff means teff_avg of SDSS_DR17_APOGEE_Allstarmerge.
        logg measn logg_avg of SDSS_DR17_APOGEE_Allstarmerge.
        For priority:
        IF Teff < 4500 AND logg > 4.0 THEN priority = 2510
        ELSE IF 3.5 <= logg <= 4.0 THEN priority = 2520
        ELSE IF logg < 3.5 THEN priority = 2530
        ELSE priority = 2540
        """

        default_priority = 2540

        # teff_avg and logg_avg are from SDSS_DR17_APOGEE_Allstarmerge
        # old name was teff, logg
        cursor = self.database.execute_sql(
            "select catalogid, nvisits, h, teff_avg, logg_avg from " +
            " sandbox.temp_mwm_bin_rv_long ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_nvisits = output[i][1]
            current_h = output[i][2]
            current_teff_avg = output[i][3]
            current_logg_avg = output[i][4]

            nn = 3 * math.ceil((18 - current_nvisits) / 3)
            if (nn < 6):
                nn = 6

            if (current_h > 10.8):
                current_cadence = 'bright_' + str(nn) + 'x2'
            else:
                current_cadence = 'bright_' + str(nn) + 'x1'

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_bin_rv_long " +
                    " set cadence = '" + current_cadence + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if (current_logg_avg is not None):
                if ((current_teff_avg is not None) and (current_teff_avg < 4500) and
                   (current_logg_avg > 4.0)):
                    current_priority = 2510
                elif ((3.5 <= current_logg_avg) and (current_logg_avg <= 4.0)):
                    current_priority = 2520
                elif (current_logg_avg < 3.5):
                    current_priority = 2530
                else:
                    current_priority = default_priority
            else:
                current_priority = default_priority

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_bin_rv_long " +
                    " set priority = " + str(current_priority) +
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
# JOIN gaia_dr3_source, twomass_psc, allWise
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
                          Gaia_DR3.parallax.is_null(False))


class MWM_bin_rv_short_Carton(BaseCarton):
    # old class name MWM_RV_Short_FPS_Carton(BaseCarton):
    """3.2.2.3. Short Baseline Targets for FPS

    Shorthand name: mwm_bin_rv_short (old name = 'mwm_rv_short_fps')

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

    name = 'mwm_bin_rv_short'  # old name = 'mwm_rv_short_fps'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'bright_18x1'
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = None  # priority is set in post_process()

    # This carton i.e. mwm_bin_rv_short does not use SDSS_DR17_APOGEE_Allstarmerge.
    # There is gaia_dr3 in the below query so we have done
    # major modification of the old query.
    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.pmra.alias('gaia_dr3_pmra'),
                         Gaia_DR3.pmdec.alias('gaia_dr3_pmdec'),
                         TwoMassPSC.h_m.alias('twomass_h_m'),
                         TIC_v8.teff,
                         TIC_v8.logg)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .switch(Catalog)
                 .join(CatalogToGaia_DR3,
                       on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .switch(Catalog)
                 .join(CatalogToTwoMassPSC,
                       on=(Catalog.catalogid == CatalogToTwoMassPSC.catalogid))
                 .join(TwoMassPSC,
                       on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
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

    def post_process(self, model):
        """
        For priority:
        IF Teff < 4500 AND logg > 4.0 THEN priority = 2515
        ELSE IF 3.5 <= logg <= 4.0 THEN priority = 2525
        ELSE IF logg < 3.5 THEN priority = 2535
        ELSE priority = 2545
        """

        default_priority = 2545

        # teff and logg are from TIC_v8
        cursor = self.database.execute_sql(
            "select catalogid, teff, logg from " +
            " sandbox.temp_mwm_bin_rv_short ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_teff = output[i][1]
            current_logg = output[i][2]

            if (current_logg is not None):
                if ((current_teff is not None) and (current_teff < 4500) and
                   (current_logg > 4.0)):
                    current_priority = 2515
                elif ((3.5 <= current_logg) and (current_logg <= 4.0)):
                    current_priority = 2525
                elif (current_logg < 3.5):
                    current_priority = 2535
                else:
                    current_priority = default_priority
            else:
                current_priority = default_priority

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_bin_rv_short " +
                    " set priority = " + str(current_priority) +
                    " where catalogid = " + str(current_catalogid) + ";")
