#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_rv.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from astropy.coordinates import Angle

from sdssdb.peewee.sdss5db.catalogdb import (AllWise, Catalog, CatalogToTIC_v8,
                                             Gaia_DR2,
                                             SDSS_APOGEE_AllStarMerge_r13,
                                             TIC_v8, TwoMassPSC)

from target_selection import log
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


class MWM_RV_Long_RM_Carton(BaseCarton):
    """2.2.1.1. Long Baseline Targets for RM Plates

    Shorthand name: mwm_rv_long_rm

    Simplified Description of selection criteria:
    Select from long-baseline targets (above) within
    the 3 BHM-RM fields observable from APO
    (sexagesimal coordinates given in Pseudo SQL below)

    Wiki page: Binaries and Substellar Companions

    Additional source catalogs needed: Select from mwm_rv_long (2.2.1)

    Additional cross-matching needed: (None)

    Return columns: apogee_id, nvisits, ra, dec, pmra, pmdec, h,
    baseline, fields (Same as 2.2.1 above)

    cadence options for these targets (list all options,
    even though no single target will receive more than one):
    mwm_rv_<nn>x8 (where <nn> = number of visits to RM fields)

    Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
    pmra, pmdec, h, baseline, fields FROM mwm_rv_long
    WHERE [(ra,dec) < 3 degrees on sky from
    (14:14:49 +53:05:00) OR (10:00:00 +02:12:00) OR (02:23:30 -04:15:00)]

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
    """
    name = 'mwm_rv_long_rm'
    category = 'science'
    cadence = None
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2571

    # peewee Model name ---> postgres table name
    # SDSS_APOGEE_AllStarMerge_r13(CatalogdbModel)--->'sdss_apogeeallstarmerge_r13'

    def build_query(self, version_id, query_region=None):
        # WHERE [(ra,dec) < 3 degrees on sky from
        # (14:14:49 +53:05:00) OR (10:00:00 +02:12:00) OR (02:23:30 -04:15:00)]

        a_ra = Angle('14h14m49s').degree
        a_dec = Angle('53d05m00s').degree

        b_ra = Angle('10h00m00s').degree
        b_dec = Angle('02d12m00s').degree

        c_ra = Angle('02h23m30s').degree
        c_dec = Angle('-04d15m00s').degree

        # In the query for cartons in mwmw_yso.py,
        # we do not need a join with Catalog
        # because the query does not use ra and dec in the WHERE clause.
        # Here, for RV cartons, we are using ra and dec in the WHERE clause in
        # q3c_radial_query() . Hence we need a join with Catalog because
        # when using coordinates in the cartons we want to always
        # use Catalog.ra and Catalog.dec.
        # This is because those coordinates have
        # all been put in a common epoch 2015.5.

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
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
                        peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, a_ra, a_dec, 3) |
                        peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, b_ra, b_dec, 3) |
                        peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, c_ra, c_dec, 3)))
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


class MWM_RV_Long_Bplates_Carton(BaseCarton):
    """ 3.2.1.2. Long Baseline Targets for Bright Time Plates

    Shorthand name: mwm_rv_long_bplates

    Simplified Description of selection criteria:
    Select from long-baseline targets (above) observable
     from APO with H brighter than 12.2 with at least 6 APOGEE
      visits within selected APOGEE fields observable during 6-month program
       (field centers in decimal degrees given below).

    Wiki page: Binaries and Substellar Companions (Page Under Construction)

    Additional source catalogs needed: Select from mwm_rv_long (2.2.1)

    Additional cross-matching needed: (None)

    Return columns: apogee_id, nvisits, ra, dec, pmra, pmdec,
     h, baseline, fields (Same as 2.2.1 above)

    cadence options for these targets (list all options,
     even though no single target will receive more than one):
      RV-6, RV-12 (plate cadences)

    Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
     pmra, pmdec, h, baseline, fields FROM mwm_rv_long
      WHERE h<12.2 AND nvisits>=6 AND
       [(ra,dec)< 3 degrees on sky
       from (11.83254, 85.251) OR (47.39417, 39.50294) OR ...
        <rest of the field centers listed below>]
    RA (Degrees) DEC(Degrees)
    11.83254    85.251
    47.39417    39.50294
    72.4084    63.603
    75.35158    22.23817
    94.8905    41.2531
    111.7472    23.6726
    114.5958    37.38599
    114.5958    21
    114.7236    17.1076
    170.3706    9.9099
    195.5556    57.06683
    200.6    -13.7
    206.4    -5.5
    207.6239    2.7267
    217.9883    8.799444
    226.0916    32.3294
    232.0711    6.8184
    234.0497    86.1432
    240.575    28.09
    248.3458    -0.5336

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
    """
    name = 'mwm_rv_long_bplates'
    category = 'science'
    cadence = None
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2526

    # peewee Model name ---> postgres table name
    # SDSS_APOGEE_AllStarMerge_r13(CatalogdbModel)--->'sdss_apogeeallstarmerge_r13'

    def build_query(self, version_id, query_region=None):
        ra = [0] * 21
        dec = [0] * 21
        ra[1] = 11.83254;  dec[1] = 85.251     # noqa: E702, E241
        ra[2] = 47.39417;  dec[2] = 39.50294   # noqa: E702, E241
        ra[3] = 72.4084;   dec[3] = 63.603     # noqa: E702, E241
        ra[4] = 75.35158;  dec[4] = 22.23817   # noqa: E702, E241
        ra[5] = 94.8905;   dec[5] = 41.2531    # noqa: E702, E241
        ra[6] = 111.7472;  dec[6] = 23.6726    # noqa: E702, E241
        ra[7] = 114.5958;  dec[7] = 37.38599   # noqa: E702, E241
        ra[8] = 114.5958;  dec[8] = 21         # noqa: E702, E241
        ra[9] = 114.7236;  dec[9] = 17.1076    # noqa: E702, E241
        ra[10] = 170.3706; dec[10] = 9.9099    # noqa: E702, E241
        ra[11] = 195.5556; dec[11] = 57.06683  # noqa: E702, E241
        ra[12] = 200.6;    dec[12] = -13.7     # noqa: E702, E241
        ra[13] = 206.4;    dec[13] = -5.5      # noqa: E702, E241
        ra[14] = 207.6239; dec[14] = 2.7267    # noqa: E702, E241
        ra[15] = 217.9883; dec[15] = 8.799444  # noqa: E702, E241
        ra[16] = 226.0916; dec[16] = 32.3294   # noqa: E702, E241
        ra[17] = 232.0711; dec[17] = 6.8184    # noqa: E702, E241
        ra[18] = 234.0497; dec[18] = 86.1432   # noqa: E702, E241
        ra[19] = 240.575;  dec[19] = 28.09     # noqa: E702, E241
        ra[20] = 248.3458; dec[20] = -0.5336   # noqa: E702, E241

        ra_dec_condition = (
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[1], dec[1], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[2], dec[2], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[3], dec[3], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[4], dec[4], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[5], dec[5], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[6], dec[6], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[7], dec[7], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[8], dec[8], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[9], dec[9], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[10], dec[10], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[11], dec[11], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[12], dec[12], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[13], dec[13], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[14], dec[14], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[15], dec[15], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[16], dec[16], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[17], dec[17], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[18], dec[18], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[19], dec[19], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[20], dec[20], 3))

# We use *mwm_rv_long_condition to unpack the tuple mwm_rv_long_condition.
# However, ra_dec_condition is not a tuple so it does not have a * in the front.
        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
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
                        ra_dec_condition,
                        SDSS_APOGEE_AllStarMerge_r13.h < 12.2,
                        SDSS_APOGEE_AllStarMerge_r13.nvisits >= 6))
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

    Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
     pmra, pmdec, h, baseline, fields FROM mwm_rv_long WHERE h<11.5

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
    """
    name = 'mwm_rv_long_fps'
    category = 'science'
    cadence = None
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2500

    # peewee Model name ---> postgres table name
    # SDSS_APOGEE_AllStarMerge_r13(CatalogdbModel)--->'sdss_apogeeallstarmerge_r13'

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
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

# Notes for above mwm_rv_short_condition:
#
# (a)Above Gaia_DR2.parallax.is_null(False) is not in the pseudo SQL
# but is in the text description as in the below line from the text:
# "with the additional requirement of having a Gaia parallax.
#
# (b) TwoMassPSC.gal_contam has a type smallint.
# Hence I compare it above with 0 instead of '000'.
# (The pseudo SQL compares gal_contam with '000')
#
# (c) There are two tables with a column called ph_qual:
#
# allwise: ph_qual       | character(4)   |
# twomass_psc: ph_qual       | character(3)     |
#
# The column twomass_psc.ph_qual has type character(3) which
# is the same as the 3 character strings in the pseudo SQL
# Hence, I have used TwoMassPSC.ph_qual above.
#
# (d) The Pseudo SQL has the below condition:
# cc_flag='000'
#
# However there is no column called cc_flag.
# There are below columns with similar names:
#
# allwise: cc_flags      | character(4)   |
# twomass_psc: cc_flg        | character(3)     |
#
# The column twomass_psc.cc_flg has type character(3) which
# is the same as the 3 character strings in the pseudo SQL
# Hence, I have used TwoMassPSC.cc_flg above.
#
#


class MWM_RV_Short_RM_Carton(BaseCarton):
    """ 2.2.2.1. Short Baseline Targets for RM Plates

    Shorthand name: mwm_rv_short_rm

    Simplified Description of selection criteria: Select from short-baseline targets
     (above) within the 3 BHM-RM fields observable from APO
     (sexagesimal coordinates given in Pseudo SQL below)

    Wiki page: Binaries and Substellar Companions (Page Under Construction)

    Additional source catalogs needed: Select from mwm_rv_short (2.2.2)

    Additional cross-matching needed: (None)

    Return columns: ra, dec, pmra, pmdec, h_m

    cadence options for these targets (list all options,
     even though no single target will receive more than one):
       mwm_rv_<nn>x8 (where <nn> = number of visits to RM fields)

    Pseudo SQL (optional): SELECT catalogid, ra, dec, pmra, pmdec,
     epoch, h_m FROM mwm_rv_short
      WHERE [(ra,dec) < 3 degrees on sky from
       (14:14:49 +53:05:00) OR (10:00:00 +02:12:00) OR (02:23:30 -04:15:00)]

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
    """
    name = 'mwm_rv_short_rm'
    category = 'science'
    cadence = None
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2571

    def setup_transaction(self):
        pass

    def build_query(self, version_id, query_region=None):

        # WHERE [(ra,dec) < 3 degrees on sky from
        # (14:14:49 +53:05:00) OR (10:00:00 +02:12:00) OR (02:23:30 -04:15:00)]

        a_ra = Angle('14h14m49s').degree
        a_dec = Angle('53d05m00s').degree

        b_ra = Angle('10h00m00s').degree
        b_dec = Angle('02d12m00s').degree

        c_ra = Angle('02h23m30s').degree
        c_dec = Angle('-04d15m00s').degree

        log.debug('Creating temporary table for radial query...')

        radial_query = (Catalog
                        .select(Catalog.catalogid)
                        .where(peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                          a_ra, a_dec, 3) |
                               peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                          b_ra, b_dec, 3) |
                               peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                          c_ra, c_dec, 3))
                        .cte('radial_query', materialized=True))

        RadialQuery = peewee.Table('mwm_rv_short_rm_temp')

        (CatalogToTIC_v8
         .select(radial_query.c.catalogid)
         .join(radial_query,
               on=(radial_query.c.catalogid == CatalogToTIC_v8.catalogid))
         .where(CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True)
         .with_cte(radial_query)
         .create_table(RadialQuery.__name__, temporary=True))

        query = (CatalogToTIC_v8
                 .select(RadialQuery.c.catalogid)
                 .join(RadialQuery,
                       on=(RadialQuery.c.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .join_from(TIC_v8, TwoMassPSC,
                            on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join_from(TIC_v8, AllWise,
                            on=(TIC_v8.allwise == AllWise.designation))
                 .where(*mwm_rv_short_condition))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(radial_query.c.ra,
                                                       radial_query.c.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class MWM_RV_Short_Bplates_Carton(BaseCarton):
    """ 3.2.2.2. Short Baseline Targets for Bright Time Plates

    Shorthand name: mwm_rv_short_bplates

    Simplified Description of selection criteria:
     Select from short-baseline targets
     (above) observable from APO with H brighter than 12.2 with
      at least 100 other targets within 3 degrees that meet the other criteria.

    Wiki page: Binaries and Substellar Companions (Page Under Construction)

    Additional source catalogs needed:
     Select from mwm_rv_short (2.2.2), but need to reference mwm_rv_long_bplates

    Additional cross-matching needed: (None)

    Return columns: ra, dec, pmra, pmdec, h_m;

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
     RV-12 (plate cadence)

    Pseudo SQL (optional): SELECT catalogid, ra, dec, pmra,
     pmdec, epoch, h_m FROM mwm_rv_short
      WHERE h_m<12.2 AND [(ra,dec)< 3 degrees on sky from
       (11.83254, 85.251) OR (47.39417, 39.50294) OR ...
       <rest of the field centers listed below>]
    RA (Degrees)    DEC(Degrees)
    11.83254    85.251
    47.39417    39.50294
    72.4084    63.603
    75.35158    22.23817
    94.8905 41.2531
    111.7472    23.6726
    114.5958    37.38599
    114.5958    21
    114.7236    17.1076
    170.3706    9.9099
    195.5556    57.06683
    200.6    -13.7
    206.4    -5.5
    207.6239    2.7267
    217.9883    8.799444
    226.0916    32.3294
    232.0711    6.8184
    234.0497    86.1432
    240.575    28.09
    248.3458    -0.5336

    Implementation:

    Non-SQL implementation:

    lead contact: Nicholas Troup

    Target selection final for v0?: No
"""
    name = 'mwm_rv_short_bplates'
    category = 'science'
    cadence = None
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2576

    def build_query(self, version_id, query_region=None):

        ra = [0] * 21
        dec = [0] * 21
        ra[1] = 11.83254;  dec[1] = 85.251     # noqa: E702, E241
        ra[2] = 47.39417;  dec[2] = 39.50294   # noqa: E702, E241
        ra[3] = 72.4084;   dec[3] = 63.603     # noqa: E702, E241
        ra[4] = 75.35158;  dec[4] = 22.23817   # noqa: E702, E241
        ra[5] = 94.8905;   dec[5] = 41.2531    # noqa: E702, E241
        ra[6] = 111.7472;  dec[6] = 23.6726    # noqa: E702, E241
        ra[7] = 114.5958;  dec[7] = 37.38599   # noqa: E702, E241
        ra[8] = 114.5958;  dec[8] = 21         # noqa: E702, E241
        ra[9] = 114.7236;  dec[9] = 17.1076    # noqa: E702, E241
        ra[10] = 170.3706; dec[10] = 9.9099    # noqa: E702, E241
        ra[11] = 195.5556; dec[11] = 57.06683  # noqa: E702, E241
        ra[12] = 200.6;    dec[12] = -13.7     # noqa: E702, E241
        ra[13] = 206.4;    dec[13] = -5.5      # noqa: E702, E241
        ra[14] = 207.6239; dec[14] = 2.7267    # noqa: E702, E241
        ra[15] = 217.9883; dec[15] = 8.799444  # noqa: E702, E241
        ra[16] = 226.0916; dec[16] = 32.3294   # noqa: E702, E241
        ra[17] = 232.0711; dec[17] = 6.8184    # noqa: E702, E241
        ra[18] = 234.0497; dec[18] = 86.1432   # noqa: E702, E241
        ra[19] = 240.575;  dec[19] = 28.09     # noqa: E702, E241
        ra[20] = 248.3458; dec[20] = -0.5336   # noqa: E702, E241

        ra_dec_condition = (
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[1], dec[1], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[2], dec[2], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[3], dec[3], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[4], dec[4], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[5], dec[5], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[6], dec[6], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[7], dec[7], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[8], dec[8], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[9], dec[9], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[10], dec[10], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[11], dec[11], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[12], dec[12], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[13], dec[13], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[14], dec[14], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[15], dec[15], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[16], dec[16], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[17], dec[17], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[18], dec[18], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[19], dec[19], 3) |
            peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec, ra[20], dec[20], 3))

# We use *mwm_rv_short_condition to unpack the tuple mwm_rv_short_condition.
# However, ra_dec_condition is not a tuple so it does not have a * in the front.
        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
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
                        ra_dec_condition,
                        TwoMassPSC.h_m < 12.2))
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
    cadence = None
    program = 'mwm_rv'
    mapper = 'MWM'
    priority = 2510

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
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
