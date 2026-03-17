#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_rv.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import math

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTIC_v8,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    SDSS_DR17_APOGEE_Allstarmerge,
    TIC_v8,
    TwoMassPSC,
)

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

# We associate each APOGEE target with a 2MASS source since the
# vast majority of APOGEE targets were selected from the 2MASS PSC.
# APOGEE_ID is essentially the same as
# the “designation” column of the 2MASS PSC;
# For sdss_dr17_apogee_allstarmerge we
# have to strip off the “2M” in the APOGEE_ID.
# However, some apogee_id values do not have 2M in the front.
# (old: for sdss_apogeeallstarmerge_r13 we also had to strip off the
#  “2M” in the APOGEE_ID.)
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# sdss5db=# select apogee_id from catalogdb.SDSS_DR17_APOGEE_Allstarmerge
# where apogee_id like '2M%';
#      apogee_id
# ---------------------
# 2M00000002+7417074
# 2M00000019-1924498
# etc.
# sdss5db=# select apogee_id from catalogdb.SDSS_DR17_APOGEE_Allstarmerge
# where apogee_id not like '2M%';
#    apogee_id
# --------------------
# 19140272-1554055
# 19155129-1617591
#
# #### start comment for old sdss_apogeeallstarmerge_r13 ###################
#
# sdss5db=# select replace(apogee_id,'2M', '') from
#  catalogdb.sdss_apogeeallstarmerge_r13 limit 2;
#       replace
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

mwm_rv_long_condition = (
    SDSS_DR17_APOGEE_Allstarmerge.h < 11.5,  # old 12.8,
    SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 6,  # old 3,
    peewee.fn.trim(SDSS_DR17_APOGEE_Allstarmerge.dist_src) == "GaiaEDR3",
    (SDSS_DR17_APOGEE_Allstarmerge.targflags % "%APOGEE_SHORT%")
    | (SDSS_DR17_APOGEE_Allstarmerge.targflags % "%APOGEE_INTERMEDIATE%")
    | (SDSS_DR17_APOGEE_Allstarmerge.targflags % "%APOGEE_LONG%")
    | (SDSS_DR17_APOGEE_Allstarmerge.targflags % "%APOGEE2_%BIN_%"),
)


# old class name MWM_RV_Long_FPS_Carton(BaseCarton):
class MWM_bin_rv_long_apogee_Carton(BaseCarton):
    """3.2.1.3. Long Baseline Targets for FPS

    Shorthand name: mwm_bin_rv_long_apogee (old name mwm_rv_long_fps)

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

    """

    # cadence must be None here so that
    # it can be set in post_process().
    # If cadence is not None here then
    # it cannot be set in post_process().
    name = "mwm_bin_rv_long_apogee"  # old name = 'mwm_rv_long_fps'
    category = "science"
    instrument = "APOGEE"
    cadence = None  # cadence is set in post_process()
    program = "mwm_bin"
    mapper = "MWM"
    priority = 1810
    can_offset = True

    # peewee Model name ---> postgres table name
    # SDSS_DR17_APOGEE_Allstarmerge(CatalogdbModel)--->'sdss_dr17_apogee_allstarmerge'
    # There is no gaia_dr3 in the below query so for v1.0 we do not have to do
    # a major modification of the v0.5 query.
    # In the below query, we use replace() instead of ltrim() since
    # ltrim('2M20', '2M') will also trim the second 2.

    def build_query(self, version_id, query_region=None):
        query = (
            Catalog.select(
                CatalogToTIC_v8.catalogid,
                SDSS_DR17_APOGEE_Allstarmerge.apogee_id,
                SDSS_DR17_APOGEE_Allstarmerge.nvisits,
                SDSS_DR17_APOGEE_Allstarmerge.ra.alias("allstarmerge_ra"),
                SDSS_DR17_APOGEE_Allstarmerge.dec.alias("allstarmerge_dec"),
                SDSS_DR17_APOGEE_Allstarmerge.pmra.alias("allstarmerge_pmra"),
                SDSS_DR17_APOGEE_Allstarmerge.pmdec.alias("allstarmerge_pmdec"),
                SDSS_DR17_APOGEE_Allstarmerge.h,
                SDSS_DR17_APOGEE_Allstarmerge.baseline,
                SDSS_DR17_APOGEE_Allstarmerge.fields,
                SDSS_DR17_APOGEE_Allstarmerge.teff_avg,  # old teff
                SDSS_DR17_APOGEE_Allstarmerge.logg_avg,  # old logg
                SDSS_DR17_APOGEE_Allstarmerge.dist,
                SDSS_DR17_APOGEE_Allstarmerge.dist_src,
            )
            .join(CatalogToTIC_v8, on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .join(
                SDSS_DR17_APOGEE_Allstarmerge,
                on=(
                    TwoMassPSC.designation
                    == peewee.fn.replace(SDSS_DR17_APOGEE_Allstarmerge.apogee_id, "2M", "")
                ),
            )
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                *mwm_rv_long_condition,
                SDSS_DR17_APOGEE_Allstarmerge.h < 11.5,
            )
        )
        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = query.where(
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
        For cadence:
        If H>10.8 then use bright_<nn>x2, otherwise use bright_<nn>x1,
        where <nn> = 3*ceiling((18-nvisits)/3)
        if <nn> is less than 6 then
            set <nn> = 6

        """

        # teff_avg and logg_avg are from SDSS_DR17_APOGEE_Allstarmerge
        # old name was teff, logg
        cursor = self.database.execute_sql(
            "select catalogid, nvisits, h, teff_avg, logg_avg from "
            + " sandbox.temp_mwm_bin_rv_long_apogee ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_nvisits = output[i][1]
            current_h = output[i][2]
            # current_teff_avg = output[i][3]
            # current_logg_avg = output[i][4]

            nn = 3 * math.ceil((18 - current_nvisits) / 3)
            if nn < 6:
                nn = 6

            if current_h > 10.8:
                current_cadence = "bright_" + str(nn) + "x2"
            else:
                current_cadence = "bright_" + str(nn) + "x1"

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_bin_rv_long_apogee "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )


# Below is from the comments for
# MWM_bin_rv_short_Base_Carton further below.
# It is put here for reference.
#
# v1.0
# AND tm.j_msigcom <= 0.1
# AND tm.h_msigcom <= 0.1
# AND tm.k_msigcom <= 0.1
# AND (   tm.ph_qual='AAA'
#      OR tm.ph_qual='AAB'
#      OR tm.ph_qual='ABA'
#      OR tm.ph_qual='BAA'
#      OR tm.ph_qual='ABB'
#      OR tm.ph_qual='BAB'
#      OR tm.ph_qual='BBA'
#      OR tm.ph_qual='BBB')
# AND tm.prox >= 6
# AND tm.cc_flg='000'
# AND tm.gal_contam='000'
# AND  (   tm.rd_flg='111'
#       OR tm.rd_flg='112'
#       OR tm.rd_flg='121'
#       OR tm.rd_flg='211'
#       OR tm.rd_flg='122'
#       OR tm.rd_flg='212'
#       OR tm.rd_flg='221'
#       OR tm.rd_flg='222')
# AND g.parallax_error/g.parallax < 0.2

mwm_rv_short_condition = (
    TwoMassPSC.j_msigcom <= 0.1,
    TwoMassPSC.h_msigcom <= 0.1,
    TwoMassPSC.k_msigcom <= 0.1,
    (
        (TwoMassPSC.ph_qual == "AAA")
        | (TwoMassPSC.ph_qual == "AAB")
        | (TwoMassPSC.ph_qual == "ABA")
        | (TwoMassPSC.ph_qual == "BAA")
        | (TwoMassPSC.ph_qual == "ABB")
        | (TwoMassPSC.ph_qual == "BAB")
        | (TwoMassPSC.ph_qual == "BBA")
        | (TwoMassPSC.ph_qual == "BBB")
    ),
    TwoMassPSC.prox >= 6,
    TwoMassPSC.cc_flg == "000",
    TwoMassPSC.gal_contam == 0,
    (
        (TwoMassPSC.rd_flg == "111")
        | (TwoMassPSC.rd_flg == "112")
        | (TwoMassPSC.rd_flg == "121")
        | (TwoMassPSC.rd_flg == "211")
        | (TwoMassPSC.rd_flg == "122")
        | (TwoMassPSC.rd_flg == "212")
        | (TwoMassPSC.rd_flg == "221")
        | (TwoMassPSC.rd_flg == "222")
    ),
)


class MWM_bin_rv_short_Base_Carton(BaseCarton):
    # old class name MWM_RV_Short_FPS_Carton(BaseCarton):
    """
    This base carton contains the first part of the below query.
    which is common to all the mmw_bin_rv_short* cartons
    The middle part (i.e. the second part) is different for the different cartons.
    The third part is not in mwm_bin_rv_short_mdwarf.
    The third part is in mwm_bin_rv_short_subgiant and mwm_bin_rv_short_rgb.
    The third part corresponds to mwm_rv_short_condition above.

    The below information is for 5.1.25. mwm_bin_rv_short_mdwarf.
    It is put here for reference for MWM_bin_rv_short_Base_Carton.
    The below information is again put further below in 5.1.25. mwm_bin_rv_short_mdwarf.

    5.1.25. mwm_bin_rv_short_mdwarf
    This carton contains the M Dwarf stars originally in mwm_bin_rv_short

    Shorthand name:  mwm_bin_rv_short_mdwarf  (Formally mwm_rv_short_fps/mwm_bin_rv_short)

    Existing carton code :
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_rv.py

    Simplified Description of selection criteria :

    SELECT catalog.catalogid, catalog.ra, catalog.dec, catalog.pmra, catalog.pmdec,
    g.ref_epoch, g.source_id, g.parallax, g.parallax_error,
    g.phot_g_mean_mag,g.phot_bp_mean_mag,
    g.phot_rp_mean_mag,g.logg_gspphot,g.teff_gspphot,tm.j_m,tm.h_m,tm.k_m
    FROM catalog
    JOIN gaia_dr3_source AS g,twomass_psc AS tm

    WHERE tm.h_m < 11.5
    AND tm.h_m > 7
    AND g.parallax>0
    AND g.phot_rp_mean_mag-tm.k_m>2 and
    AND (g.phot_rp_mean_mag-tm.k_m)*2+0.6<tm.h_m-5*log10(1000/g.parallax)+5

    Gaia DR2 parameters to be converted to Gaia DR3 :
    coordinates, proper motions, parallax

    Return columns (again ok to copy from v0.5 version) :  (See pseudo-SQL above)

    Metadata:

    can_offset = False
    Cadence options: bright_18x1
    priority: priority = 2515
    Lead contact:  Nicholas Troup Nathan De Lee
    """

    # Below query does not use SDSS_DR17_APOGEE_Allstarmerge.
    # There is gaia_dr3 in the below query so we have done
    # major modification of the old v0.5 query.

    def build_query(self, version_id, query_region=None):
        query = (
            Catalog.select(
                Catalog.catalogid,
                Catalog.ra,
                Catalog.dec,
                Catalog.pmra,
                Catalog.pmdec,
                Gaia_DR3.ref_epoch,
                Gaia_DR3.source_id,
                Gaia_DR3.parallax,
                Gaia_DR3.parallax_error,
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag,
                Gaia_DR3.logg_gspphot,
                Gaia_DR3.teff_gspphot,
                TwoMassPSC.j_m,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.pmra.alias("gaia_dr3_pmra"),
                Gaia_DR3.pmdec.alias("gaia_dr3_pmdec"),
            )
            .join(CatalogToGaia_DR3, on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(Catalog)
            .join(CatalogToTwoMassPSC, on=(Catalog.catalogid == CatalogToTwoMassPSC.catalogid))
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                TwoMassPSC.ext_key >> None,
                Gaia_DR3.parallax.is_null(False),
                Gaia_DR3.parallax > 0,
            )
        )

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )
        return query


class MWM_bin_rv_short_mdwarf_Base_Carton(MWM_bin_rv_short_Base_Carton):
    """5.1.25. mwm_bin_rv_short_mdwarf
    This base carton contains the M Dwarf stars originally in mwm_bin_rv_short

    Shorthand name:  mwm_bin_rv_short_mdwarf  (Formally mwm_rv_short_fps/mwm_bin_rv_short)

    Existing carton code :
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_rv.py

    Simplified Description of selection criteria :

    SELECT catalog.catalogid, catalog.ra, catalog.dec, catalog.pmra, catalog.pmdec,
    g.ref_epoch, g.source_id, g.parallax, g.parallax_error,
    g.phot_g_mean_mag,g.phot_bp_mean_mag,
    g.phot_rp_mean_mag,g.logg_gspphot,g.teff_gspphot,tm.j_m,tm.h_m,tm.k_m
    FROM catalog
    JOIN gaia_dr3_source AS g,twomass_psc AS tm

    WHERE tm.h_m < 11.5
    AND tm.h_m > 7
    AND g.parallax>0
    AND g.phot_rp_mean_mag-tm.k_m>2 and
    AND (g.phot_rp_mean_mag-tm.k_m)*2+0.6<tm.h_m-5*log10(1000/g.parallax)+5

    The mwm_bin_rv_short_mdwarf carton has the above query with only two parts.
    The mwm_bin_rv_short_subgiant and mwm_bin_rv_short_rgb cartons
    have queries with three parts
    The third part corresponds to mwm_rv_short_condition above.

    Gaia DR2 parameters to be converted to Gaia DR3 :
    coordinates, proper motions, parallax

    Return columns (again ok to copy from v0.5 version) : (See pseudo-SQL above)

    Metadata:

    can_offset = False
    Cadence options: bright_18x1
    priority: priority = 2515
    Lead contact:  Nicholas Troup Nathan De Lee
    """

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            TwoMassPSC.h_m > 7,
            TwoMassPSC.h_m < 11.5,
            Gaia_DR3.phot_rp_mean_mag - TwoMassPSC.k_m > 2,
            ((Gaia_DR3.phot_rp_mean_mag - TwoMassPSC.k_m) * 2 + 0.6)
            < (TwoMassPSC.h_m - 5 * peewee.fn.log10(1000 / Gaia_DR3.parallax) + 5),
        )
        return query


class MWM_bin_rv_short_mdwarf_apogee_18epoch(MWM_bin_rv_short_mdwarf_Base_Carton):
    name = "mwm_bin_rv_short_mdwarf_apogee_18epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_18x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 1305
    can_offset = False


class MWM_bin_rv_short_mdwarf_apogee_12epoch(MWM_bin_rv_short_mdwarf_Base_Carton):
    name = "mwm_bin_rv_short_mdwarf_apogee_12epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_12x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 1306
    can_offset = False


class MWM_bin_rv_short_mdwarf_apogee_08epoch(MWM_bin_rv_short_mdwarf_Base_Carton):
    name = "mwm_bin_rv_short_mdwarf_apogee_08epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_8x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 1307
    can_offset = False


class MWM_bin_rv_short_subgiant_apogee_Carton(MWM_bin_rv_short_Base_Carton):
    """5.1.26. mwm_bin_rv_short_subgiant_apogee
    This carton contains subgiant stars and lower red giant branch stars
    originally in mwm_bin_rv_short

    Shorthand name:  mwm_bin_rv_short_subgiant  (Formally mwm_rv_short_fps/mwm_bin_rv_short)

    Existing carton code :
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_rv.py

    Simplified Description of selection criteria :

    SELECT catalog.catalogid, catalog.ra, catalog.dec, catalog.pmra, catalog.pmdec,
    g.ref_epoch, g.source_id, g.parallax, g.parallax_error,
    g.phot_g_mean_mag,g.phot_bp_mean_mag,
    g.phot_rp_mean_mag,g.logg_gspphot,g.teff_gspphot,tm.j_m,tm.h_m,tm.k_m
    FROM catalog
    JOIN gaia_dr3_source AS g,twomass_psc AS tm

    WHERE tm.h_m <10.8
    AND tm.h_m > 7
    AND g.logg_gspphot > 3.0
    AND g.logg_gspphot < 4.1
    AND g.teff_gspphot < 7000
    AND g.teff_gspphot > 4500

    AND (tm.j_msigcom <= 0.1 AND tm.h_msigcom<=0.1 AND tm.k_msigcom <= 0.1)
    AND(tm.ph_qual= 'AAA' OR tm.ph_qual='AAB' OR tm.ph_qual='ABA'
        OR tm.ph_qual='BAA' OR tm.ph_qual= 'ABB' OR tm.ph_qual='BAB'
        OR tm.ph_qual='BBA' OR tm.ph_qual ='BBB')
    AND tm.prox >= 6  AND tm.cc_flg='000' AND tm.gal_contam='000'
    AND  (tm.rd_flg='111' OR tm.rd_flg='112'  OR tm.rd_flg='121'
          OR tm.rd_flg='211'   OR tm.rd_flg='122'  OR tm.rd_flg='212'
          OR tm.rd_flg='221'  OR tm.rd_flg='222')
    AND g.parallax_error/g.parallax < 0.2

    Gaia DR2 parameters to be converted to Gaia DR3 :
    coordinates, proper motions, parallax

    Return columns (again ok to copy from v0.5 version) : (See pseudo-SQL above)

    Metadata:

    can_offset = False
    Cadence options: bright_18x1
    priority: priority = 2525
    Lead contact:  Nicholas Troup Nathan De Lee
    """

    name = "mwm_bin_rv_short_subgiant_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_18x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2525
    can_offset = False

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            TwoMassPSC.h_m > 7,
            TwoMassPSC.h_m < 10.8,
            Gaia_DR3.logg_gspphot > 3.0,
            Gaia_DR3.logg_gspphot < 4.1,
            Gaia_DR3.teff_gspphot < 7000,
            Gaia_DR3.teff_gspphot > 4500,
            Gaia_DR3.parallax_error / Gaia_DR3.parallax < 0.2,
            *mwm_rv_short_condition,
        )
        return query


class MWM_bin_rv_short_rgb_apogee_Carton(MWM_bin_rv_short_Base_Carton):
    """5.1.27. mwm_bin_rv_short_rgb_apogee
    This carton contains the red clump and higher red giant stars
    originally in mwm_bin_rv_short

    Shorthand name:  mwm_bin_rv_short_rgb  (Formally mwm_rv_short_fps/mwm_bin_rv_short)

    Existing carton code :
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_rv.py

    Simplified Description of selection criteria :

    SELECT catalog.catalogid, catalog.ra, catalog.dec, catalog.pmra, catalog.pmdec,
    g.ref_epoch, g.source_id, g.parallax, g.parallax_error,
    g.phot_g_mean_mag,g.phot_bp_mean_mag,
    g.phot_rp_mean_mag,g.logg_gspphot,g.teff_gspphot,tm.j_m,tm.h_m,tm.k_m
    FROM catalog
    JOIN gaia_dr3_source AS g,twomass_psc AS tm

    WHERE tm.h_m<10.8
    AND tm.h_m > 7
    AND g.logg_gspphot > 1.2
    AND g.logg_gspphot <= 3.0
    AND g.teff_gspphot < 5500

    AND (tm.j_msigcom <= 0.1 AND tm.h_msigcom<=0.1 AND tm.k_msigcom <= 0.1)
    AND(tm.ph_qual= 'AAA' OR tm.ph_qual='AAB' OR tm.ph_qual='ABA'
        OR tm.ph_qual='BAA' OR tm.ph_qual= 'ABB' OR tm.ph_qual='BAB'
        OR tm.ph_qual='BBA' OR tm.ph_qual ='BBB')
    AND tm.prox >= 6  AND tm.cc_flg='000' AND tm.gal_contam='000'
    AND  (tm.rd_flg='111' OR tm.rd_flg='112'  OR tm.rd_flg='121'
          OR tm.rd_flg='211'   OR tm.rd_flg='122'  OR tm.rd_flg='212'
          OR tm.rd_flg='221'  OR tm.rd_flg='222')
    AND g.parallax_error/g.parallax < 0.2

    Gaia DR2 parameters to be converted to Gaia DR3 :
    coordinates, proper motions, parallax

    Return columns (again ok to copy from v0.5 version) : (See pseudo-SQL above)

    Metadata:

    can_offset = False
    Cadence options: bright_18x1
    priority: priority = 2535
    Lead contact:  Nicholas Troup Nathan De Lee
    """

    name = "mwm_bin_rv_short_rgb_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_18x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2535
    can_offset = False

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            TwoMassPSC.h_m > 7,
            TwoMassPSC.h_m < 10.8,
            Gaia_DR3.logg_gspphot > 1.2,
            Gaia_DR3.logg_gspphot <= 3.0,
            Gaia_DR3.teff_gspphot < 5500,
            Gaia_DR3.parallax_error / Gaia_DR3.parallax < 0.2,
            *mwm_rv_short_condition,
        )
        return query


class MWM_bin_rv_short_rgb_apogee_08epoch_Carton(MWM_bin_rv_short_rgb_apogee_Carton):
    """mwm_bin_rv_short_rgb_apogee_08epoch -
    same as mwm_bin_rv_short_rgb_apogee but
    with cadence=bright_8x1 and priority=2537
    """

    name = "mwm_bin_rv_short_rgb_apogee_08epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_8x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2537
    can_offset = True


class MWM_bin_rv_short_rgb_apogee_12epoch_Carton(MWM_bin_rv_short_rgb_apogee_Carton):
    """mwm_bin_rv_short_rgb_apogee_12epoch -
    same as mwm_bin_rv_short_rgb_apogee but
    with cadence=bright_12x1 and priority=2536
    """

    name = "mwm_bin_rv_short_rgb_apogee_12epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_12x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2536
    can_offset = True


class MWM_bin_rv_short_subgiant_apogee_08epoch_Carton(MWM_bin_rv_short_subgiant_apogee_Carton):
    """mwm_bin_rv_short_subgiant_apogee_08epoch -
    same as mwm_bin_rv_short_subgiant_apogee but
    with cadence=bright_8x1 and priority=2527
    """

    name = "mwm_bin_rv_short_subgiant_apogee_08epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_8x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2527
    can_offset = True


class MWM_bin_rv_short_subgiant_apogee_12epoch_Carton(MWM_bin_rv_short_subgiant_apogee_Carton):
    """mwm_bin_rv_short_subgiant_apogee_12epoch -
    same as mwm_bin_rv_short_subgiant_apogee but
    with cadence=bright_12x1 and priority=2526
    """

    name = "mwm_bin_rv_short_subgiant_apogee_12epoch"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_12x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2526
    can_offset = True
