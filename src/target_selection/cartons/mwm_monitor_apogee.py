#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-24
# @Filename: mwm_monitor_apogee.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    SDSS_DR17_APOGEE_Allstarmerge,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee model : postgres table name
# Catalog : catalogdb.catalog
# Gaia_DR3 : catalogdb.gaia_dr3_source
# SDSS_DR17_APOGEE_Allstarmerge : catalogdb.sdss_dr17_apogee_allstarmerge
# TwoMassPSC : catalogdb.twomass_psc
#
#
# We associate each APOGEE target with a 2MASS source since the
# vast majority of APOGEE targets were selected from the 2MASS PSC.
# APOGEE_ID is the same as
# the “designation” column of the 2MASS PSC;
#
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
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
#


class MWM_monitor_apogee_Base_Carton(BaseCarton):
    """
    MWM_monitor_apogee_Base_Carton is a base carton.

    Actual cartons are implemented as subclasses of MWM_monitor_apogee_Base_Carton
    for different conditions in the WHERE clause.

    5.1.25. mwm_monitor_*apogee_*

    Shorthand name: mwm_monitor_apogee

    Existing carton code: None

    Simplified Description of selection criteria:
    These are targets selected from globular and open clusters that
    have been observed many times over many years with the apogee instrument.

    mwm_monitor_n188_apogee_long:
    SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%N188%'
    AND baseline >= 3000 AND nvisits >= 12  high priority

    mwm_monitor_n188_apogee_short:
    SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%N188%'
    AND baseline < 3000 AND baseline >= 1800 AND nvisits >= 12  lower priority

    mwm_monitor_m67_apogee_long:
    SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M67%'
    AND baseline >= 2000 AND nvisits >= 12  high priority

    mwm_monitor_m67_apogee_short:
    SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M67%'
    AND baseline < 2000 AND baseline >= 1300 AND nvisits >= 12  lower priority

    mwm_monitor_m15_apogee_long:
    SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M15%'
    AND baseline >= 1300 AND nvisits >= 12  high priority

    mwm_monitor_m15_apogee_short:
    SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M15%'
    AND baseline < 1300 AND baseline >= 900 AND nvisits >= 12  lower priority

    Gaia DR2 parameters to be converted to Gaia DR3: Nothing to be converted.

    Return columns: RA, DEC, proper motion, hmag, 2mass name.

    Metadata:

    Priority:
    High priority=1300
    Lower Priority =1302

    Cadence: bright_1x4

    Instrument: APOGEE

    can_offset=True

    Lead contact: Nathan De Lee
    """

    # SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
    # FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%N188%'
    # AND baseline >= 3000 AND nvisits >= 12  high priority

    # In the below query, we use replace() instead of ltrim() since
    # ltrim('2M20', '2M') will also trim the second 2.

    def build_query(self, version_id, query_region=None):
        query = (
            Catalog.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag,
                TwoMassPSC.pts_key,
                TwoMassPSC.designation,
                SDSS_DR17_APOGEE_Allstarmerge.apogee_id,
                SDSS_DR17_APOGEE_Allstarmerge.ra,
                SDSS_DR17_APOGEE_Allstarmerge.dec,
                SDSS_DR17_APOGEE_Allstarmerge.nvisits,
                SDSS_DR17_APOGEE_Allstarmerge.baseline,
                SDSS_DR17_APOGEE_Allstarmerge.j,
                SDSS_DR17_APOGEE_Allstarmerge.j_err,
                SDSS_DR17_APOGEE_Allstarmerge.h,
                SDSS_DR17_APOGEE_Allstarmerge.h_err,
                SDSS_DR17_APOGEE_Allstarmerge.k,
                SDSS_DR17_APOGEE_Allstarmerge.k_err,
                SDSS_DR17_APOGEE_Allstarmerge.fields,
            )
            .join(CatalogToGaia_DR3, on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .join(
                SDSS_DR17_APOGEE_Allstarmerge,
                on=(
                    TwoMassPSC.designation
                    == peewee.fn.replace(SDSS_DR17_APOGEE_Allstarmerge.apogee_id, "2M", "")
                ),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
            )
        )

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

    def post_process(self, model):
        """
        Select the brightest sources (i.e. smallest phot_g_mean_mag)
        for duplicate apogee_id (but different gaia dr3 source id)
        """

        # self.name is the carton name
        temp_table_name = "sandbox.temp_" + self.name
        self.database.execute_sql("update " + temp_table_name + " set selected = false")

        cursor = self.database.execute_sql(
            "select catalogid, apogee_id, phot_g_mean_mag from "
            + temp_table_name
            + " order by apogee_id asc, phot_g_mean_mag asc;"
        )

        output = cursor.fetchall()
        list_of_catalog_id = [0] * len(output)

        # We use the name counter instead of count
        # since count is a python built-in function.
        # For counter we use a dictionary
        # since current_apogee_id is a string.
        counter = {}

        current_target = 0
        for i in range(len(output)):
            current_apogee_id = output[i][1]
            # If two targets have same apogee id then we only keep the
            # first source. Since the sources are
            # ordered by apogee_id (ascending) and phot_g_mean_mag (ascending)
            # this means that we will keep the sources with
            # smallest phot_g_mean_mag (i.e. brightest sources).
            if current_apogee_id not in counter.keys():
                counter[current_apogee_id] = True
                list_of_catalog_id[current_target] = output[i][0]
                current_target = current_target + 1

        max_target = current_target
        for k in range(max_target + 1):
            self.database.execute_sql(
                " update "
                + temp_table_name
                + " set selected = true "
                + " where catalogid = "
                + str(list_of_catalog_id[k])
                + ";"
            )


# Below are sub classes of the above base class MWM_monitor_apogee_Base_Carton.
# Note that there is no __init__() in below sub classes since
# they use the __init__() from the BaseCarton class in base.py.


class MWM_monitor_n188_apogee_long_Carton(MWM_monitor_apogee_Base_Carton):
    name = "mwm_monitor_n188_apogee_long"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x4"
    program = "mwm_monitor"
    mapper = "MWM"
    priority = 1300
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            SDSS_DR17_APOGEE_Allstarmerge.baseline >= 3000,
            SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 12,
            SDSS_DR17_APOGEE_Allstarmerge.fields.contains("N188"),
        )

        return query

    def post_process(self, model):
        # Note that self is a parameter of post_process() above but
        # self is not a parameter of post_process() below
        # since we are using super().
        super().post_process(model)


class MWM_monitor_n188_apogee_short_Carton(MWM_monitor_apogee_Base_Carton):
    name = "mwm_monitor_n188_apogee_short"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x4"
    program = "mwm_monitor"
    mapper = "MWM"
    priority = 1302
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            SDSS_DR17_APOGEE_Allstarmerge.baseline >= 1800,
            SDSS_DR17_APOGEE_Allstarmerge.baseline < 3000,
            SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 12,
            SDSS_DR17_APOGEE_Allstarmerge.fields.contains("N188"),
        )

        return query

    def post_process(self, model):
        super().post_process(model)


class MWM_monitor_m67_apogee_long_Carton(MWM_monitor_apogee_Base_Carton):
    name = "mwm_monitor_m67_apogee_long"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x4"
    program = "mwm_monitor"
    mapper = "MWM"
    priority = 1300
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            SDSS_DR17_APOGEE_Allstarmerge.baseline >= 2000,
            SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 12,
            SDSS_DR17_APOGEE_Allstarmerge.fields.contains("M67"),
        )

        return query

    def post_process(self, model):
        super().post_process(model)


class MWM_monitor_m67_apogee_short_Carton(MWM_monitor_apogee_Base_Carton):
    name = "mwm_monitor_m67_apogee_short"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x4"
    program = "mwm_monitor"
    mapper = "MWM"
    priority = 1302
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            SDSS_DR17_APOGEE_Allstarmerge.baseline >= 1300,
            SDSS_DR17_APOGEE_Allstarmerge.baseline < 2000,
            SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 12,
            SDSS_DR17_APOGEE_Allstarmerge.fields.contains("M67"),
        )

        return query

    def post_process(self, model):
        super().post_process(model)


class MWM_monitor_m15_apogee_long_Carton(MWM_monitor_apogee_Base_Carton):
    name = "mwm_monitor_m15_apogee_long"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x4"
    program = "mwm_monitor"
    mapper = "MWM"
    priority = 1300
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            SDSS_DR17_APOGEE_Allstarmerge.baseline >= 1300,
            SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 12,
            SDSS_DR17_APOGEE_Allstarmerge.fields.contains("M15"),
        )

        return query

    def post_process(self, model):
        super().post_process(model)


class MWM_monitor_m15_apogee_short_Carton(MWM_monitor_apogee_Base_Carton):
    name = "mwm_monitor_m15_apogee_short"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x4"
    program = "mwm_monitor"
    mapper = "MWM"
    priority = 1302
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        query = query.where(
            SDSS_DR17_APOGEE_Allstarmerge.baseline >= 900,
            SDSS_DR17_APOGEE_Allstarmerge.baseline < 1300,
            SDSS_DR17_APOGEE_Allstarmerge.nvisits >= 12,
            SDSS_DR17_APOGEE_Allstarmerge.fields.contains("M15"),
        )

        return query

    def post_process(self, model):
        super().post_process(model)
