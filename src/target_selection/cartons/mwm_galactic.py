#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta
# @Date: 2023-02-28
# @Filename: mwm_galactic.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import random

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


# The below MWM_Galactic_Core_apogee_Carton carton is obsolete.
# Do not run this carton for v1.0.
# However, do not delete this carton since
# this carton is referred to by the MWM_Galactic_Core_Dist_apogee_Carton
class MWM_Galactic_Core_apogee_Carton(BaseCarton):
    """Galactic Genesis carton.

    Definition: Selection of all IR-bright, red stars - vast majority are red
    giants; follows the density distribution of the MW (concentrated in the
    plane and bulge). Select sources brighter than H<11 AND ((G-H) > 3.5 OR
    Gaia non-detection).

    v0.5
    Pseudo-query:

        SELECT catalog.catalogid, catalog.ra, catalog.ra, catalog.dec,
            catalog.pmra, catalog.pmdec, catalog.epoch, twomass_psc.h_m,
            gaia_dr2_source.phot_g_mean_mag
        FROM catalog
        JOIN catalog_to_TIC, TIC, gaia_dr2_source, gaia_clean, twomass_psc
        WHERE h_m < 11
            AND (ph_qual[1] = 'A' OR ph_qual[1] = 'B')
            AND gal_contam=0
            AND cc_flg[1]=0
            AND (rd_flg[1] > 0 AND rd_flg[1] <= 3)
            AND [(phot_g_mean_mag-h_m) > 3.5  OR NOT_EXISTS(phot_g_mean_mag)]

    v1.0
    Shorthand name: mwm_galactic_core_apogee
    Existing carton code:
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_galactic.py

    Simplified Description of selection criteria:
    Selection of all IR-bright, red stars – vast majority are red giants;
     follows the density distribution of the MW (concentrated in the plane and bulge).
     Select sources brighter than H<11 AND ((G-H) > 3.5 OR Gaia non-detection).
    Gaia DR2 parameters to be converted to Gaia DR3: position information, phot_g_mean_mag,

    Pseudo-code: Same as v0.5

    SELECT catalog.catalogid, gaia_dr3.source_id, gaia_dr3.ra,  gaia_dr3.dec,
        gaia_dr3.pmra, gaia_dr3.pmdec,  twomass_psc.h_m, twomass_psc.designation,
        gaia_dr3.phot_g_mean_mag
    FROM catalog JOIN catalog_to_twomass_psc, catalog_to_gaia_dr3,  gaia_dr3,  twomass_psc
    WHERE
    h_m < 11 AND
    (ph_qual[1] = 'A' OR ph_qual[1] = 'B') AND
    gal_contam=0 AND
    cc_flg[1]=0 AND
    (rd_flg[1] > 0 AND rd_flg[1] <= 3)
    AND [(phot_g_mean_mag-h_m) > 3.5 OR NOT_EXISTS(phot_g_mean_mag)]

    Joins: TwoMassPSC, catalog_to_twomass_psc, catalog_to_gaia_dr3, gaia_dr3

    Return columns: Gaia E/DR3 source_id, ra, dec, parallax,
     Gaia E/DR3 G mag, Gaia E/DR3 BP mag, Gaia E/DR3 RP mag, J mag,
      H mag K mag, catalogid, 2MASS designation, 2MASS ph_qual, cc_flg, rd_flg, gal_contam

    Metadata:
    unchanged
    can_offset=TRUE
    Lead contact: Jonathan Bird

    """

    name = "mwm_galactic_core_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2710
    program = "mwm_galactic"
    mapper = "MWM"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToTwoMassPSC.select(
                CatalogToTwoMassPSC.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra,
                Gaia_DR3.dec,
                Gaia_DR3.pmra,
                Gaia_DR3.pmdec,
                Gaia_DR3.phot_g_mean_mag,
                TwoMassPSC.h_m,
                TwoMassPSC.designation,
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
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
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                TwoMassPSC.h_m < 11,
                (peewee.fn.substr(TwoMassPSC.ph_qual, 2, 1) == "A")
                | (peewee.fn.substr(TwoMassPSC.ph_qual, 2, 1) == "B"),
                TwoMassPSC.gal_contam == 0,
                peewee.fn.substr(TwoMassPSC.cc_flg, 2, 1) == "0",
                (peewee.fn.substr(TwoMassPSC.rd_flg, 2, 1).cast("integer") > 0)
                & (peewee.fn.substr(TwoMassPSC.rd_flg, 2, 1).cast("integer") <= 3),
                ((Gaia_DR3.phot_g_mean_mag - TwoMassPSC.h_m) > 3.5)
                | (Gaia_DR3.phot_g_mean_mag >> None),
            )
        )
        # above condition (Gaia_DR3.phot_g_mean_mag >> None) ensures that
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


class MWM_Galactic_Core_Dist_apogee_Carton(BaseCarton):
    """Galactic Genesis carton.

    Definition: Selection of all IR-bright, red stars - vast majority are red
    giants; follows the density distribution of the MW (concentrated in the
    plane and bulge). Select sources brighter than H<11 AND ((G-H) > 5 OR
    Gaia non-detection).

    v0.5 only had mwm_galactic_core
    (there was no mwm_galactic_core_dist_apogee in v0.5)

    v1.0
    Shorthand name: mwm_galactic_core_dist_apogee
    This is a new carton. Code should be put in mwm_galactic.py.
    https://github.com/sdss/target_selection/blob/main/python/target_selection/cartons/mwm_galactic.py

    Simplified Description of selection criteria:
    Selection of all IR-bright, red stars – vast majority are red giants;
     follows the density distribution of the MW (concentrated in the plane and bulge).
     Select sources brighter than H<11 AND ((G-H) > 5 OR Gaia non-detection).
    Gaia DR2 parameters to be converted to Gaia DR3: position information, phot_g_mean_mag,

    Pseudo-code: Same as v1.0 mwm_galactic_core.

    SELECT catalog.catalogid, gaia_dr3.source_id, gaia_dr3.ra,  gaia_dr3.dec,
        gaia_dr3.pmra, gaia_dr3.pmdec,  twomass_psc.h_m, twomass_psc.designation,
        gaia_dr3.phot_g_mean_mag
    FROM catalog JOIN catalog_to_twomass_psc, catalog_to_gaia_dr3,  gaia_dr3,  twomass_psc
    WHERE
    h_m < 11 AND
    (ph_qual[1] = 'A' OR ph_qual[1] = 'B') AND
    gal_contam=0 AND
    cc_flg[1]=0 AND
    (rd_flg[1] > 0 AND rd_flg[1] <= 3)

    Above  is same as for v1.0 mwm_galactic_core_apogee.
    Below is additional for v1.0 mwm_galactic_core_dist_apogee.
    The below condition on g.phot_g_mean_mag is different from the condition
    in the mwm_galactic_core carton.

    tmass.h_m < 11 AND
    (

    ( gaia.parallax < 100 * power(10, 0.2 * (-2.5 - tmass.h_m)) ) OR
    (g.phot_g_mean_mag - tmass.h_m > 5) OR
    (NO gaia.phot_g_mean_mag)
    )

    Joins: TwoMassPSC, catalog_to_twomass_psc, catalog_to_gaia_dr3, gaia_dr3

    Return columns: Gaia E/DR3 source_id, ra, dec, parallax,
     Gaia E/DR3 G mag, Gaia E/DR3 BP mag, Gaia E/DR3 RP mag, J mag,
      H mag K mag, catalogid, 2MASS designation, 2MASS ph_qual, cc_flg, rd_flg, gal_contam

    Metadata:
    unchanged
    can_offset=TRUE
    Lead contact: Jonathan Bird

    """

    name = "mwm_galactic_core_dist_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2710
    program = "mwm_galactic"
    mapper = "MWM"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToTwoMassPSC.select(
                CatalogToTwoMassPSC.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra,
                Gaia_DR3.dec,
                Gaia_DR3.pmra,
                Gaia_DR3.pmdec,
                Gaia_DR3.parallax,
                Gaia_DR3.phot_g_mean_mag,
                TwoMassPSC.h_m,
                TwoMassPSC.designation,
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
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
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                TwoMassPSC.h_m < 11,
                (peewee.fn.substr(TwoMassPSC.ph_qual, 2, 1) == "A")
                | (peewee.fn.substr(TwoMassPSC.ph_qual, 2, 1) == "B"),
                TwoMassPSC.gal_contam == 0,
                peewee.fn.substr(TwoMassPSC.cc_flg, 2, 1) == "0",
                (peewee.fn.substr(TwoMassPSC.rd_flg, 2, 1).cast("integer") > 0)
                & (peewee.fn.substr(TwoMassPSC.rd_flg, 2, 1).cast("integer") <= 3),
                (Gaia_DR3.parallax < 100 * peewee.fn.power(10, 0.2 * (-2.5 - TwoMassPSC.h_m)))
                | ((Gaia_DR3.phot_g_mean_mag - TwoMassPSC.h_m) > 5)
                | (Gaia_DR3.phot_g_mean_mag >> None),
            )
        )

        # above condition (Gaia_DR3.phot_g_mean_mag >> None) ensures that
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


class MWM_Galactic_Core_Dist_apogee_sparse_Carton(MWM_Galactic_Core_Dist_apogee_Carton):
    """mwm_galactic_core_dist_apogee_sparse
    Short description: Sparse-sampled Galactic Genesis sample.
    First do an initial selection as in mwm_galactic_core_dist_apogee and order by catalogid.
    Then randomly select exactly 2/3 of the sources using a fixed random seed.
    Metadata:
    Priority: 2710
    Cadence: bright_1x1
    Instrument: APOGEE
    can_offset = True
    Lead contact: Michael Blanton
    """

    name = "mwm_galactic_core_dist_apogee_sparse"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2710
    program = "mwm_galactic"
    mapper = "MWM"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        return query

    def post_process(self, model):
        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_galactic_core_dist_apogee_sparse " + "set selected = false;"
        )

        cursor = self.database.execute_sql(
            "select count(1) from " + " sandbox.temp_mwm_galactic_core_dist_apogee_sparse ;"
        )

        output = cursor.fetchall()

        total_num_rows = output[0][0]

        # This selected_fraction must be the same as in
        # mwm_galactic_core_dist_apogee_extra
        selected_fraction = 2.0 / 3.0
        selected_num_rows = int(selected_fraction * total_num_rows)

        b = [True] * selected_num_rows
        c = [False] * (total_num_rows - selected_num_rows)
        is_selected = b + c

        # This random seed must be the same as in
        # mwm_galactic_core_dist_apogee_extra
        random.seed(6789)
        random.shuffle(is_selected)

        # The below "order by catalogid" ensures that the random selection
        # further below gives the same result every time we run this carton.
        cursor = self.database.execute_sql(
            "select catalogid from "
            + " sandbox.temp_mwm_galactic_core_dist_apogee_sparse "
            + " order by catalogid;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]

            # This condition is opposite of the condition in
            # mwm_galactic_core_dist_apogee_extra
            if is_selected[i] is True:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_galactic_core_dist_apogee_sparse "
                    + " set selected = true "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )


class MWM_Galactic_Core_Dist_apogee_extra_Carton(MWM_Galactic_Core_Dist_apogee_Carton):
    """mwm_galactic_core_dist_apogee_extra
    Short description: Complement of above carton mwm_galactic_core_dist_apogee_sparse.
    First do an initial selection as in mwm_galactic_core_dist_apogee and order by catalogid.
    Then select those targets which were NOT selected in
    above carton mwm_galactic_core_dist_apogee_sparse.
    Metadata:
    Priority: 2711
    Cadence: bright_1x1
    Instrument: APOGEE
    can_offset = True
    Lead contact: Michael Blanton
    """

    name = "mwm_galactic_core_dist_apogee_extra"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2711
    program = "mwm_galactic"
    mapper = "MWM"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        return query

    def post_process(self, model):
        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_galactic_core_dist_apogee_extra " + "set selected = false;"
        )

        cursor = self.database.execute_sql(
            "select count(1) from " + " sandbox.temp_mwm_galactic_core_dist_apogee_extra ;"
        )

        output = cursor.fetchall()

        total_num_rows = output[0][0]

        # This selected_fraction must be the same as in
        # mwm_galactic_core_dist_apogee_sparse
        selected_fraction = 2.0 / 3.0
        selected_num_rows = int(selected_fraction * total_num_rows)

        b = [True] * selected_num_rows
        c = [False] * (total_num_rows - selected_num_rows)
        is_selected = b + c

        # This random seed must be the same as in
        # mwm_galactic_core_dist_apogee_sparse
        random.seed(6789)
        random.shuffle(is_selected)

        # The below "order by catalogid" ensures that the random selection
        # further below gives the same result every time we run this carton.
        cursor = self.database.execute_sql(
            "select catalogid from "
            + " sandbox.temp_mwm_galactic_core_dist_apogee_extra "
            + " order by catalogid;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]

            # select those sources which were not selected in
            # mwm_galactic_core_dist_apogee_sparse.
            # Hence, this condition is opposite of the condition in
            # mwm_galactic_core_dist_apogee_sparse
            if is_selected[i] is False:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_galactic_core_dist_apogee_extra "
                    + " set selected = true "
                    + " where catalogid = "
                    + str(current_catalogid)
                    + ";"
                )
