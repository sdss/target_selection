#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta
# @Date: 2023-02-28
# @Filename: mwm_galactic.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

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
