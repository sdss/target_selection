#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-29
# @Filename: mwm_galactic.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             CatalogToTwoMassPSC,
                                             Gaia_DR3, TwoMassPSC)

from target_selection.cartons import BaseCarton


class MWM_Galactic_Core_Carton(BaseCarton):
    """Galactic Genesis carton.

    Definition: Selection of all IR-bright, red stars – vast majority are red
    giants; follows the density distribution of the MW (concentrated in the
    plane and bulge). Select sources brighter than H<11 AND ((G-H) > 3.5 OR
    Gaia non-detection). Approximately 5 million stars.

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
    Shorthand name: mwm_galactic_core
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

    name = 'mwm_galactic_core'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 2710
    program = 'mwm_galactic'
    mapper = 'MWM'
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra,
                         Gaia_DR3.dec,
                         Gaia_DR3.pmra,
                         Gaia_DR3.pmdec,
                         TwoMassPSC.h_m,
                         TwoMassPSC.designation,
                         Gaia_DR3.phot_g_mean_mag)
                 .join(CatalogToGaia_DR3,
                       on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
                 .join(Gaia_DR3,
                       on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .switch(Catalog)
                 .join(CatalogToTwoMassPSC,
                       on=(Catalog.catalogid == CatalogToTwoMassPSC.catalogid))
                 .join(TwoMassPSC,
                       on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        TwoMassPSC.h_m < 11,
                        (TwoMassPSC.ph_qual[1] == 'A') | (TwoMassPSC.ph_qual[1] == 'B'),
                        TwoMassPSC.gal_contam == 0,
                        TwoMassPSC.cc_flg[1] == 0,
                        (TwoMassPSC.rd_flg[1] > 0) & (TwoMassPSC.rd_flg[1] <= 3),
                        ((Gaia_DR3.phot_g_mean_mag - TwoMassPSC.h_m) > 3.5) |
                        (Gaia_DR3.phot_g_mean_mag >> None)))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
