#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-29
# @Filename: mwm_galactic.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             TIC_v8, TwoMassPSC)

from . import BaseCarton


class MWM_Galactic_Genesis_Carton(BaseCarton):
    """Galactic Genesis carton.

    Definition: Selection of all IR-bright, red stars – vast majority are red
    giants; follows the density distribution of the MW (concentrated in the
    plane and bulge). Select sources brighter than H<11 AND ((G-H) > 3.5 OR
    Gaia non-detection). Approximately 5 million stars.

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

    """

    name = 'mwm_gg_core'
    category = 'science'
    cadence = 'mwm_galactic_1x1'
    priority = 2710
    program = 'mwm_gg'
    mapper = 'MWM'

    def build_query(self, version_id, query_region=None):

        Hmag = TIC_v8.hmag
        Gmag = TIC_v8.gaiamag

        ph_qual = TwoMassPSC.ph_qual
        cc_flg = TwoMassPSC.cc_flg
        rd_flg = TwoMassPSC.rd_flg
        rd_flag_1 = peewee.fn.substr(rd_flg, 2, 1).cast('integer')

        gal_contam = TwoMassPSC.gal_contam

        h_max = self.parameters['h_max']
        g_h = self.parameters['g_h']

        gg = (TIC_v8
              .select(CatalogToTIC_v8.catalogid,
                      Hmag, Gmag,
                      ph_qual, cc_flg, rd_flg, gal_contam)
              .join(TwoMassPSC)
              .join_from(TIC_v8, CatalogToTIC_v8)
              .where(Hmag < h_max,
                     ph_qual.regexp('.(A|B).'),
                     gal_contam == 0,
                     peewee.fn.substr(cc_flg, 2, 1) == '0',
                     rd_flag_1 > 0, rd_flag_1 <= 3,
                     ((Gmag - Hmag) > g_h) | (TIC_v8.gaia >> None))
              .where(CatalogToTIC_v8.version_id == version_id,
                     CatalogToTIC_v8.best >> True))

        if query_region:
            gg = (gg
                  .join_from(CatalogToTIC_v8, Catalog)
                  .where(peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                    query_region[0],
                                                    query_region[1],
                                                    query_region[2])))

        return gg
