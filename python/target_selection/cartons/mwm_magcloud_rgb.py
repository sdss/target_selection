#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2023-04-28
# @Filename: mwm_magcloud_rgb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             CatalogToTwoMassPSC,
                                             Gaia_DR3, TwoMassPSC)

from target_selection.cartons import BaseCarton


class MWM_MagCloud_RGB_Base(BaseCarton):
    """MWM Magellanic clouds RGBs.

    Definition:

    Select AGB targets in the Magellanic Clouds using Gaia DR3 photometry and astrometry.
    Selection based on color cuts, proper motion cuts and parallax cuts. See wiki
    for David Nidever's code.

    """

    mapper = 'MWM'
    category = 'science'
    program = 'mwm_magcloud_agb'
    can_offset = True

    def build_query(self, version_id, query_region=None):

        # Parallax cut
        parallax_cut = (Gaia_DR3.parallax > 0) & (Gaia_DR3.parallax / Gaia_DR3.parallax_error > 5)

        # Rough cuts. Just to reduce the number of rows returned for post-process.
        colour_cut = ((Gaia_DR3.phot_g_mean_mag <= self.parameters['gaiag_rough']) &
                      (Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag >=
                       self.parameters['bp_rp_rough']))

        pm_cut = ((Gaia_DR3.pmra >= self.parameters['pmra_rough'][0]) &
                  (Gaia_DR3.pmra <= self.parameters['pmra_rough'][1]) &
                  (Gaia_DR3.pmdec >= self.parameters['pmdec_rough'][0]) &
                  (Gaia_DR3.pmdec <= self.parameters['pmdec_rough'][1]))

        astro_cut = (Gaia_DR3.ra >= self.parameters['ra_rough'][0] &
                     Gaia_DR3.ra <= self.parameters['ra_rough'][1] &
                     Gaia_DR3.dec >= self.parameters['dec_rough'][0] &
                     Gaia_DR3.dec <= self.parameters['dec_rough'][1])

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.ra,
                         Gaia_DR3.dec,
                         Gaia_DR3.l,
                         Gaia_DR3.b,
                         Gaia_DR3.pmra,
                         Gaia_DR3.pmdec,
                         Gaia_DR3.parallax,
                         Gaia_DR3.parallax_error,
                         TwoMassPSC.h_m,
                         TwoMassPSC.j_m,
                         TwoMassPSC.k_m)
                 .join(Gaia_DR3)
                 .join(CatalogToGaia_DR3, CatalogToTwoMassPSC,
                       on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        CatalogToTwoMassPSC.best >> True,
                        parallax_cut,
                        colour_cut,
                        pm_cut,
                        astro_cut))

        return query


class MWM_MagCloud_RGB_APOGEE(MWM_MagCloud_RGB_Base):
    """MWM Magellanic clouds RGBs. APOGEE carton."""

    name = 'mwm_magcloud_agb_apogee'
    instrument = 'APOGEE'
    priority = 2818

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(TwoMassPSC.h_m < self.parameters['h_lim'])

        return query


class MWM_MagCloud_RGB_BOSS(MWM_MagCloud_RGB_Base):
    """MWM Magellanic clouds RGBs. BOSS carton."""

    name = 'mwm_magcloud_agb_boss'
    instrument = 'BOSS'
    priority = 2819

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3 < self.parameters['g_lim'])

        return query
