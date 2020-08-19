#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-07-26
# @Filename: ops_skies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import CatalogToSkies_v1, Skies_v1

from . import BaseCarton


class OPS_BOSS_Sky_Carton(BaseCarton):
    """Skies for the BOSS spectrograph.

    Definition:
        Select sky positions from catalogdb.skies_v1 that don't have a
        nearby Gaia, LS8, Tycho2, or 2MASS source.

    """

    name = 'ops_sky_boss'
    cadence = None
    category = 'ops_sky'
    program = 'SKY'
    mapper = None
    priority = 5000

    load_magnitudes = False

    def build_query(self, version_id, query_region):

        min_separation = 10 + ((12. - Skies_v1.mag_neighbour_gaia) / 0.2)

        query = (Skies_v1
                 .select(CatalogToSkies_v1.catalogid,
                         Skies_v1.ra,
                         Skies_v1.dec,
                         Skies_v1.pix_32768,
                         Skies_v1.tile_32)
                 .join(CatalogToSkies_v1)
                 .where(Skies_v1.sep_neighbour_gaia > min_separation)
                 .where(CatalogToSkies_v1.version_id == version_id,
                        CatalogToSkies_v1.best >> True)
                 .where(Skies_v1.gaia_sky >> True,
                        Skies_v1.ls8_sky >> True,
                        Skies_v1.tmass_sky >> True,
                        Skies_v1.tycho2_sky >> True,
                        Skies_v1.tmass_xsc_sky >> True))

        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Skies_v1.ra,
                                                       Skies_v1.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class OPS_APOGEE_Sky_Carton(BaseCarton):
    """Skies for the APOGEE spectrograph.

    Definition:
        Select sky positions from catalogdb.skies_v1 that don't have a
        nearby Gaia, Tycho2, or 2MASS source.

    """

    name = 'ops_sky_apogee'
    cadence = None
    category = 'ops_sky'
    program = 'SKY'
    mapper = None
    priority = 5200

    load_magnitudes = False

    def build_query(self, version_id, query_region):

        query = (Skies_v1
                 .select(CatalogToSkies_v1.catalogid,
                         Skies_v1.ra, Skies_v1.dec,
                         Skies_v1.pix_32768, Skies_v1.tile_32)
                 .join(CatalogToSkies_v1)
                 .where(CatalogToSkies_v1.version_id == version_id,
                        CatalogToSkies_v1.best >> True)
                 .where(Skies_v1.gaia_sky >> True,
                        Skies_v1.tmass_sky >> True,
                        Skies_v1.tycho2_sky >> True,
                        Skies_v1.tmass_xsc_sky >> True))

        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Skies_v1.ra,
                                                       Skies_v1.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
