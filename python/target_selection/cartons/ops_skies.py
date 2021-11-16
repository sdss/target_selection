#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu), Updated Nov 2021 by Tom Dwelly
# @Date: 2020-07-26
# @Filename: ops_skies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from target_selection.cartons.base import BaseCarton

from sdssdb.peewee.sdss5db.catalogdb import CatalogToSkies_v2, Skies_v2


class OPS_Sky_Boss_Best_Carton(BaseCarton):
    """Top quality skies for the BOSS spectrograph.

    Definition:
        Select sky positions from catalogdb.skies_v2 that don't have a
        nearby Gaia, Tycho2, or 2MASS source and that are
        valid in either lsdr8 or ps1dr2.
        Take care near edges of ls8 and ps1dr2 footprints by selecting
        only sky locations that have at least one neighbour (not too close)
        in their originating catalogue.

    """

    name = 'ops_sky_boss_best'
    cadence = None
    category = 'ops_sky'
    program = 'SKY'
    mapper = None
    instrument = 'BOSS'
    inertial = True
    priority = 5000

    load_magnitudes = False

    def build_query(self, version_id, query_region=None):
        pars = self.parameters

        query = (
            Skies_v2
            .select(
                CatalogToSkies_v2.catalogid,
                Skies_v2.ra,
                Skies_v2.dec,
            )
            .join(CatalogToSkies_v2)
            .where(CatalogToSkies_v2.version_id == version_id,
                   CatalogToSkies_v2.best >> True)
            .where(Skies_v2.valid_gaia >> True,
                   Skies_v2.valid_tmass >> True,
                   Skies_v2.valid_tycho2 >> True,
                   Skies_v2.valid_tmass_xsc >> True)
            .where(
                fn.COALESCE(Skies_v2.sep_neighbour_ls8, 1e30) > pars['min_sep_ls8'],
                fn.COALESCE(Skies_v2.sep_neighbour_ps1dr2, 1e30) > pars['min_sep_ps1dr2'],
                (
                    (Skies_v2.valid_ls8 >> True) &
                    (Skies_v2.sep_neighbour_ls8 < pars['max_sep_ls8'])
                ) |
                (
                    (Skies_v2.valid_ps1dr2 >> True) &
                    (Skies_v2.sep_neighbour_ps1dr2 < pars['max_sep_ps1dr2'])
                )
            )
        )

        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Skies_v2.ra,
                                                       Skies_v2.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class OPS_Sky_Boss_Fallback_Carton(BaseCarton):
    """Reasonable quality skies for the BOSS spectrograph.

    Definition:
        Select sky positions from catalogdb.skies_v2 that don't have a
        nearby Gaia, Tycho2, or 2MASS source

    """

    name = 'ops_sky_boss_fallback'
    cadence = None
    category = 'ops_sky'
    program = 'SKY'
    mapper = None
    instrument = 'BOSS'
    inertial = True
    priority = 5001

    load_magnitudes = False

    def build_query(self, version_id, query_region=None):
        pars = self.parameters

        query = (
            Skies_v2
            .select(
                CatalogToSkies_v2.catalogid,
                Skies_v2.ra,
                Skies_v2.dec,
            )
            .join(CatalogToSkies_v2)
            .where(Skies_v2.sep_neighbour_gaia > pars['min_sep_gaia'])
            .where(CatalogToSkies_v2.version_id == version_id,
                   CatalogToSkies_v2.best >> True)
            .where(Skies_v2.valid_gaia >> True,
                   Skies_v2.valid_tmass >> True,
                   Skies_v2.valid_tycho2 >> True,
                   Skies_v2.valid_tmass_xsc >> True)
        )

        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Skies_v2.ra,
                                                       Skies_v2.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
