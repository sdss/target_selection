#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu), Updated Nov 2021 by Tom Dwelly and Felipe Santana
# @Date: 2020-07-26
# @Filename: ops_skies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import CatalogToSkies_v2, Skies_v2

from target_selection.cartons.base import BaseCarton


# ############################################
# ############################################
# ############################################
# ############################################
# This module provides the following cartons in v0.5:
#  *  ops_sky_boss_best
#  *  ops_sky_boss_good
#  *  #### ops_sky_boss_fallback   <- not sure how to implement this


class OPS_Sky_Boss_Best_Carton(BaseCarton):
    """Top quality skies for the BOSS spectrograph.

    Definition:
        Select sky positions from catalogdb.skies_v2 that don't have a
        nearby Gaia, Tycho2, or 2MASS source and that are
        valid in at least lsdr8 OR ps1dr2 (and not obviously invalid in the other).
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
        a_large_value = 1e30
        pars = self.parameters

        query = (
            Skies_v2
            .select(
                CatalogToSkies_v2.catalogid,
                Skies_v2.ra,
                Skies_v2.dec,
                Skies_v2.pix_32768,
                Skies_v2.tile_32,
            )
            .join(CatalogToSkies_v2)
            .where(CatalogToSkies_v2.version_id == version_id,
                   CatalogToSkies_v2.best >> True)
            .where(Skies_v2.valid_gaia >> True,
                   Skies_v2.valid_tmass >> True,
                   Skies_v2.valid_tycho2 >> True,
                   Skies_v2.valid_tmass_xsc >> True)
            .where(
                fn.COALESCE(Skies_v2.sep_neighbour_ls8, a_large_value) > pars['min_sep_ls8'],
                fn.COALESCE(Skies_v2.sep_neighbour_ps1dr2, a_large_value) > pars['min_sep_ps1dr2'],
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


class OPS_Sky_Boss_Good_Carton(BaseCarton):
    """Reasonable quality skies for the BOSS spectrograph.
       Use these skies when there are not enough OPS_Sky_Boss_Best_Carton skies available

    Definition:
        Select sky positions from catalogdb.skies_v2 that don't have a
        nearby Gaia, Tycho2, or 2MASS source.
        Results in a catalogue that covers the vast majority of the sky,
        excluding only the Galactic bulge and inner parts of the Magellanic clouds.

    """

    name = 'ops_sky_boss_good'
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
                Skies_v2.pix_32768,
                Skies_v2.tile_32,
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


# class OPS_Sky_Boss_Fallback_Carton(BaseCarton):
#    """Unconstrained skies for the BOSS spectrograph.
#       Use these skies when there are not enough OPS_Sky_Boss_Best_Carton,
#       or OPS_Sky_Boss_Good_Carton skies available
#
#    Definition:
#        Select tile_32 locations that have few best/good skies available
#        Take all sky locations in those pixels
#
#    """
#
#    name = 'ops_sky_boss_fallback'
#    cadence = None
#    category = 'ops_sky'
#    program = 'SKY'
#    mapper = None
#    instrument = 'BOSS'
#    inertial = True
#    priority = 5002
#
#    load_magnitudes = False

    '''
    Here is the SQL to generate the sort of query I want -
       However I have no idea how to implement this in peewee
       within the BaseCarton framework

# Algorithm:
# 1) First find all nside=32 pixels with fewer than 1000 skies/pix in ops_sky_boss_good
#    Here I use the temp table, but would be better to read carton from targetdb or
#    just to re-issue the ops_sky_boss_good query
#    There are 12288 healpixels in NSIDE=32, and so we do the
#    generate_series(0,12287) step to make sure we catch them all
#
# 2) Now select sky candidates (from skies_v2) that land in those pixels,
#    and filter them with some less rigorous criteria than used for good+best cartons

DROP TABLE IF EXISTS sandbox.temp_ops_sky_boss_good_missing_pix ;

SELECT p.tile_32,COALESCE(b.nsky,0) as nsky
INTO sandbox.temp_ops_sky_boss_good_missing_pix
FROM (SELECT generate_series(0,12287) AS tile_32) AS p
     LEFT OUTER JOIN
     (SELECT tile_32,count(*) AS nsky FROM sandbox.temp_ops_sky_boss_good GROUP BY tile_32) AS b
     ON p.tile_32 = b.tile_32
     WHERE COALESCE(b.nsky,0) < 1000 ;

CREATE INDEX ON sandbox.temp_ops_sky_boss_good_missing_pix (tile_32);
ANALYZE sandbox.temp_ops_sky_boss_good_missing_pix;
DROP TABLE IF EXISTS sandbox.temp_ops_sky_boss_good_missing_pix_skies;

SELECT p.nsky,s.*
INTO sandbox.temp_ops_sky_boss_good_missing_pix_skies
FROM sandbox.temp_ops_sky_boss_good_missing_pix AS p
JOIN skies_v2 AS s
ON p.tile_32 = s.tile_32
WHERE selected_gaia is true
  AND COALESCE(sep_neighbour_gaia,1e30) > 3.0
  AND COALESCE(sep_neighbour_ps1dr2,1e30) > 3.0
  AND COALESCE(sep_neighbour_tycho2,1e30) > 15.0
  AND COALESCE(sep_neighbour_tmass,1e30) > 5.0;

    '''
class OPS_Sky_APOGEE_Best_Carton(BaseCarton):
    """First option of skies for the APOGEE spectrograph.

    Definition:
        Select sky positions from catalogdb.skies_v2 that don't have a
        nearby Gaia, Tycho2, 2MASS or 2MASS_XSC source.

    """

    name = 'ops_sky_apogee_best'
    cadence = None
    category = 'ops_sky'
    program = 'SKY'
    mapper = None
    instrument='APOGEE'
    interial=True
    priority = 5200

    load_magnitudes = False

    def build_query(self, version_id, query_region=None):

        query = (Skies_v2
                 .select(CatalogToSkies_v2.catalogid,
                         Skies_v2.ra, Skies_v2.dec,
                         Skies_v2.pix_32768, Skies_v2.tile_32)
                 .join(CatalogToSkies_v2)
                 .where(CatalogToSkies_v2.version_id == version_id,
                        CatalogToSkies_v2.best >> True)
                 .where(Skies_v2.valid_gaia >> True,
                        Skies_v2.valid_tmass >> True,
                        Skies_v2.valid_tycho2 >> True,
                        Skies_v2.valid_tmass_xsc >> True))

        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Skies_v2.ra,
                                                       Skies_v2.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class OPS_Sky_APOGEE_Good_Carton(BaseCarton):
    """Last resource skies for the APOGEE spectrograph to be used when carton 1 is not available.

    Definition:
        Select sky positions from catalogdb.skies_v2 that 
        are invalid in 2MASS but selected.

    """

    name = 'ops_sky_apogee_good'
    cadence = None
    category = 'ops_sky'
    program = 'SKY'
    mapper = None
    instrument='APOGEE'
    inertial=True
    priority = 5201

    load_magnitudes = False

    def build_query(self, version_id, query_region=None):

        query = (Skies_v2
                 .select(CatalogToSkies_v2.catalogid,
                         Skies_v2.ra, Skies_v2.dec,
                         Skies_v2.pix_32768, Skies_v2.tile_32)
                 .join(CatalogToSkies_v2)
                 .where(CatalogToSkies_v2.version_id == version_id,
                        CatalogToSkies_v2.best >> True)
                 .where(Skies_v2.selected_tmass >> True,
                        Skies_v2.valid_tmass >> False)
                 .where(Skies_v2.mag_neighbour_tmass>10.0))
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Skies_v2.ra,
                                                       Skies_v2.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


 
