#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
#          Updated for v1 by dwelly@mpe.mpg.de
# @Date: 2021-01-12
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    EROSITASupersetv1Compactobjects,
    EROSITASupersetv1Stars,
    Gaia_DR3,
)

from target_selection.cartons import BaseCarton


# from target_selection.exceptions import TargetSelectionError

# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# For example
# peewee Model name ---> postgres table name
# AllWise ---> 'catalogdb.allwise'
# EROSITASupersetAGN --->'catalogdb.erosita_superset_agn'
# EROSITASupersetClusters ---> 'catalogdb.erosita_superset_clusters'
# EROSITASupersetStars ---> 'catalogdb.erosita_superset_stars'
# EROSITASupersetCompactobjects ---> 'catalogdb.erosita_superset_compactobjects'
# Gaia_DR2 --->'catalogdb.gaia_dr2_source'
# TwoMassPSC --->'catalogdb.twomass_psc'


class MWM_EROSITA_Stars_Boss_Carton(BaseCarton):
    """MWM eROSITA Stars
    Owner: Lead by MWM (with assistance from BHM?)?? Jennifer Johnson, Tom Dwelly

    Shorthand name: mwm_erosita_stars_boss

    What is it?: Optical counterparts to eROSITA sources with
    a high likelihood to be stars. Stellar identification are
    based on Bayesian scheme using the Gaia edr3 and Tycho catalogs.
    Those objects shall be observed to increase the legacy value
    of eROSITA / SDSS; possible science questions relate to
    Halpha emission although EW are expected to be around 1 Angstrom or less.

    Simplified Description of selection criteria:
    "Select all objects with xmatch_metric>0.0 and
    xmatch_method = 'HamStar' and xmatch_version = 'v1.1.1'
    from catalogdb.erosita_superset_v1_stars that are
    in a magnitude range amenable for SDSS-V BOSS spectroscopy
    (suitability estimated from lo < Gmag < hi)"

    Wiki page: N/A

    Additional source catalogs needed: None

    Additional cross-matching needed: No

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    boss_bright_1xN or boss_dark_1xN
    (moon phase and number of exposures depends on target magnitude).
    No detailed time spacing required.

    Pseudo SQL (optional):

    SELECT c.catalogid, c.ra, c.dec,
    gaia.phot_g_mean_mag, gaia.phot_bp_mean_mag,
    gaia.phot_rp_mean_mag, gaia.parallax, gaia.pmra, gaia.pmdec,
    estars.target_priority, estars.xmatch_metric, estars.ero_flux,
    estars.ero_ra, estars.ero_dec, estars.opt_ra, estars.opt_dec
     FROM catalog c
     INNER JOIN catalog_to_gaia_dr3_source c2g3 USING (catalogid)
     INNER JOIN gaia_dr3_source g3 ON g3.source_id = c2g3.target_id
     INNER JOIN erosita_superset_v1_stars estars ON estars.gaia_dr3_source_id = g3.source_id
    WHERE (c2g3.version_id = 21) AND /* control version! */
     (c2g3.best is true) AND /* and enforce unique-ish crossmatch */
     estars.xmatch_method = 'HamStar' AND
     estars.xmatch_version = 'v1.1.1' AND
     estars.xmatch_metric > 0.0 AND
     g3.source_phot_g_mean_mag BETWEEN 13.0 and 22.0
    ;

    See cadence logic in post_process() below.

    xmatch_metric now spans the range [0.2:1.0]

    In the query below we use xmatch_metric > 0.0 - i.e. no cut on this param
    """

    # name = 'mwm_erosita_stars'
    name = "mwm_erosita_stars_boss"
    category = "science"
    instrument = "BOSS"
    program = "mwm_erosita"
    mapper = "MWM"
    can_offset = True

    faintest_cadence = "dark_flexible_2x1"
    # faintest_priority = 1910
    flat_priority = 1920

    def build_query(self, version_id, query_region=None):
        bright_bright_limit = 13
        cadence = peewee.Case(
            None,
            (
                (Gaia_DR3.phot_g_mean_mag < 17.0, "bright_1x1"),
                (Gaia_DR3.phot_g_mean_mag < 18.0, "dark_flexible_2x1"),
            ),
            self.faintest_cadence,
        ).cast("text")

        # priority = (
        #     EROSITASupersetv1Stars.xmatch_flags - 1 +
        #     peewee.Case(
        #         None,
        #         (
        #             (Gaia_DR3.phot_g_mean_mag < 17.0, 2400),
        #         ),
        #         self.faintest_priority
        #     ).cast('integer')
        # )
        priority = peewee.Value(self.flat_priority).cast("integer")

        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag.alias("gaia_g"),
                Gaia_DR3.phot_bp_mean_mag.alias("bp"),
                Gaia_DR3.phot_rp_mean_mag.alias("rp"),
                Gaia_DR3.parallax,
                Gaia_DR3.pmra,
                Gaia_DR3.pmdec,
                cadence.alias("cadence"),
                priority.alias("priority"),
                EROSITASupersetv1Stars.target_priority,
                EROSITASupersetv1Stars.ero_detuid,
                EROSITASupersetv1Stars.xmatch_metric,
                EROSITASupersetv1Stars.xmatch_flags,
                EROSITASupersetv1Stars.ero_flux,
                EROSITASupersetv1Stars.ero_ra,
                EROSITASupersetv1Stars.ero_dec,
                EROSITASupersetv1Stars.opt_ra,
                EROSITASupersetv1Stars.opt_dec,
            )
            .distinct(CatalogToGaia_DR3.catalogid)
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(
                EROSITASupersetv1Stars,
                on=(Gaia_DR3.source_id == EROSITASupersetv1Stars.gaia_dr3_source_id),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                EROSITASupersetv1Stars.xmatch_method == "HamStar",
                EROSITASupersetv1Stars.xmatch_version == "v1.1.1",
                EROSITASupersetv1Stars.xmatch_metric > 0.0,
                Gaia_DR3.phot_g_mean_mag > bright_bright_limit,
            )
        )

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.
        #

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


class MWM_EROSITA_Compact_Boss_Carton(BaseCarton):
    """MWM eROSITA Compact-Objects (unified for v1)

    Owner: eROSITA Compact objects WG Axel Schwope

    Shorthand name: mwm_erosita_compact_boss

    What is it?: eROSITA selected pointlike X-ray sources with
    likely Gaia counterpart (currently drawn from DR3),
    with high likelihood of being a galactic object
    according to proper motion and parallax,
    and not being an ordinary coronal emitter,
    also being (optically) reasonably faint,
    taking together being a candidate of a compact binary,
    i.e. accreting compact object.
    The original selection criteria that were applied
    to all X-ray/Gaia pairs are described on page Compact Binaries.
    These form the erosita_superset_v1_compactobjects.

    Simplified Description of selection criteria:
    "Select all targets from erosita_superset_compactobjects
    that have minimal distance to eROSITA eRASS:3 position
    and G > 13mag. All selected targets will receive BOSS spectroscopy"

    Wiki page: link or N/A

    Additional source catalogs needed: None

    Additional cross-matching needed: None

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    minimum exposure is 1x1 for bright 16<g < 17 mag),
    1x2 for medium (17<g<18)
    Fainter targets have 1x2 in the normal version of the cvarton,
    but have 2x2 in the _deep version of this carton.

    Pseudo SQL (optional):

    select * from erosita_superset_v1_compactobjects
    where xmatch_method = 'NWAY_CV' AND
    xmatch_version = '20Oct2022'

    reveals 11113 targets for 11087 X-ray sources

    Implementation:
    SELECT c.catalogid, c.ra, c.dec, gaia.phot_g_mean_mag,
    gaia.phot_bp_mean_mag, gaia.phot_rp_mean_mag, gaia.parallax,
    gaia.pmra, gaia.pmdec,
    estars.target_priority, estars.xmatch_metric,
    estars.ero_flux, estars.ero_ra, estars.ero_dec,
    estars.opt_ra, estars.opt_dec
    FROM catalog c
    INNER JOIN catalog_to_gaia_dr3_source c2g3 USING (catalogid)
    INNER JOIN gaia_dr3_source gaia ON gaia.source_id = c2g3.target_id
    INNER JOIN erosita_superset_v1_compactobjects estars
        ON estars.gaia_dr3_source_id = gaia.source_id
    WHERE (c2g3.version_id = 31) AND /* control version! */
    (c2g3.best is true) AND /* and enforce unique-ish crossmatch */
    estars.xmatch_method = 'NWAY_CV' AND
    estars.xmatch_version = '20Oct2022'
    limit 10;
    """

    name = "mwm_erosita_compact_boss"
    category = "science"
    instrument = "BOSS"
    program = "mwm_erosita"
    mapper = "MWM"
    can_offset = True
    faintest_cadence = "dark_flexible_2x2"
    flat_priority = 1810

    def build_query(self, version_id, query_region=None):
        cadence = peewee.Case(
            None,
            (
                (Gaia_DR3.phot_g_mean_mag < 17.0, "bright_1x1"),
                (Gaia_DR3.phot_g_mean_mag < 18.0, "dark_flexible_2x1"),
            ),
            self.faintest_cadence,
        ).cast("text")

        # priority = peewee.Case(
        #     None,
        #     (
        #         (Gaia_DR3.phot_g_mean_mag < 17.0, 2400),
        #     ),
        #     self.faintest_priority
        # ).cast('integer')
        priority = peewee.Value(self.flat_priority).cast("integer")

        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag.alias("gaia_g"),
                Gaia_DR3.phot_bp_mean_mag.alias("bp"),
                Gaia_DR3.phot_rp_mean_mag.alias("rp"),
                Gaia_DR3.parallax,
                Gaia_DR3.pmra,
                Gaia_DR3.pmdec,
                cadence.alias("cadence"),
                priority.alias("priority"),
                EROSITASupersetv1Compactobjects.target_priority,
                EROSITASupersetv1Compactobjects.ero_detuid,
                EROSITASupersetv1Compactobjects.xmatch_metric,
                EROSITASupersetv1Compactobjects.ero_flux,
                EROSITASupersetv1Compactobjects.ero_ra,
                EROSITASupersetv1Compactobjects.ero_dec,
                EROSITASupersetv1Compactobjects.opt_ra,
                EROSITASupersetv1Compactobjects.opt_dec,
            )
            .distinct(CatalogToGaia_DR3.catalogid)
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(
                EROSITASupersetv1Compactobjects,
                on=(Gaia_DR3.source_id == EROSITASupersetv1Compactobjects.gaia_dr3_source_id),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                EROSITASupersetv1Compactobjects.xmatch_method == "NWAY_CV",
                EROSITASupersetv1Compactobjects.xmatch_version == "20Oct2022",
            )
        )

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.
        #
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


class MWM_EROSITA_Compact_Boss_Shallow_Carton(MWM_EROSITA_Compact_Boss_Carton):
    name = "mwm_erosita_compact_boss_shallow"
    faintest_cadence = "dark_flexible_2x1"
    faintest_priority = 1812


class MWM_EROSITA_Compact_Boss_D3_Carton(MWM_EROSITA_Compact_Boss_Carton):
    name = "mwm_erosita_compact_boss_d3"
    faintest_cadence = "dark_flexible_3x1"
    faintest_priority = 1811
