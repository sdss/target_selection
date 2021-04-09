#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-02
# @Filename: mwm_snc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2, TIC_v8)

from target_selection.cartons import BaseCarton


class MWM_SNC_100pc_Carton(BaseCarton):
    """MWM Solar Neighbourhood Census 100pc.

    Definition:

        - Divide the sky in crowded/not crowded regions via simple geometric
          cuts in (l,b).
        - Select all stars that fall within their one-sigma uncertainty within
          100pc in the not crowded regions.
        - Select all stars that fall within their one-sigma uncertainty
          within 100pc and have an astrometric_excess_noise < 2 in the
          crowded regions.

    Pseudo SQL:

        parallax-parallax_error>10 &&
        ( (astrometric_excess_noise < 2 &&
           ( (l <= 180 && b < -0.139 * l + 25 && b > 0.139 * l - 25) ||
             (l > 180 && b > -0.139 * l + 25 && b < 0.139 * l - 25) ||
             (sqrt(pow(l - 303.2, 2) + 2 * pow(b + 44.4, 2)) < 5) ||
             (sqrt(pow(l - 280.3, 2) + 2 * pow(b + 33.0, 2)) < 8) )
        ) ||
        !( (l <= 180 && b < -0.139 * l + 25 && b > 0.139 * l - 25) ||
           (l > 180 && b > -0.139 * l + 25 && b < 0.139 * l - 25) ||
           (sqrt(pow(l - 303.2, 2) + 2 * pow(b + 44.4, 2)) < 5) ||
           (sqrt(pow(l - 280.3, 2) + 2 * pow(b + 33.0, 2)) < 8) )
        )

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):

        l = Gaia_DR2.l  # noqa
        b = Gaia_DR2.b

        # Dense regions (Galactic plane, SMC, LMC).
        gal_cut = (((l <= 180) & (b < (-0.139 * l + 25)) & (b > (0.139 * l - 25))) |  # noqa
                   ((l > 180) & (b > (-0.139 * l + 25)) & (b < (0.139 * l - 25))) |
                   (fn.sqrt(fn.pow(l - 303.2, 2) + 2 * fn.pow(b + 44.4, 2)) < 5) |
                   (fn.sqrt(fn.pow(l - 280.3, 2) + 2 * fn.pow(b + 33.0, 2)) < 8))

        query = (Gaia_DR2
                 .select(CatalogToTIC_v8.catalogid,
                         TIC_v8.gaia_int.alias('gaia_source_id'),
                         Gaia_DR2.phot_g_mean_mag,
                         l, b)
                 .join(TIC_v8)
                 .join(CatalogToTIC_v8)
                 .where((Gaia_DR2.parallax - Gaia_DR2.parallax_error) > 10,
                        ((Gaia_DR2.astrometric_excess_noise < 2) & gal_cut) | ~(gal_cut))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(fn.q3c_radial_query(Catalog.ra,
                                                Catalog.dec,
                                                query_region[0],
                                                query_region[1],
                                                query_region[2])))

        return query


class MWM_SNC_100pc_APOGEE_Carton(MWM_SNC_100pc_Carton):
    """MWM SNC 100pc targets to be observed with APOGEE."""

    name = 'mwm_snc_100pc_apogee'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 1805

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region=query_region)

        return query.where(TIC_v8.hmag < 11)


class MWM_SNC_100pc_BOSS_Carton(MWM_SNC_100pc_Carton):
    """MWM SNC 100pc targets to be observed with BOSS."""

    name = 'mwm_snc_100pc_boss'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'BOSS'
    priority = 1800

    def post_process(self, model, **kwargs):

        # G > 16 => cadence = dark_2x1
        model.update(cadence='dark_2x1').where(model.phot_g_mean_mag > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence='bright_2x1').where(model.phot_g_mean_mag < 16).execute()

        return model


class MWM_SNC_250pc_Carton(BaseCarton):
    """MWM Solar Neighbourhood Census 250pc.

    Definition:

        - Divide the sky in crowded/not crowded regions via simple geometric
          cuts in (l,b)
        - Apply a general cut of absolute G-band magnitude < 6mag,
          corresponding to ~G-type stars and earlier
        - Select all stars that fall within their one-sigma uncertainty
          within 250pc in the not crowded regions
        - Select all stars that fall within their one-sigma uncertainty within
          250pc and have an astrometric_excess_noise < 2 in the crowded
          regions.

    Pseudo SQL:

        (phot_g_mean_mag + 5 * log10(parallax / 1000) + 5 < 6) &&
        parallax - parallax_error > 4 &&
        ( (astrometric_excess_noise < 2 &&
           ( (l <= 180 && b < -0.139 * l + 25 && b > 0.139 * l - 25) ||
             (l > 180 && b > -0.139 * l + 25 && b < 0.139 * l - 25) ||
             (sqrt(pow(l - 303.2, 2) + 2 * pow(b + 44.4, 2)) < 5) ||
             (sqrt(pow(l - 280.3, 2) + 2 * pow(b + 33.0, 2)) < 8))) ||
          !( (l <= 180 && b < -0.139 * l + 25 && b > 0.139 * l - 25) ||
             (l > 180 && b > -0.139 * l + 25 && b < 0.139 * l - 25) ||
             (sqrt(pow(l - 303.2, 2) + 2 * pow(b + 44.4, 2)) < 5) ||
             (sqrt(pow(l - 280.3, 2) + 2 * pow(b + 33.0, 2)) < 8)
            )
        )

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):

        l = Gaia_DR2.l  # noqa
        b = Gaia_DR2.b

        colour_cut = (Gaia_DR2.phot_g_mean_mag +
                      5 * fn.log(Gaia_DR2.parallax / 1000) + 5) < 6

        gal_cut = (((l <= 180) & (b < (-0.139 * l + 25)) & (b > (0.139 * l - 25))) |  # noqa
                   ((l > 180) & (b > (-0.139 * l + 25)) & (b < (0.139 * l - 25))) |
                   (fn.sqrt(fn.pow(l - 303.2, 2) + 2 * fn.pow(b + 44.4, 2)) < 5) |
                   (fn.sqrt(fn.pow(l - 280.3, 2) + 2 * fn.pow(b + 33.0, 2)) < 8))

        query = (Gaia_DR2
                 .select(CatalogToTIC_v8.catalogid,
                         TIC_v8.gaia_int.alias('gaia_source_id'),
                         Gaia_DR2.phot_g_mean_mag,
                         l, b)
                 .join(TIC_v8)
                 .join(CatalogToTIC_v8)
                 .where(colour_cut,
                        (Gaia_DR2.parallax - Gaia_DR2.parallax_error) > 4,
                        ((Gaia_DR2.astrometric_excess_noise < 2) & gal_cut) | ~(gal_cut))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(fn.q3c_radial_query(Catalog.ra,
                                                Catalog.dec,
                                                query_region[0],
                                                query_region[1],
                                                query_region[2])))

        return query


class MWM_SNC_250pc_APOGEE_Carton(MWM_SNC_250pc_Carton):
    """MWM SNC 250pc targets to be observed with APOGEE."""

    name = 'mwm_snc_250pc_apogee'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 1815

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region=query_region)

        return query.where(TIC_v8.hmag < 11)


class MWM_SNC_250pc_BOSS_Carton(MWM_SNC_250pc_Carton):
    """MWM SNC 250pc targets to be observed with BOSS."""

    name = 'mwm_snc_250pc_boss'
    category = 'science'
    program = 'mwm_snc'
    mapper = 'MWM'
    instrument = 'BOSS'
    priority = 1810

    def post_process(self, model, **kwargs):

        # G > 16 => cadence = dark_2x1
        model.update(cadence='dark_2x1').where(model.phot_g_mean_mag > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence='bright_2x1').where(model.phot_g_mean_mag < 16).execute()

        return model
