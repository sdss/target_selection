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

from . import BaseCarton


__all__ = ['MWM_SNC_100pc', 'MWM_SNC_250pc']


class MWM_SNC_100pc(BaseCarton):
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

    """

    name = 'mwm_100pc'
    category = 'science'
    program = 'SNC'
    mapper = 'MWM'

    def build_query(self, version_id, query_region=None):

        l = TIC_v8.gallong  # noqa
        b = TIC_v8.gallat

        gal_cut = (((l <= 180) & (b < (-0.139 * l + 25)) & (b > (0.139 * l - 25))) |  # noqa
                   ((l > 180) & (b > (-0.139 * l + 25)) & (b < (0.139 * l - 25))) |
                   (fn.sqrt(fn.pow(l - 303.2, 2) + 2 * fn.pow(b + 44.4, 2)) < 5) |
                   (fn.sqrt(fn.pow(l - 280.3, 2) + 2 * fn.pow(b + 33.0, 2)) < 8))

        query = (Gaia_DR2
                 .select(Catalog.catalogid,
                         TIC_v8.gaia_int.alias('gaia_source_id'),
                         l, b)
                 .join(TIC_v8)
                 .join(CatalogToTIC_v8)
                 .join(Catalog)
                 .where((TIC_v8.plx - TIC_v8.e_plx) > 10,
                        ((Gaia_DR2.astrometric_excess_noise < 2) & gal_cut) |
                        ~(gal_cut))
                 .where(CatalogToTIC_v8.version_id == version_id))

        if query_region:
            query = query.where(fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                    query_region[0],
                                                    query_region[1],
                                                    query_region[2]))

        return query


class MWM_SNC_250pc(BaseCarton):
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

    """

    name = 'mwm_250pc'
    category = 'science'
    program = 'SNC'
    mapper = 'MWM'

    def build_query(self, version_id, query_region=None):

        l = TIC_v8.gallong  # noqa
        b = TIC_v8.gallat

        colour_cut = (TIC_v8.gaiamag + 5 * fn.log(TIC_v8.plx / 1000) + 5) < 6
        gal_cut = (((l <= 180) & (b < (-0.139 * l + 25)) & (b > (0.139 * l - 25))) |  # noqa
                   ((l > 180) & (b > (-0.139 * l + 25)) & (b < (0.139 * l - 25))) |
                   (fn.sqrt(fn.pow(l - 303.2, 2) + 2 * fn.pow(b + 44.4, 2)) < 5) |
                   (fn.sqrt(fn.pow(l - 280.3, 2) + 2 * fn.pow(b + 33.0, 2)) < 8))

        query = (Gaia_DR2
                 .select(Catalog.catalogid,
                         TIC_v8.gaia_int.alias('gaia_source_id'),
                         l, b)
                 .join(TIC_v8)
                 .join(CatalogToTIC_v8)
                 .join(Catalog)
                 .where(colour_cut,
                        (TIC_v8.plx - TIC_v8.e_plx) > 4,
                        ((Gaia_DR2.astrometric_excess_noise < 2) & gal_cut) |
                        ~(gal_cut))
                 .where(Catalog.version_id == version_id,
                        CatalogToTIC_v8.version_id == version_id))

        if query_region:
            query = query.where(fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                    query_region[0],
                                                    query_region[1],
                                                    query_region[2]))

        return query
