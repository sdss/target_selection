#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2021-11-02
# @Filename: bhm_galaxies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToLegacy_Survey_DR10,
    Legacy_Survey_DR10,
)

from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import AB2nMgy


"""
Details: Start here
https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms
  see particularly: https://wiki.sdss.org/display/BHM/BHM+Cartons+of+Last+Resort

This module provides the following BHM cartons:
bhm_colr_galaxies_lsdr10
bhm_colr_galaxies_lsdr10_d3

"""


# used by cartons that need to compute Galactic latitude:
north_gal_pole_ra = 192.85948  # deg, J2000
north_gal_pole_dec = +27.12825  # deg, J2000


class BhmColrGalaxiesLsdr10Carton(BaseCarton):
    """
    A sample of bright galaxies selected from legacysurvey/dr10
    photometry+astrometry+morphology
    """

    name = "bhm_colr_galaxies_lsdr10"
    category = "science"
    mapper = "BHM"
    program = "bhm_filler"
    tile = False
    # we do not set cadence here since
    # cadence is set later below
    # cadence = 'dark_1x1'
    instrument = "BOSS"
    can_offset = False
    only_faintest_cadence = False

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()
        ls = Legacy_Survey_DR10.alias()

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get("priority", 10000)))
        value = peewee.Value(self.parameters.get("value", 0.0)).cast("float")
        inertial = peewee.Value(True)
        instrument = peewee.Value(self.instrument)
        # cadence = peewee.Value(self.parameters.get('cadence', self.cadence))

        dered_flux_z_min = AB2nMgy(self.parameters["dered_mag_z_max"])
        dered_fiberflux_z_min = AB2nMgy(self.parameters["dered_fibermag_z_max"])
        fiberflux_z_min = AB2nMgy(self.parameters["fibermag_z_max"])
        fiberflux_z_max = AB2nMgy(self.parameters["fibermag_z_min"])
        fiberflux_r_min = AB2nMgy(self.parameters["fibermag_r_max"])
        fiberflux_r_max = AB2nMgy(self.parameters["fibermag_r_min"])
        fiberflux_g_min = AB2nMgy(self.parameters["fibermag_g_max"])
        fiberflux_g_max = AB2nMgy(self.parameters["fibermag_g_min"])

        fiberflux_r_min_for_cadence1 = AB2nMgy(self.parameters["fibermag_r_for_cadence1"])
        fiberflux_r_min_for_cadence2 = AB2nMgy(self.parameters["fibermag_r_for_cadence2"])

        # compute transformed SDSS mags uniformly
        # transform the legacysurvey grz into sdss fiber2mag griz
        # https://wiki.sdss.org/display/BHM/BHM+magnitude+transformations+for+v0.5

        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_colr_galaxies_lsdr8/lsdr8_fibermag_to_sdss_fiber2mag_?_results.log   # noqa

        coeffs = {
            "g2_e": -0.611168,
            "g1_e": 1.637641,
            "g0_e": -0.684489,
            "i2_e": -1.821603,
            "i1_e": 2.410118,
            "i0_e": -0.895875,
            "r2_e": -0.337699,
            "r1_e": 1.038298,
            "r0_e": -0.371602,
            "z2_e": -1.609400,
            "z1_e": 2.678952,
            "z0_e": -0.821672,
        }

        nMgy_min = 1e-3  # equiv to AB=30
        # extended - start from ls10 fiberfluxes
        # TODO investigate using real lsdr10 i-band where available
        g0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g))
        r0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r))
        i0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_i))
        z0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_z))
        g_r_e = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.fiberflux_g)
            / peewee.fn.greatest(nMgy_min, ls.fiberflux_r)
        )
        r_z_e = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.fiberflux_r)
            / peewee.fn.greatest(nMgy_min, ls.fiberflux_z)
        )

        g_e = g0_e + coeffs["g0_e"] + coeffs["g1_e"] * g_r_e + coeffs["g2_e"] * g_r_e * g_r_e
        r_e = r0_e + coeffs["r0_e"] + coeffs["r1_e"] * g_r_e + coeffs["r2_e"] * g_r_e * g_r_e
        i_e = r0_e + coeffs["i0_e"] + coeffs["i1_e"] * r_z_e + coeffs["i2_e"] * r_z_e * r_z_e
        z_e = z0_e + coeffs["z0_e"] + coeffs["z1_e"] * r_z_e + coeffs["z2_e"] * r_z_e * r_z_e

        valid_colour_min_g_r = 0.0
        valid_colour_max_g_r = 2.0
        valid_colour_min_r_z = 0.2
        valid_colour_max_r_z = 1.2

        # select_flux_ratio_min_g_r = 10**(-0.4 * self.parameters['select_max_g_r'])
        # select_flux_ratio_max_g_r = 10**(-0.4 * self.parameters['select_min_g_r'])
        # select_flux_ratio_min_r_z = 10**(-0.4 * self.parameters['select_max_r_z'])
        # select_flux_ratio_max_r_z = 10**(-0.4 * self.parameters['select_min_r_z'])

        # magnitude validity checks
        valid = (
            g0_e.between(0.1, 29.9)
            & r0_e.between(0.1, 29.9)
            & z0_e.between(0.1, 29.9)
            & g_r_e.between(valid_colour_min_g_r, valid_colour_max_g_r)
            & r_z_e.between(valid_colour_min_r_z, valid_colour_max_r_z)
        )

        opt_prov = peewee.Case(None, ((valid, "sdss_fiber2mag_from_lsdr10"),), "undefined")
        magnitude_g = peewee.Case(None, ((valid, g_e),), "NaN")
        magnitude_r = peewee.Case(None, ((valid, r_e),), "NaN")
        magnitude_i = peewee.Case(None, ((valid, i_e),), "NaN")
        magnitude_z = peewee.Case(None, ((valid, z_e),), "NaN")
        magnitude_gaia_g = peewee.Case(
            None, ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),), "NaN"
        )
        magnitude_gaia_bp = peewee.Case(
            None, ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),), "NaN"
        )
        magnitude_gaia_rp = peewee.Case(
            None, ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),), "NaN"
        )

        opt_prov = peewee.Case(None, ((valid, "sdss_fiber2mag_from_lsdr10"),), "undefined")

        cadence1 = self.parameters["cadence1"]
        cadence2 = self.parameters["cadence2"]
        cadence3 = self.parameters["cadence3"]
        cadence = peewee.Case(
            None,
            (
                (ls.fiberflux_r > fiberflux_r_min_for_cadence1, cadence1),
                (ls.fiberflux_r > fiberflux_r_min_for_cadence2, cadence2),
            ),
            cadence3,
        )

        # compute the abs(Galactic latitude):
        gal_lat = peewee.fn.abs(
            90.0 - peewee.fn.q3c_dist(north_gal_pole_ra, north_gal_pole_dec, c.ra, c.dec)
        )

        # https://www.legacysurvey.org/dr10/bitmasks/
        maskbits_mask = 2**1 + 2**13
        # maskbits_mask = 2**0 + 2**1 + 2**4 + 2**7 + 2**8 + 2**10 + 2**11 + 2**13
        # fitbits_mask = 2**5 + 2**6 + 2**7 + 2**8 + 2**12
        query = (
            c.select(
                c.catalogid.alias("catalogid"),
                ls.ls_id.alias("ls_id"),  # extra
                c.ra.alias("ra"),  # extra
                c.dec.alias("dec"),  # extra
                priority.alias("priority"),
                value.cast("real").alias("value"),
                cadence.alias("cadence"),
                instrument.alias("instrument"),
                opt_prov.alias("optical_prov"),
                magnitude_g.cast("real").alias("g"),
                magnitude_r.cast("real").alias("r"),
                magnitude_i.cast("real").alias("i"),
                magnitude_z.cast("real").alias("z"),
                magnitude_gaia_g.alias("gaia_g"),
                magnitude_gaia_bp.alias("bp"),
                magnitude_gaia_rp.alias("rp"),
                inertial.alias("inertial"),
                g0_e.cast("real").alias("ls10_fibermag_g"),  # extra
                r0_e.cast("real").alias("ls10_fibermag_r"),  # extra
                i0_e.cast("real").alias("ls10_fibermag_i"),  # extra
                z0_e.cast("real").alias("ls10_fibermag_z"),  # extra
                g_r_e.cast("real").alias("ls10_g_r_e"),  # extra
                r_z_e.cast("real").alias("ls10_r_z_e"),  # extra
                ls.flux_g.alias("ls10_flux_g"),  # extra
                ls.flux_r.alias("ls10_flux_r"),  # extra
                ls.flux_i.alias("ls10_flux_i"),  # extra
                ls.flux_z.alias("ls10_flux_z"),  # extra
                ls.ebv.alias("ls10_ebv"),  # extra
                ls.maskbits.alias("ls10_maskbits"),  # extra
                ls.fitbits.alias("ls10_fitbits"),  # extra
                ls.mw_transmission_z.alias("ls10_mw_transmission_z"),  # extra
                gal_lat.alias("abs_gal_lat"),  # extra
            )
            .join(c2ls)
            .join(ls)
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                ls.type != "PSF",
                ls.parallax <= 0.0,
                ls.flux_z > dered_flux_z_min * ls.mw_transmission_z,
                ls.fiberflux_g > fiberflux_g_min,
                ls.fiberflux_r > fiberflux_r_min,
                ls.fiberflux_z > fiberflux_z_min,
                ls.fiberflux_z > dered_fiberflux_z_min * ls.mw_transmission_z,
                ls.fiberflux_g < fiberflux_g_max,
                ls.fiberflux_r < fiberflux_r_max,
                ls.fiberflux_z < fiberflux_z_max,
                # safety check using Gaia mags to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters["gaia_g_mag_limit"])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters["gaia_rp_mag_limit"])),
                ls.shape_r >= self.parameters["shape_r_min"],
                gal_lat > self.parameters["min_gal_lat"],
                ls.ebv < self.parameters["max_ebv"],
                (ls.maskbits.bin_and(maskbits_mask) == 0),  # avoid bad ls data
                # (ls.fitbits.bin_and(fitbits_mask) == 0),  # avoid bad ls fits
            )
            .distinct(c.catalogid)
        )

        if self.only_faintest_cadence:
            query = query.where(cadence == cadence3)

        # query_region[0] is ra of center of the region, degrees
        # query_region[1] is dec of center of the region, degrees
        # query_region[2] is radius of the region, degrees
        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    c.ra, c.dec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


class BhmColrGalaxiesLsdr10D3Carton(BhmColrGalaxiesLsdr10Carton):
    name = "bhm_colr_galaxies_lsdr10_d3"
    only_faintest_cadence = True
