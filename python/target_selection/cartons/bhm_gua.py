#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_gua.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR2,
    Gaia_DR2,
    CatalogFromSDSS_DR19p_Speclite,
    SDSS_DR19p_Speclite,
    Gaia_unWISE_AGN,
)
# DEBUG STUFF TO USE TEMP TABLE
# CatalogToSDSS_DR19p_Speclite._meta.table_name = 'temp_catalog_to_sdss_dr19p_speclite'
# CatalogToSDSS_DR19p_Speclite._meta._schema = 'sandbox'

from target_selection.cartons.base import BaseCarton

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms
#
# This module provides the following BHM cartons:
# bhm_gua_dark
# bhm_gua_dark_d3
# bhm_gua_bright
#
# Updated by TD on 27/01/2021 to add v0.5 improvements
# switch to a q3c join to sdss_dr16_specobj

"""
  pseudo-SQL

    both Cartons:
        SELECT catalogid,{derived fields}
        FROM catalog AS c
        JOIN catalog_to_tic AS c2tic ON ...
        JOIN tic_v8 AS tic ON ...
        JOIN gaia_dr2_source AS g ON ...
        JOIN gaia_unwise_agn AS t ON g.source_id = t.gaia_sourceid
        LEFT OUTER JOIN sdss_dr16_specobj AS s
          ON ( q3c_join(s.ra,s.dec,c.ra,c.dec,{match_radius_spectro})
               AND s.zwarning = 0
               AND s.snmedianl > 2.0
               AND s.zerr < 0.01
               AND s.scienceprimary > 0
             )
        WHERE
                c.version_id = {version_id}
            AND c2tic.version_id = {version_id}
            AND c2tic.best = True
            AND t.prob_rf > 0.8
            AND s.specobjid = NULL

    # bhm_gua_dark carton:
            AND ( t.g > 16.5 AND
                  t.rp > 16.5 AND
                  (t.g < 21.2 OR t.rp < 21.0 )
                )

    # bhm_gua_bright carton:
            AND ( t.g > 13.0 AND
                  t.rp > 13.5 AND
                  (t.g < 18.5 OR t.rp < 18.5)
                )
"""


class BhmGuaBaseCarton(BaseCarton):
    """
    Parent class that provides the basic selections for both Gaia UnWISE AGN cartons
    To be sub-classed, not to be called directly.

    To get from Catalog to GUA we join tables via :
    Catalog -> CatalogToGaia_DR2 -> Gaia_DR2 -> Gaia_unWISE_AGN

    Cadence+priority per target depends on brightness
    Keep the two bhm_gua_* cartons separate to enable overlap in brightness ranges
    """

    name = "bhm_gua_base"
    category = "science"
    mapper = "BHM"
    program = "bhm_filler"
    tile = False
    priority = None
    cadence = None
    instrument = "BOSS"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2g2 = CatalogToGaia_DR2.alias()
        g2 = Gaia_DR2.alias()
        t = Gaia_unWISE_AGN.alias()

        spec_sn_thresh = self.parameters["spec_sn_thresh"]
        spec_z_err_thresh = self.parameters["spec_z_err_thresh"]

        # #########################################################################
        # prepare the spectroscopy catalogue

        # SDSS DR19p
        # downslect only 'good' spectra
        c2s19 = CatalogFromSDSS_DR19p_Speclite.alias()
        ss19 = SDSS_DR19p_Speclite.alias()
        s19 = (
            ss19.select(
                ss19.pk.alias("s19_pk"),
            )
            .where(
                ss19.sn_median_all >= spec_sn_thresh,
                ss19.zwarning == 0,
                ss19.z_err <= spec_z_err_thresh,
                ss19.z_err > 0.0,
                ss19.specprimary > 0,
            )
            .alias("s19")
        )

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get("priority", 10000)))
        value = peewee.Value(self.parameters.get("value", 1.0)).cast("float")
        inertial = peewee.Value(True)
        cadence = peewee.Value(self.parameters["cadence"])
        instrument = peewee.Value(self.instrument)

        # compute transformed SDSS mags for pointlike and extended sources separately
        # transform the Gaia dr2 G,BP,RP into sdss psfmag griz

        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_gua/gdr2_mag_to_sdss_psfmag_?_results.log  # noqa
        coeffs = {
            "g3_p": 0.184158,
            "g2_p": -0.457316,
            "g1_p": 0.553505,
            "g0_p": -0.029152,
            "i3_p": 0.709818,
            "i2_p": -2.207549,
            "i1_p": 1.520957,
            "i0_p": -0.417666,
            "r3_p": 0.241611,
            "r2_p": -0.803702,
            "r1_p": 0.599944,
            "r0_p": -0.119959,
            "z3_p": 0.893988,
            "z2_p": -2.759177,
            "z1_p": 1.651668,
            "z0_p": -0.440676,
        }

        #  bp_rp = t.bp - t.rp
        #  g = (t.g + coeffs['g0_p'] + coeffs['g1_p'] * bp_rp + coeffs['g2_p'] * bp_rp * bp_rp +
        #       coeffs['g3_p'] * bp_rp * bp_rp * bp_rp)
        #  r = (t.g + coeffs['r0_p'] + coeffs['r1_p'] * bp_rp + coeffs['r2_p'] * bp_rp * bp_rp +
        #       coeffs['r3_p'] * bp_rp * bp_rp * bp_rp)
        #  i = (t.g + coeffs['i0_p'] + coeffs['i1_p'] * bp_rp + coeffs['i2_p'] * bp_rp * bp_rp +
        #       coeffs['i3_p'] * bp_rp * bp_rp * bp_rp)
        #  z = (t.g + coeffs['z0_p'] + coeffs['z1_p'] * bp_rp + coeffs['z2_p'] * bp_rp * bp_rp +
        #       coeffs['z3_p'] * bp_rp * bp_rp * bp_rp)

        g = (
            g2.phot_g_mean_mag
            + coeffs["g0_p"]
            + coeffs["g1_p"] * g2.bp_rp
            + coeffs["g2_p"] * g2.bp_rp * g2.bp_rp
            + coeffs["g3_p"] * g2.bp_rp * g2.bp_rp * g2.bp_rp
        )
        r = (
            g2.phot_g_mean_mag
            + coeffs["r0_p"]
            + coeffs["r1_p"] * g2.bp_rp
            + coeffs["r2_p"] * g2.bp_rp * g2.bp_rp
            + coeffs["r3_p"] * g2.bp_rp * g2.bp_rp * g2.bp_rp
        )
        i = (
            g2.phot_g_mean_mag
            + coeffs["i0_p"]
            + coeffs["i1_p"] * g2.bp_rp
            + coeffs["i2_p"] * g2.bp_rp * g2.bp_rp
            + coeffs["i3_p"] * g2.bp_rp * g2.bp_rp * g2.bp_rp
        )
        z = (
            g2.phot_g_mean_mag
            + coeffs["z0_p"]
            + coeffs["z1_p"] * g2.bp_rp
            + coeffs["z2_p"] * g2.bp_rp * g2.bp_rp
            + coeffs["z3_p"] * g2.bp_rp * g2.bp_rp * g2.bp_rp
        )

        # validity checks - set limits semi-manually
        bp_rp_min = 0.0
        bp_rp_max = 1.8
        valid = (
            g2.phot_g_mean_mag.between(0.1, 29.9)
            & g2.phot_bp_mean_mag.between(0.1, 29.9)
            & g2.phot_rp_mean_mag.between(0.1, 29.9)
            & g2.bp_rp.between(bp_rp_min, bp_rp_max)
        )

        opt_prov = peewee.Case(None, ((valid, "sdss_psfmag_from_gaiadr2"),), "undefined")
        magnitude_g = peewee.Case(None, ((valid, g),), "NaN")
        magnitude_r = peewee.Case(None, ((valid, r),), "NaN")
        magnitude_i = peewee.Case(None, ((valid, i),), "NaN")
        magnitude_z = peewee.Case(None, ((valid, z),), "NaN")

        # Create temporary tables for the base query and the Q3C cross-match
        # tables.

        query = (
            c.select(
                c.catalogid,
                c.ra,  # extra
                c.dec,  # extra
                t.gaia_sourceid.alias("gaia_dr2_source_id"),  # extra
                t.unwise_objid,  # extra
                priority.alias("priority"),
                value.alias("value"),
                inertial.alias("inertial"),
                cadence.alias("cadence"),
                instrument.alias("instrument"),
                opt_prov.alias("optical_prov"),
                magnitude_g.alias("g"),
                magnitude_r.alias("r"),
                magnitude_i.alias("i"),
                magnitude_z.alias("z"),
                g2.phot_g_mean_mag.alias("gaia_g"),
                g2.phot_bp_mean_mag.alias("bp"),
                g2.phot_rp_mean_mag.alias("rp"),
                t.w1.alias("gua_w1"),  # extra
                t.w2.alias("gua_w2"),  # extra
                t.prob_rf.alias("gua_prob_rf"),  # extra
                t.phot_z.alias("gua_phot_z"),  # extra
                s19.c.s19_pk.alias("s19_pk"),  # extra
                # rely on the centralised magnitude routines for 'real' griz, bp,rp,gaia_g
            )
            .join(c2g2)
            .join(g2, on=(c2g2.target_id == g2.source_id))
            .join(t, on=(g2.source_id == t.gaia_sourceid))
            # joining to spectroscopy
            .switch(c)
            .join(c2s19, JOIN.LEFT_OUTER)
            .join(s19, JOIN.LEFT_OUTER, on=(s19.c.s19_pk == c2s19.target_id))
            .where(
                c.version_id == version_id,
                c2g2.version_id == version_id,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
                c2g2.best >> True,
                # fn.coalesce(c2s19.best, True) >> True,
            )
            .where(
                (t.prob_rf >= self.parameters["prob_rf_min"]),
                (g2.phot_g_mean_mag >= self.parameters["mag_g_min"]),
                (g2.phot_rp_mean_mag >= self.parameters["mag_rp_min"]),
                (
                    (g2.phot_g_mean_mag < self.parameters["mag_g_max"])
                    | (g2.phot_rp_mean_mag < self.parameters["mag_rp_max"])
                ),
            )
            # then reject any GUA targets with existing good DR19p spectroscopy
            .where(s19.c.s19_pk.is_null(True))
            # avoid duplicates - trust the gaia ids in the GUA parent sample
            .distinct([t.gaia_sourceid])
        )

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    c.ra, c.dec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


class BhmGuaDarkCarton(BhmGuaBaseCarton):
    """
    -------  bhm_gua_dark   ------
    SQL as above plus:

        AND ( gua.g > 16.x AND gua.rp > 16.x)
    """

    name = "bhm_gua_dark"


class BhmGuaDarkD3Carton(BhmGuaBaseCarton):
    name = "bhm_gua_dark_d3"


class BhmGuaBrightCarton(BhmGuaBaseCarton):
    """
    ------  bhm_gua_bright   ------
    SQL as above plus:

        AND ( gua.g < 18.x OR gua.rp < 18.x)

    """

    name = "bhm_gua_bright"
