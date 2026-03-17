#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_spiders_agn.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

# ### flake8: noqa
# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn

from target_selection.cartons.base import BaseCarton

# general catalogdb imports
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    EROSITASupersetv1AGN,
)

# imports of existing spectro catalogue
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogFromSDSS_DR19p_Speclite,
    SDSS_DR19p_Speclite,
)

# additional imports required by bhm_spiders_agn_lsdr10
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToLegacy_Survey_DR10,
    Legacy_Survey_DR10,
)

# additional imports required by bhm_spiders_agn_gaiadr3
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToGaia_DR3,
    Gaia_DR3,
)

# # additional imports required by bhm_spiders_agn_ps1dr2
# from sdssdb.peewee.sdss5db.catalogdb import (
#     Panstarrs1,
#     CatalogToPanstarrs1,
# )

from target_selection.mag_flux import AB2nMgy
# from target_selection.mag_flux import AB2Jy

# DEBUG STUFF TO USE TEMP TABLE
# CatalogToSDSS_DR19p_Speclite._meta.table_name = 'temp_catalog_to_sdss_dr19p_speclite'
# CatalogToSDSS_DR19p_Speclite._meta._schema = 'sandbox'

# used by cartons that need to compute Galactic latitude:
north_gal_pole_ra = 192.85948  # deg, J2000
north_gal_pole_dec = +27.12825  # deg, J2000

# maskbits to determine if X-ray detection is potentially problematic
# https://wiki.mpe.mpg.de/eRosita/EroFollowup/SDSSV_data_model#Notes_on_ero_flags_column_carrying_eRASS_quality_flagging_information # noqa
ero_flags_mask = (
    0
    + 2**0  # FLAG_SP_SNR - eRASS source is located in/near a supernova remnant # noqa
    + 2**1  # FLAG_SP_BPS - eRASS source is located in/near a bright (X-ray?) point source # noqa
    + 2**2  # FLAG_SP_SCL - eRASS source is located in/near a star cluster # noqa
    + 2**3  # FLAG_SP_LGA - eRASS source is located in/near a local galaxy # noqa
    + 2**4  # FLAG_SP_GC - eRASS source is located in/near a globular cluster # noqa
    + 2
    ** 5  # FLAG_SP_GC_CONS - eRASS source is located in/near a globular cluster (conservative criteria) # noqa
    # + 2**6   # FLAG_NO_RADEC_ERR - eRASS source is missing an estimate of uncertainty on its sky position # noqa
    # + 2**7   # FLAG_NO_EXT_ERR - eRASS source is missing an estimate of uncertainty on its X-ray extent # noqa
    # + 2**8   # FLAG_NO_CTS_ERR - eRASS source is missing an estimate of uncertainty on the number of detected X-ray counts # noqa
    # + 2**9   # FLAG_HARD_SRC   - eRASS source is significantly detected in the 2.3-5keV band (in either eRASS:1 3B or eRASS:3 3B, with DET_LIKE_3 > 8, 30arcsec radius ) # noqa
)

# ############################################
# ############################################
# ############################################
# ############################################
# This provides the following BHM SPIDERS AGN cartons in v1.0:
#  *  bhm_spiders_agn_lsdr10
#  *  bhm_spiders_agn_gaiadr3
#  *  bhm_spiders_agn_sep
#  *  bhm_spiders_agn_tda
#  *  bhm_spiders_agn_hard
#  *  bhm_spiders_agn_lsdr10_d3
#  *  bhm_spiders_agn_gaiadr3_d3
#  *  bhm_spiders_agn_sep_d3
#  *  bhm_spiders_agn_tda_d3
#  *  bhm_spiders_agn_hard_d3
# ############################################
# ############################################
# ############################################
# ############################################

# some reference points for AB->nMgy conversions
# 30.0 AB = 1e-3 nMgy
# 22.5 AB = 1.0 nMgy
# 22.0 AB = 1.58489 nMgy
# 21.5 AB = 2.51189 nMgy
# 21.0 AB = 3.98107 nMgy
# 20.0 AB = 10.0 nMgy
# 18.5 AB = 39.8107 nMgy
# 16.5 AB = 251.189 nMgy
# 14.0 AB = 2511.89 nMgy
# 13.5 AB = 3981.07 nMgy

# some reference points for AB->Jy conversions (for ps1_dr2 fluxes)
# 30.0 AB = 3.631e-9 Jy
# 22.5 AB = 3.631e-6 Jy
# 22.0 AB = 5.754e-6 Jy
# 21.5 AB = 9.120e-6 Jy
# 21.0 AB = 1.445e-5 Jy
# 20.5 AB = 2.291e-5  Jy
# 18.5 AB = 1.445e-4 Jy
# 16.5 AB = 9.120e-4 Jy
# 14.0 AB = 9.120e-3 Jy
# 13.5 AB = 1.445e-2 Jy

# Notes on catalogdb.panstarrs1.flags aka objInfoFlag from ObjectThin
# https://outerspace.stsci.edu/display/PANSTARRS/PS1+ObjectThin+table+fields
# https://outerspace.stsci.edu/display/PANSTARRS/PS1+Object+Flags
# select objects that have the GOOD_STACK flag set:
# Flag name   value      decimal   Notes
# GOOD_STACK  0x08000000 134217728 good-quality object in the stack (> 1 good stack measurement)

# Use these two flags to decide whether to use aper mags or not
# Flag name   value      decimal   Notes
# EXT	      0x00800000   8388608 extended in our data (eg, PS)
# EXT_ALT     0x01000000  16777216 extended in external data (eg, 2MASS)

"""
# noqa
# Notes on how many targets to expect:
 sdss5db=> SELECT ero_version,xmatch_method,xmatch_version,opt_cat,ero_flux_type,count(*)
           FROM erosita_superset_v1_agn GROUP BY ero_version,xmatch_method,xmatch_version,opt_cat,ero_flux_type;
           ero_version           |  xmatch_method  |      xmatch_version      | opt_cat | ero_flux_type |  count
---------------------------------+-----------------+--------------------------+---------+---------------+---------
 eRASS_s1_3B_221031_poscorr      | XPS/NWAY        | JBJWMS_24Nov22           | lsdr10  | 2.3-5keV      |    3433
 eRASS_s3_1B_220829_poscorr_v006 | XPS/NWAY_EROTDA | JWMS_06Oct22_erotda      | lsdr10  | 0.2-2.3keV    |    9703
 eRASS_s3_1B_221007_poscorr_v007 | XPS/NWAY        | JWMS_06Feb23_nomask      | lsdr10  | 0.2-2.3keV    | 1974450
 eRASS_s3_1B_221007_poscorr_v007 | XPS/NWAY        | JWMS_21Oct22             | lsdr10  | 0.2-2.3keV    | 1895479
 eRASS_s3_1B_221007_poscorr_v007 | XPS/NWAY        | JWMS_21Oct22             | lsdr9   | 0.2-2.3keV    |   47172
 eRASS_s3_1B_221007_poscorr_v007 | XPS/NWAY        | JWMS_21Oct22_cw2020_gedr | gedr3   | 0.2-2.3keV    | 1298743
 eRASS_s3_1B_221007_poscorr_v007 | XPS/NWAY        | JWMS_24Oct22             | gedr3   | 0.2-2.3keV    | 2465166
 eRASS_s3_1B_221007_poscorr_v007 | XPS/NWAY        | JWMS_24Oct22_nomask      | lsdr10  | 0.2-2.3keV    | 1937267
 eRASS_s5_V29C                   | XPS/NWAY        | JWTL_Oct22               | gedr3   | 0.2-2.3keV    |    2007
 eRASS_s5_V29C                   | XPS/NWAY        | JWTL_Oct22               | lsdr10  | 0.2-2.3keV    |    5207
(10 rows)
"""

# Notes on avoiding saturated legacysurvey sources
# https://www.legacysurvey.org/dr8/bitmasks/
# Bit   Name       Description
# 0     NPRIMARY   touches a pixel that is outside the BRICK_PRIMARY region of a brick
# 1     BRIGHT     touches a pixel within the locus of a radius-magnitude relation for
#                  Tycho-2 stars or one for Gaia DR2 stars to G < 13
# 2     SATUR_G    touches a pixel that was saturated in at least one g-band image
# 3     SATUR_R    touches a pixel that was saturated in at least one r-band image
# 4     SATUR_Z    touches a pixel that was saturated in at least one z-band image
# 5     ALLMASK_G  touches a pixel that has any of the ALLMASK_G bits set
# 6     ALLMASK_R  touches a pixel that has any of the ALLMASK_R bits set
# 7     ALLMASK_Z  touches a pixel that has any of the ALLMASK_Z bits set
# 8     WISEM1     touches a pixel in a WISEMASK_W1 bright star mask
# 9     WISEM2     touches a pixel in a WISEMASK_W2 bright star mask
# 10    BAILOUT    touches a pixel in a blob where we "bailed out" of source fitting
# 11    MEDIUM     touches a pixel within the locus of a radius-magnitude relation
#                  for Gaia DR2 stars to G < 16
# 12    GALAXY     touches a pixel in an SGA large galaxy
# 13    CLUSTER    touches a pixel in a globular cluster
#
# so, mask to avoid saturated targets is 2**2 + 2**3 + 2**4 = 4+8+16 = 28

#
# END PREAMBLE
# ##################################################################################


class BhmSpidersAgnLsdr10Carton(BaseCarton):
    name = "bhm_spiders_agn_lsdr10"
    category = "science"
    mapper = "BHM"
    program = "bhm_spiders"
    tile = False
    instrument = "BOSS"
    can_offset = True
    only_faintest_cadence = False

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        x = EROSITASupersetv1AGN.alias()
        ls = Legacy_Survey_DR10.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()

        instrument = peewee.Value(self.instrument)

        fiberflux_r_max = AB2nMgy(self.parameters["fibermag_r_min"])
        fiberflux_r_min = AB2nMgy(self.parameters["fibermag_r_max"])
        fiberflux_i_max = AB2nMgy(self.parameters["fibermag_i_min"])
        fiberflux_i_min = AB2nMgy(self.parameters["fibermag_i_max"])
        fiberflux_z_max = AB2nMgy(self.parameters["fibermag_z_min"])
        fiberflux_z_min = AB2nMgy(self.parameters["fibermag_z_max"])
        fiberflux_r_max_for_core = AB2nMgy(self.parameters["fibermag_r_min_for_core"])
        fiberflux_r_min_for_core = AB2nMgy(self.parameters["fibermag_r_max_for_core"])
        fiberflux_i_max_for_core = AB2nMgy(self.parameters["fibermag_i_min_for_core"])
        fiberflux_i_min_for_core = AB2nMgy(self.parameters["fibermag_i_max_for_core"])
        fiberflux_z_max_for_core = AB2nMgy(self.parameters["fibermag_z_min_for_core"])
        fiberflux_z_min_for_core = AB2nMgy(self.parameters["fibermag_z_max_for_core"])

        fiberflux_r_min_for_cadence1 = AB2nMgy(self.parameters["fibermag_r_for_cadence1"])
        fiberflux_r_min_for_cadence2 = AB2nMgy(self.parameters["fibermag_r_for_cadence2"])
        fiberflux_i_min_for_cadence1 = AB2nMgy(self.parameters["fibermag_i_for_cadence1"])
        fiberflux_i_min_for_cadence2 = AB2nMgy(self.parameters["fibermag_i_for_cadence2"])
        gaia_g_max_for_cadence1 = self.parameters["gaia_g_max_for_cadence1"]
        gaia_rp_max_for_cadence1 = self.parameters["gaia_rp_max_for_cadence1"]

        # #########################################################################
        # prepare the spectroscopy catalogue
        spec_sn_thresh = self.parameters["spec_sn_thresh"]
        spec_z_err_thresh = self.parameters["spec_z_err_thresh"]

        # SDSS DR19p
        # first downslect only 'good' spectra
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
        #########################################################################

        # compute the abs(Galactic latitude):
        gal_lat = peewee.fn.abs(
            90.0 - peewee.fn.q3c_dist(north_gal_pole_ra, north_gal_pole_dec, c.ra, c.dec)
        )

        # logic is written this backwards way so that a failure to meet any core
        # criterion results in non-core status
        is_core = peewee.Case(
            None,
            (
                (gal_lat < self.parameters["min_gal_lat_for_core"], False),
                (c.dec < self.parameters["min_dec_for_core"], False),
                (x.ero_flux < self.parameters["min_ero_flux_for_core"], False),
                (x.ero_det_like < self.parameters["min_det_like_for_core"], False),
                (
                    ~(
                        (
                            ls.fiberflux_r.between(
                                fiberflux_r_min_for_core, fiberflux_r_max_for_core
                            )
                        )
                        | (
                            ls.fiberflux_i.between(
                                fiberflux_i_min_for_core, fiberflux_i_max_for_core
                            )
                        )
                        | (
                            ls.fiberflux_z.between(
                                fiberflux_z_min_for_core, fiberflux_z_max_for_core
                            )
                        )
                    ),
                    False,
                ),
                (ls.maskbits.bin_and(2**13) > 0, False),  # demote globular clusters and MCs
                (x.ero_flags.bin_and(ero_flags_mask) > 0, False),  # demote problematic X-ray data
            ),
            True,
        )

        # value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        value = peewee.Case(None, ((is_core, self.parameters.get("value", 1.0)),), 0.0).cast(
            "float"
        )

        # priority is determined from individual target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +dpriority_match_flags if target is a secondary cross-match (match_flag > 1)
        # add +dpriority_det_like if target has a low value of ero_det_like
        # add +dpriority_has_spec if target has existing good SDSS spectroscopy
        priority_1 = peewee.Case(None, ((x.xmatch_flags > 1, 1),), 0)
        priority_2 = peewee.Case(
            None, ((x.ero_det_like < self.parameters["det_like_for_priority"], 1),), 0
        )
        priority_3 = peewee.Case(None, ((s19.c.s19_pk.is_null(False), 1),), 0)
        priority_4 = peewee.Case(None, ((is_core, 0),), 1)
        priority_5 = peewee.Case(None, ((x.ero_flags.bin_and(2**9) == 0, 1),), 0)

        # priority = fn.max(
        priority = (
            self.parameters["priority_floor"]
            + priority_1 * self.parameters["dpriority_match_flags"]
            + priority_2 * self.parameters["dpriority_det_like"]
            + priority_3 * self.parameters["dpriority_has_spec"]
            + priority_4 * self.parameters["dpriority_non_core"]
            + priority_5 * self.parameters["dpriority_not_hard"]
        )

        # choose cadence based on fiber magnitude in r-band
        cadence1 = self.parameters["cadence1"]
        cadence2 = self.parameters["cadence2"]
        cadence3 = self.parameters["cadence3"]
        # cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                (
                    (
                        (ls.fiberflux_r > fiberflux_r_min_for_cadence1)
                        | (ls.fiberflux_i > fiberflux_i_min_for_cadence1)
                        | (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1))
                        | (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))
                    ),
                    cadence1,
                ),
                (
                    (ls.fiberflux_r > fiberflux_r_min_for_cadence2)
                    | (ls.fiberflux_i > fiberflux_i_min_for_cadence2),
                    cadence2,
                ),
                #  ((ls.fiberflux_r < fiberflux_r_min_for_cadence2) &
                #   (ls.fiberflux_i < fiberflux_i_min_for_cadence2),
                #    cadence3),
            ),
            cadence3,
        )
        #    cadence4)

        # compute transformed SDSS mags for pointlike and extended sources separately
        # transform the legacysurvey grz into sdss psfmag griz

        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_agn_lsdr8_*/lsdr8_*mag_to_sdss_*mag_?_results.log  # noqa
        coeffs = {
            "g2_e": -0.113816,
            "g1_e": 0.317176,
            "g0_e": 0.094145,
            "i2_e": -0.415858,
            "i1_e": 0.168922,
            "i0_e": -0.010771,
            "r2_e": 0.029398,
            "r1_e": -0.019938,
            "r0_e": 0.354042,
            "z2_e": -0.111262,
            "z1_e": 0.237656,
            "z0_e": 0.148923,
            "g2_p": 0.187193,
            "g1_p": -0.184362,
            "g0_p": 0.049492,
            "i2_p": -0.098979,
            "i1_p": -0.405518,
            "i0_p": 0.009688,
            "r2_p": -0.001935,
            "r1_p": 0.098201,
            "r0_p": 0.050321,
            "z2_p": -0.034163,
            "z1_p": 0.109878,
            "z0_p": -0.030167,
        }

        nMgy_min = 1e-3  # equiv to AB=30
        # pointlike - start from ls8 (psf)fluxes
        g0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g))
        r0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r))
        i0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_i))
        z0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z))
        g_r_p = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g) / peewee.fn.greatest(nMgy_min, ls.flux_r)
        )
        r_z_p = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r) / peewee.fn.greatest(nMgy_min, ls.flux_z)
        )

        # extended - start from ls8 fiberfluxes
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

        g_p = g0_p + coeffs["g0_p"] + coeffs["g1_p"] * g_r_p + coeffs["g2_p"] * g_r_p * g_r_p
        r_p = r0_p + coeffs["r0_p"] + coeffs["r1_p"] * g_r_p + coeffs["r2_p"] * g_r_p * g_r_p
        i_p = r0_p + coeffs["i0_p"] + coeffs["i1_p"] * r_z_p + coeffs["i2_p"] * r_z_p * r_z_p
        z_p = z0_p + coeffs["z0_p"] + coeffs["z1_p"] * r_z_p + coeffs["z2_p"] * r_z_p * r_z_p

        g_e = g0_e + coeffs["g0_e"] + coeffs["g1_e"] * g_r_e + coeffs["g2_e"] * g_r_e * g_r_e
        r_e = r0_e + coeffs["r0_e"] + coeffs["r1_e"] * g_r_e + coeffs["r2_e"] * g_r_e * g_r_e
        i_e = r0_e + coeffs["i0_e"] + coeffs["i1_e"] * r_z_e + coeffs["i2_e"] * r_z_e * r_z_e
        z_e = z0_e + coeffs["z0_e"] + coeffs["z1_e"] * r_z_e + coeffs["z2_e"] * r_z_e * r_z_e

        # validity checks - set limits semi-manually
        g_r_p_min = -0.25
        g_r_p_max = 1.75
        r_z_p_min = -0.5
        r_z_p_max = 2.5
        g_r_e_min = 0.0
        g_r_e_max = 1.75
        r_z_e_min = 0.2
        r_z_e_max = 1.6
        valid_p = (
            g0_p.between(0.1, 29.9)
            & r0_p.between(0.1, 29.9)
            & z0_p.between(0.1, 29.9)
            & g_r_p.between(g_r_p_min, g_r_p_max)
            & r_z_p.between(r_z_p_min, r_z_p_max)
        )
        valid_e = (
            g0_e.between(0.1, 29.9)
            & r0_e.between(0.1, 29.9)
            & z0_e.between(0.1, 29.9)
            & g_r_e.between(g_r_e_min, g_r_e_max)
            & r_z_e.between(r_z_e_min, r_z_e_max)
        )

        # We want to switch between psfmags and fibermags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, "sdss_psfmag_from_lsdr10"),
                ((ls.type != "PSF") & valid_e, "sdss_fiber2mag_from_lsdr10"),
            ),
            "undefined",
        )

        magnitude_g = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, g_p.cast("float")),
                ((ls.type != "PSF") & valid_e, g_e.cast("float")),
            ),
            "NaN",
        )
        magnitude_r = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, r_p.cast("float")),
                ((ls.type != "PSF") & valid_e, r_e.cast("float")),
            ),
            "NaN",
        )
        magnitude_i = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, i_p.cast("float")),
                ((ls.type != "PSF") & valid_e, i_e.cast("float")),
            ),
            "NaN",
        )
        magnitude_z = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, z_p.cast("float")),
                ((ls.type != "PSF") & valid_e, z_e.cast("float")),
            ),
            "NaN",
        )

        magnitude_gaia_g = peewee.Case(
            None, ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),), "NaN"
        )
        magnitude_gaia_bp = peewee.Case(
            None, ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),), "NaN"
        )
        magnitude_gaia_rp = peewee.Case(
            None, ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),), "NaN"
        )

        query = (
            c.select(
                c.catalogid.alias("catalogid"),
                priority.alias("priority"),
                value.alias("value"),
                cadence.alias("cadence"),
                instrument.alias("instrument"),
                opt_prov.alias("optical_prov"),
                magnitude_g.alias("g"),
                magnitude_r.alias("r"),
                magnitude_i.alias("i"),
                magnitude_z.alias("z"),
                magnitude_gaia_g.alias("gaia_g"),
                magnitude_gaia_bp.alias("bp"),
                magnitude_gaia_rp.alias("rp"),
                ls.ls_id.alias("ls_id"),  # extra
                ls.gaia_dr3_source_id.alias("gaia_dr3_source_id"),  # extra
                ls.gaia_dr2_source_id.alias("gaia_dr2_source_id"),  # extra
                x.ero_detuid.alias("ero_detuid"),  # extra
                x.ero_flux.alias("ero_flux"),  # extra
                x.ero_det_like.alias("ero_det_like"),  # extra
                x.ero_flags.alias("ero_flags"),  # extra
                x.xmatch_flags.alias("xmatch_flags"),  # extra
                x.xmatch_metric.alias("xmatch_metric"),  # extra
                s19.c.s19_pk.alias("sdss_dr19p_speclite_pk"),  # extra
                c.ra.alias("ra"),  # extra
                c.dec.alias("dec"),  # extra
                g0_p.alias("ls10_mag_g"),  # extra
                r0_p.alias("ls10_mag_r"),  # extra
                i0_p.alias("ls10_mag_i"),  # extra
                z0_p.alias("ls10_mag_z"),  # extra
                g0_e.alias("ls10_fibermag_g"),  # extra
                r0_e.alias("ls10_fibermag_r"),  # extra
                i0_e.alias("ls10_fibermag_i"),  # extra
                z0_e.alias("ls10_fibermag_z"),  # extra
                ls.nobs_g.alias("ls10_nobs_g"),
                ls.nobs_r.alias("ls10_nobs_r"),
                ls.nobs_i.alias("ls10_nobs_i"),
                ls.nobs_z.alias("ls10_nobs_z"),
                ls.type.alias("ls10_type"),  # extra
                ls.shape_r.alias("ls10_shape_r"),  # extra
                is_core.alias("is_core"),  # extra
                ls.ref_cat.alias("ls_ref_cat"),  # extra
                ls.ref_id.alias("ls_ref_id"),  # extra
                ls.maskbits.alias("ls_maskbits"),  # extra
                ls.fitbits.alias("ls_fitbits"),  # extra
                gal_lat.alias("abs_gal_lat"),  # extra
            )
            .join(c2ls)
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s19, JOIN.LEFT_OUTER)
            .join(s19, JOIN.LEFT_OUTER, on=(s19.c.s19_pk == c2s19.target_id))
            # finished joining the spectroscopy
            # admin criteria
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
                c2ls.best >> True,
                # fn.coalesce(c2s19.best, True) >> True,
            )
            # science criteria
            .where(
                (x.ero_version == self.parameters["ero_version"]),
                (x.xmatch_method == self.parameters["xmatch_method"]),
                (
                    (x.xmatch_version == self.parameters["xmatch_version1"])
                    | (x.xmatch_version == self.parameters["xmatch_version2"])
                ),
                (
                    (x.opt_cat == self.parameters["opt_cat1"])
                    | (x.opt_cat == self.parameters["opt_cat2"])
                ),
                (x.xmatch_metric >= self.parameters["p_any_min"]),
                (
                    (ls.fiberflux_r.between(fiberflux_r_min, fiberflux_r_max))
                    | (ls.fiberflux_i.between(fiberflux_i_min, fiberflux_i_max))
                    | (ls.fiberflux_z.between(fiberflux_z_min, fiberflux_z_max))
                ),
                (x.ero_det_like > self.parameters["det_like_min"]),
                # (ls.maskbits.bin_and(2**2 + 2**3 + 2**4) == 0),  # avoid saturated sources
                # avoid bright stars and globular clusters:
                # (ls.maskbits.bin_and(2**1 + 2**13) == 0),
                # always avoid very bright stars, see https://www.legacysurvey.org/dr10/bitmasks/
                (ls.maskbits.bin_and(2**1) == 0),
                # (ls.nobs_r > 0),                        # always require r-band coverage
                # ((ls.nobs_g > 0) | (ls.nobs_z > 0)),    # plus at least one other optical band
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters["gaia_g_mag_limit"])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters["gaia_rp_mag_limit"])),
            )
            # .group_by(ls)   # avoid duplicates - we trust the legacy survey entries
            .distinct([ls.ls_id])  # avoid duplicates - we trust the legacy survey entries
        )

        if self.only_faintest_cadence:
            query = query.where(cadence == cadence3)

        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    c.ra, c.dec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


#
# END BhmSpidersAgnLsdr10Carton
# ##################################################################################
class BhmSpidersAgnLsdr10D3Carton(BhmSpidersAgnLsdr10Carton):
    name = "bhm_spiders_agn_lsdr10_d3"
    only_faintest_cadence = True


# we can get away with just inheriting the selection code from
# the lsdr10 hemisphere match and adjusting the parameters only
# ##################################################################################
class BhmSpidersAgnHardCarton(BhmSpidersAgnLsdr10Carton):
    name = "bhm_spiders_agn_hard"


class BhmSpidersAgnHardD3Carton(BhmSpidersAgnLsdr10Carton):
    name = "bhm_spiders_agn_hard_d3"
    only_faintest_cadence = True


# # Testing of the North part of lsdr10 (i.e. dr9)
# # we can get away with just inheriting the selection code from
# # the lsdr10 hemisphere match and adjusting the parameters only
# # ##################################################################################
# class BhmSpidersAgnLsdr10NorthCarton(BhmSpidersAgnLsdr10Carton):
#     name = 'bhm_spiders_agn_lsdr10_north'


class BhmSpidersAgnGaiadr3Carton(BaseCarton):
    name = "bhm_spiders_agn_gaiadr3"
    category = "science"
    mapper = "BHM"
    program = "bhm_spiders"
    tile = False
    instrument = "BOSS"
    can_offset = True
    only_faintest_cadence = False

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        x = EROSITASupersetv1AGN.alias()
        g3 = Gaia_DR3.alias()
        c2g3 = CatalogToGaia_DR3.alias()

        instrument = peewee.Value(self.instrument)

        gaia_g_max_for_cadence1 = self.parameters["gaia_g_max_for_cadence1"]
        gaia_rp_max_for_cadence1 = self.parameters["gaia_rp_max_for_cadence1"]
        gaia_g_max_for_cadence2 = self.parameters["gaia_g_max_for_cadence2"]
        gaia_rp_max_for_cadence2 = self.parameters["gaia_rp_max_for_cadence2"]

        # these control matching to spectroscopy
        spec_sn_thresh = self.parameters["spec_sn_thresh"]
        spec_z_err_thresh = self.parameters["spec_z_err_thresh"]

        # #########################################################################
        # prepare the spectroscopy catalogues

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
        # #########################################################################

        # compute the abs(Galactic latitude):
        gal_lat = peewee.fn.abs(
            90.0 - peewee.fn.q3c_dist(north_gal_pole_ra, north_gal_pole_dec, c.ra, c.dec)
        )

        is_core = peewee.Case(
            None,
            (
                (gal_lat < self.parameters["min_gal_lat_for_core"], False),
                (c.dec < self.parameters["min_dec_for_core"], False),
                (x.ero_flux < self.parameters["min_ero_flux_for_core"], False),
                (x.ero_det_like < self.parameters["min_det_like_for_core"], False),
                (g3.phot_g_mean_mag < self.parameters["min_gaia_g_for_core"], False),
                (g3.phot_g_mean_mag > self.parameters["max_gaia_g_for_core"], False),
            ),
            True,
        )

        # value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        value = peewee.Case(
            None,
            #  ((gal_lat > self.parameters['in_plane_lat_cut'],
            #    self.parameters.get('value', 1.0)),),
            ((is_core, self.parameters.get("value", 1.0)),),
            0.0,
        ).cast("float")

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +dpriority_match_flags if target is a secondary cross-match (match_flag > 1)
        # add +dpriority_det_like if target has a low value of ero_det_like
        # add +dpriority_has_spec if target has existing good SDSS spectroscopy

        priority_1 = peewee.Case(None, ((x.xmatch_flags > 1, 1),), 0)
        priority_2 = peewee.Case(
            None, ((x.ero_det_like < self.parameters["det_like_for_priority"], 1),), 0
        )
        priority_3 = peewee.Case(None, ((s19.c.s19_pk.is_null(False), 1),), 0)
        priority_4 = peewee.Case(None, ((is_core, 0),), 1)
        priority_5 = peewee.Case(None, ((x.ero_flags.bin_and(2**9) == 0, 1),), 0)

        priority = (
            self.parameters["priority_floor"]
            + priority_1 * self.parameters["dpriority_match_flags"]
            + priority_2 * self.parameters["dpriority_det_like"]
            + priority_3 * self.parameters["dpriority_has_spec"]
            + priority_4 * self.parameters["dpriority_non_core"]
            + priority_5 * self.parameters["dpriority_not_hard"]
        )

        # choose cadence based on magnitude in Gaia G and RP-bands
        cadence1 = self.parameters["cadence1"]
        cadence2 = self.parameters["cadence2"]
        cadence3 = self.parameters["cadence3"]
        cadence4 = "unknown_cadence"
        cadence = peewee.Case(
            None,
            (
                (
                    (g3.phot_g_mean_mag < gaia_g_max_for_cadence1)
                    | (g3.phot_rp_mean_mag < gaia_rp_max_for_cadence1),
                    cadence1,
                ),
                (
                    (g3.phot_g_mean_mag < gaia_g_max_for_cadence2)
                    | (g3.phot_rp_mean_mag < gaia_rp_max_for_cadence2),
                    cadence2,
                ),
                (
                    (g3.phot_g_mean_mag >= gaia_g_max_for_cadence2)
                    & (g3.phot_rp_mean_mag >= gaia_rp_max_for_cadence2),
                    cadence3,
                ),
            ),
            cadence4,
        )

        # compute transformed SDSS mags
        # transform the Gaia dr3 G,BP,RP into sdss psfmag griz
        # piecewise transformation either side of BP-RP=1.8
        # fit to blue end is cubic, fit to red end is quadratic
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if(FILENAME~/_red/){pe="red"} else if (FILENAME~/_blue/){pe="blue"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_agn_gaiadr2_red/gdr2_*mag_to_sdss_*mag_?_results.log bhm_spiders_agn_gaiadr2_blue/gdr2_*mag_to_sdss_*mag_?_results.log  # noqa
        coeffs = {
            "g2_red": 0.081178,
            "g1_red": 0.355677,
            "g0_red": 0.510306,
            "i2_red": 0.048864,
            "i1_red": -0.287475,
            "i0_red": -0.336712,
            "r2_red": 0.028080,
            "r1_red": 0.542331,
            "r0_red": -1.055168,
            "z2_red": -0.131385,
            "z1_red": 0.302555,
            "z0_red": -1.381648,
            "g3_blue": 0.639054,
            "g2_blue": -1.739187,
            "g1_blue": 1.420330,
            "g0_blue": -0.194071,
            "i3_blue": 0.780585,
            "i2_blue": -2.549848,
            "i1_blue": 1.489880,
            "i0_blue": -0.241381,
            "r3_blue": 0.575494,
            "r2_blue": -2.077000,
            "r1_blue": 1.573302,
            "r0_blue": -0.295026,
            "z3_blue": 1.064986,
            "z2_blue": -3.162969,
            "z1_blue": 1.493750,
            "z0_blue": -0.199582,
        }

        g_blue = (
            g3.phot_g_mean_mag
            + coeffs["g0_blue"]
            + coeffs["g1_blue"] * g3.bp_rp
            + coeffs["g2_blue"] * g3.bp_rp * g3.bp_rp
            + coeffs["g3_blue"] * g3.bp_rp * g3.bp_rp * g3.bp_rp
        )
        r_blue = (
            g3.phot_g_mean_mag
            + coeffs["r0_blue"]
            + coeffs["r1_blue"] * g3.bp_rp
            + coeffs["r2_blue"] * g3.bp_rp * g3.bp_rp
            + coeffs["r3_blue"] * g3.bp_rp * g3.bp_rp * g3.bp_rp
        )
        i_blue = (
            g3.phot_g_mean_mag
            + coeffs["i0_blue"]
            + coeffs["i1_blue"] * g3.bp_rp
            + coeffs["i2_blue"] * g3.bp_rp * g3.bp_rp
            + coeffs["i3_blue"] * g3.bp_rp * g3.bp_rp * g3.bp_rp
        )
        z_blue = (
            g3.phot_g_mean_mag
            + coeffs["z0_blue"]
            + coeffs["z1_blue"] * g3.bp_rp
            + coeffs["z2_blue"] * g3.bp_rp * g3.bp_rp
            + coeffs["z3_blue"] * g3.bp_rp * g3.bp_rp * g3.bp_rp
        )

        g_red = (
            g3.phot_g_mean_mag
            + coeffs["g0_red"]
            + coeffs["g1_red"] * g3.bp_rp
            + coeffs["g2_red"] * g3.bp_rp * g3.bp_rp
        )
        r_red = (
            g3.phot_g_mean_mag
            + coeffs["r0_red"]
            + coeffs["r1_red"] * g3.bp_rp
            + coeffs["r2_red"] * g3.bp_rp * g3.bp_rp
        )
        i_red = (
            g3.phot_g_mean_mag
            + coeffs["i0_red"]
            + coeffs["i1_red"] * g3.bp_rp
            + coeffs["i2_red"] * g3.bp_rp * g3.bp_rp
        )
        z_red = (
            g3.phot_g_mean_mag
            + coeffs["z0_red"]
            + coeffs["z1_red"] * g3.bp_rp
            + coeffs["z2_red"] * g3.bp_rp * g3.bp_rp
        )

        # validity checks - set limits semi-manually
        bp_rp_min = 0.0
        bp_rp_max = 3.0
        valid = (
            g3.phot_g_mean_mag.between(0.1, 29.9)
            & g3.phot_bp_mean_mag.between(0.1, 29.9)
            & g3.phot_rp_mean_mag.between(0.1, 29.9)
            & g3.bp_rp.between(bp_rp_min, bp_rp_max)
        )
        opt_prov = peewee.Case(None, ((valid, "sdss_psfmag_from_gdr3"),), "undefined")
        magnitude_g = peewee.Case(
            None,
            (
                (valid & (g3.bp_rp < 1.8), g_blue),
                (valid & (g3.bp_rp > 1.8), g_red),
            ),
            "NaN",
        )
        magnitude_r = peewee.Case(
            None,
            (
                (valid & (g3.bp_rp < 1.8), r_blue),
                (valid & (g3.bp_rp > 1.8), r_red),
            ),
            "NaN",
        )
        magnitude_i = peewee.Case(
            None,
            (
                (valid & (g3.bp_rp < 1.8), i_blue),
                (valid & (g3.bp_rp > 1.8), i_red),
            ),
            "NaN",
        )
        magnitude_z = peewee.Case(
            None,
            (
                (valid & (g3.bp_rp < 1.8), z_blue),
                (valid & (g3.bp_rp > 1.8), z_red),
            ),
            "NaN",
        )

        query = (
            c.select(
                c.catalogid.alias("catalogid"),
                x.ero_detuid.alias("ero_detuid"),  # extra
                x.ero_flux.alias("ero_flux"),  # extra
                x.ero_det_like.alias("ero_det_like"),  # extra
                g3.source_id.alias("gaia_dr3_source_id"),  # extra
                s19.c.s19_pk.alias("sdss_dr19p_speclite_pk"),  # extra
                c.ra.alias("ra"),  # extra
                c.dec.alias("dec"),  # extra
                priority.alias("priority"),
                value.alias("value"),
                cadence.alias("cadence"),
                instrument.alias("instrument"),
                opt_prov.alias("optical_prov"),
                magnitude_g.alias("g"),
                magnitude_r.alias("r"),
                magnitude_i.alias("i"),
                magnitude_z.alias("z"),
                g3.phot_g_mean_mag.alias("gaia_g"),
                g3.phot_bp_mean_mag.alias("bp"),
                g3.phot_rp_mean_mag.alias("rp"),
                is_core.alias("is_core"),  # extra
                gal_lat.alias("abs_gal_lat"),  # extra
                x.xmatch_version.alias("xmatch_version"),  # extra
            )
            .join(c2g3)
            .where(
                c.version_id == version_id,
                c2g3.version_id == version_id,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
                c2g3.best >> True,
            )
            .join(g3)
            .join(x, on=(g3.source_id == x.gaia_dr3_source_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s19, JOIN.LEFT_OUTER)
            .join(s19, JOIN.LEFT_OUTER, on=(s19.c.s19_pk == c2s19.target_id))
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters["ero_version"]),
                (
                    (x.xmatch_method == self.parameters["xmatch_method1"])
                    | (x.xmatch_method == self.parameters["xmatch_method2"])
                ),
                (
                    (x.xmatch_version == self.parameters["xmatch_version1"])
                    | (x.xmatch_version == self.parameters["xmatch_version2"])
                ),
                (x.opt_cat == self.parameters["opt_cat"]),
                (x.xmatch_metric >= self.parameters["p_any_min"]),
                (g3.phot_g_mean_mag > self.parameters["gaia_g_mag_limit"]),
                (g3.phot_rp_mean_mag > self.parameters["gaia_rp_mag_limit"]),
                (x.ero_det_like > self.parameters["det_like_min"]),
            )
            .distinct([g3.source_id])  # avoid duplicates - we trust the gaia ids
        )

        if self.only_faintest_cadence:
            query = query.where(cadence == cadence3)

        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    c.ra, c.dec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


#
# END BhmSpidersAgnGaiadr3Carton


# For the version that goes via catwise202,
# we can get away with just inheriting the selection code from
# the gaia dr3 hemisphere match and adjusting the parameters only
# ##################################################################################
class BhmSpidersAgnGaiadr3viaCW2020Carton(BhmSpidersAgnGaiadr3Carton):
    name = "bhm_spiders_agn_gaiadr3_viacw2020"


# ##################################################################################
class BhmSpidersAgnGaiadr3BothCarton(BhmSpidersAgnGaiadr3Carton):
    name = "bhm_spiders_agn_gaiadr3_both"


# ##################################################################################
# we can get away with just inheriting the selection code from
# the gaia dr3 hemisphere match and adjusting the parameters only
class BhmSpidersAgnSepCarton(BhmSpidersAgnGaiadr3Carton):
    name = "bhm_spiders_agn_sep"


# ##################################################################################
class BhmSpidersAgnGaiadr3D3CartonCarton(BhmSpidersAgnGaiadr3Carton):
    name = "bhm_spiders_agn_gaiadr3_d3"
    only_faintest_cadence = True


# ##################################################################################
class BhmSpidersAgnSepD3Carton(BhmSpidersAgnGaiadr3Carton):
    name = "bhm_spiders_agn_sep_d3"
    only_faintest_cadence = True


class BhmSpidersAgnTdaCarton(BaseCarton):
    name = "bhm_spiders_agn_tda"
    category = "science"
    mapper = "BHM"
    program = "bhm_spiders"
    tile = False
    instrument = "BOSS"
    can_offset = True
    only_faintest_cadence = False

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        x = EROSITASupersetv1AGN.alias()
        ls = Legacy_Survey_DR10.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()

        instrument = peewee.Value(self.instrument)

        fiberflux_r_max = AB2nMgy(self.parameters["fibermag_r_min"])
        fiberflux_r_min = AB2nMgy(self.parameters["fibermag_r_max"])
        fiberflux_i_max = AB2nMgy(self.parameters["fibermag_i_min"])
        fiberflux_i_min = AB2nMgy(self.parameters["fibermag_i_max"])
        fiberflux_z_max = AB2nMgy(self.parameters["fibermag_z_min"])
        fiberflux_z_min = AB2nMgy(self.parameters["fibermag_z_max"])
        fiberflux_r_min_for_cadence1 = AB2nMgy(self.parameters["fibermag_r_for_cadence1"])
        fiberflux_r_min_for_cadence2 = AB2nMgy(self.parameters["fibermag_r_for_cadence2"])
        gaia_g_max_for_cadence1 = self.parameters["gaia_g_max_for_cadence1"]
        gaia_rp_max_for_cadence1 = self.parameters["gaia_rp_max_for_cadence1"]

        # choose cadence based on fiber magnitude in r-band + gaia G,RP
        cadence1 = self.parameters["cadence1"]
        cadence2 = self.parameters["cadence2"]
        cadence3 = self.parameters["cadence3"]
        cadence4 = "unknown_cadence"
        cadence = peewee.Case(
            None,
            (
                (
                    (
                        (ls.fiberflux_r > fiberflux_r_min_for_cadence1)
                        | (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1))
                        | (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))
                    ),
                    cadence1,
                ),
                (ls.fiberflux_r > fiberflux_r_min_for_cadence2, cadence2),
                (ls.fiberflux_r <= fiberflux_r_min_for_cadence2, cadence3),
            ),
            cadence4,
        )

        value = peewee.Value(self.parameters["value"]).cast("float")
        priority = peewee.Value(self.parameters["priority_floor"]).cast("integer")

        # compute transformed SDSS mags for pointlike and extended sources separately
        # transform the legacysurvey grz into sdss psfmag griz

        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_agn_lsdr8_*/lsdr8_*mag_to_sdss_*mag_?_results.log  # noqa
        coeffs = {
            "g2_e": -0.113816,
            "g1_e": 0.317176,
            "g0_e": 0.094145,
            "i2_e": -0.415858,
            "i1_e": 0.168922,
            "i0_e": -0.010771,
            "r2_e": 0.029398,
            "r1_e": -0.019938,
            "r0_e": 0.354042,
            "z2_e": -0.111262,
            "z1_e": 0.237656,
            "z0_e": 0.148923,
            "g2_p": 0.187193,
            "g1_p": -0.184362,
            "g0_p": 0.049492,
            "i2_p": -0.098979,
            "i1_p": -0.405518,
            "i0_p": 0.009688,
            "r2_p": -0.001935,
            "r1_p": 0.098201,
            "r0_p": 0.050321,
            "z2_p": -0.034163,
            "z1_p": 0.109878,
            "z0_p": -0.030167,
        }

        nMgy_min = 1e-3  # equiv to AB=30
        # pointlike - start from ls8 (psf)fluxes
        g0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g))
        r0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r))
        i0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_i))
        z0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z))
        g_r_p = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g) / peewee.fn.greatest(nMgy_min, ls.flux_r)
        )
        r_z_p = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r) / peewee.fn.greatest(nMgy_min, ls.flux_z)
        )

        # extended - start from ls8 fiberfluxes
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

        g_p = g0_p + coeffs["g0_p"] + coeffs["g1_p"] * g_r_p + coeffs["g2_p"] * g_r_p * g_r_p
        r_p = r0_p + coeffs["r0_p"] + coeffs["r1_p"] * g_r_p + coeffs["r2_p"] * g_r_p * g_r_p
        i_p = r0_p + coeffs["i0_p"] + coeffs["i1_p"] * r_z_p + coeffs["i2_p"] * r_z_p * r_z_p
        z_p = z0_p + coeffs["z0_p"] + coeffs["z1_p"] * r_z_p + coeffs["z2_p"] * r_z_p * r_z_p

        g_e = g0_e + coeffs["g0_e"] + coeffs["g1_e"] * g_r_e + coeffs["g2_e"] * g_r_e * g_r_e
        r_e = r0_e + coeffs["r0_e"] + coeffs["r1_e"] * g_r_e + coeffs["r2_e"] * g_r_e * g_r_e
        i_e = r0_e + coeffs["i0_e"] + coeffs["i1_e"] * r_z_e + coeffs["i2_e"] * r_z_e * r_z_e
        z_e = z0_e + coeffs["z0_e"] + coeffs["z1_e"] * r_z_e + coeffs["z2_e"] * r_z_e * r_z_e

        # validity checks - set limits semi-manually
        g_r_p_min = -0.25
        g_r_p_max = 1.75
        r_z_p_min = -0.5
        r_z_p_max = 2.5
        g_r_e_min = 0.0
        g_r_e_max = 1.75
        r_z_e_min = 0.2
        r_z_e_max = 1.6
        valid_p = (
            g0_p.between(0.1, 29.9)
            & r0_p.between(0.1, 29.9)
            & z0_p.between(0.1, 29.9)
            & g_r_p.between(g_r_p_min, g_r_p_max)
            & r_z_p.between(r_z_p_min, r_z_p_max)
        )
        valid_e = (
            g0_e.between(0.1, 29.9)
            & r0_e.between(0.1, 29.9)
            & z0_e.between(0.1, 29.9)
            & g_r_e.between(g_r_e_min, g_r_e_max)
            & r_z_e.between(r_z_e_min, r_z_e_max)
        )

        # We want to switch between psfmags and fibermags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, "sdss_psfmag_from_lsdr10"),
                ((ls.type != "PSF") & valid_e, "sdss_fiber2mag_from_lsdr10"),
            ),
            "undefined",
        )

        magnitude_g = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, g_p.cast("float")),
                ((ls.type != "PSF") & valid_e, g_e.cast("float")),
            ),
            "NaN",
        )
        magnitude_r = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, r_p.cast("float")),
                ((ls.type != "PSF") & valid_e, r_e.cast("float")),
            ),
            "NaN",
        )
        magnitude_i = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, i_p.cast("float")),
                ((ls.type != "PSF") & valid_e, i_e.cast("float")),
            ),
            "NaN",
        )
        magnitude_z = peewee.Case(
            None,
            (
                ((ls.type == "PSF") & valid_p, z_p.cast("float")),
                ((ls.type != "PSF") & valid_e, z_e.cast("float")),
            ),
            "NaN",
        )

        magnitude_gaia_g = peewee.Case(
            None, ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),), "NaN"
        )
        magnitude_gaia_bp = peewee.Case(
            None, ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),), "NaN"
        )
        magnitude_gaia_rp = peewee.Case(
            None, ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),), "NaN"
        )

        query = (
            c.select(
                c.catalogid.alias("catalogid"),
                priority.alias("priority"),
                value.alias("value"),
                cadence.alias("cadence"),
                instrument.alias("instrument"),
                opt_prov.alias("optical_prov"),
                magnitude_g.alias("g"),
                magnitude_r.alias("r"),
                magnitude_i.alias("i"),
                magnitude_z.alias("z"),
                magnitude_gaia_g.alias("gaia_g"),
                magnitude_gaia_bp.alias("bp"),
                magnitude_gaia_rp.alias("rp"),
                ls.ls_id.alias("ls_id"),  # extra
                ls.gaia_dr3_source_id.alias("gaia_dr3_source_id"),  # extra
                ls.gaia_dr2_source_id.alias("gaia_dr2_source_id"),  # extra
                x.ero_detuid.alias("ero_detuid"),  # extra
                x.ero_flux.alias("ero_flux"),  # extra
                x.ero_det_like.alias("ero_det_like"),  # extra
                x.ero_flags.alias("ero_flags"),  # extra
                x.xmatch_flags.alias("xmatch_flags"),  # extra
                x.xmatch_metric.alias("xmatch_metric"),  # extra
                c.ra.alias("ra"),  # extra
                c.dec.alias("dec"),  # extra
                g0_p.alias("ls10_mag_g"),  # extra
                r0_p.alias("ls10_mag_r"),  # extra
                i0_p.alias("ls10_mag_i"),  # extra
                z0_p.alias("ls10_mag_z"),  # extra
                g0_e.alias("ls10_fibermag_g"),  # extra
                r0_e.alias("ls10_fibermag_r"),  # extra
                i0_e.alias("ls10_fibermag_i"),  # extra
                z0_e.alias("ls10_fibermag_z"),  # extra
                ls.nobs_g.alias("ls10_nobs_g"),
                ls.nobs_r.alias("ls10_nobs_r"),
                ls.nobs_i.alias("ls10_nobs_i"),
                ls.nobs_z.alias("ls10_nobs_z"),
                ls.type.alias("ls10_type"),  # extra
                ls.shape_r.alias("ls10_shape_r"),  # extra
                ls.ref_cat.alias("ls_ref_cat"),  # extra
                ls.ref_id.alias("ls_ref_id"),  # extra
                ls.maskbits.alias("ls_maskbits"),  # extra
                ls.fitbits.alias("ls_fitbits"),  # extra
                # gal_lat.alias('abs_gal_lat'),  # extra
            )
            .join(c2ls)
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            # admin criteria
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True,
            )
            # science criteria
            .where(
                (x.ero_version == self.parameters["ero_version"]),
                (x.xmatch_method == self.parameters["xmatch_method"]),
                (x.xmatch_version == self.parameters["xmatch_version"]),
                (x.opt_cat == self.parameters["opt_cat"]),
                (
                    (ls.fiberflux_r.between(fiberflux_r_min, fiberflux_r_max))
                    | (ls.fiberflux_i.between(fiberflux_i_min, fiberflux_i_max))
                    | (ls.fiberflux_z.between(fiberflux_z_min, fiberflux_z_max))
                ),
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters["gaia_g_mag_limit"])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters["gaia_rp_mag_limit"])),
            )
            # .group_by(ls)   # avoid duplicates - we trust the legacy survey entries
            .distinct([ls.ls_id])  # avoid duplicates - we trust the legacy survey entries
        )

        if self.only_faintest_cadence:
            query = query.where(cadence == cadence3)

        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    c.ra, c.dec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


# ##################################################################################
class BhmSpidersAgnTdaD3Carton(BhmSpidersAgnTdaCarton):
    name = "bhm_spiders_agn_tda_d3"
    only_faintest_cadence = True
