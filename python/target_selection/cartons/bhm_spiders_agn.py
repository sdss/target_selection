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
    EROSITASupersetAGN,
)

# imports of existing spectro catalogue
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToSDSS_DR19p_Speclite,
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

## additional imports required by bhm_spiders_agn_ps1dr2
#from sdssdb.peewee.sdss5db.catalogdb import (
#    Panstarrs1,
#    CatalogToPanstarrs1,
#)

from target_selection.mag_flux import AB2nMgy, AB2Jy


# used by cartons that need to compute Galactic latitude:
north_gal_pole_ra = 192.85948   # deg, J2000
north_gal_pole_dec = +27.12825   # deg, J2000

# ############################################
# ############################################
# ############################################
# ############################################
# This provides the following BHM SPIDERS AGN cartons in v1.0:
#  *  bhm_spiders_agn_lsdr10
#  *  bhm_spiders_agn_gaiadr3
#  *  bhm_spiders_agn_sep
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

# Notes on how many targets to expect:
# sdss5db=> SELECT ero_version,xmatch_method,xmatch_version,opt_cat,count(*)
#           FROM erosita_superset_agn GROUP BY ero_version,xmatch_method,xmatch_version,opt_cat;
#        ero_version        | xmatch_method  |      xmatch_version      |   opt_cat    | count
# --------------------------+----------------+--------------------------+--------------+--------
#  eFEDS_c001_V18C_V3_ext   | XPS/MLR        | Merged_03DEC2020         | lsdr8        |     14
#  eFEDS_c001_V18C_V3_ext   | XPS/NWAY       | Merged_03DEC2020         | lsdr8        |    248
#  eFEDS_c001_V18C_V3_main  | XPS/MLR        | Merged_03DEC2020         | lsdr8        |    794
#  eFEDS_c001_V18C_V3_main  | XPS/NWAY       | Merged_03DEC2020         | lsdr8        |  26575
#  em01_c946_201008_poscorr | XPS/NWAY       | JWMS_CW2_v_03_TDopt      | gaiadr2      | 441175
#  em01_c946_201008_poscorr | XPS/NWAY       | JWMS_CW2_v_03_TDopt      | lsdr8        | 305076
#  em01_c946_201008_poscorr | XPS/NWAY       | JWMS_CW2_v_03_TDopt      | ps1dr2       | 241150
#  em01_c946_201008_poscorr | XPS/NWAY       | JWMS_CW2_v_03_TDopt      | skymapperdr2 | 312372
#  em01_c946_201008_poscorr | XPS/NWAY       | JWMS_v_03                | catwise2020  | 740691
#  em01_c946_201008_poscorr | XPS/NWAY       | JWMS_v_40                | lsdr8        | 345189
#  em01_SEP_c946            | XPS/NWAY       | SEP_CW2_07DEC2020        | catwise2020  |  32268
#  em01_SEP_c946            | XPS/NWAY       | SEP_CW2_07DEC2020_TDopt  | gaiadr2      |    309
#  em01_SEP_c946            | XPS/NWAY       | SEP_CW2_NOV2020_MSopt    | gaiadr2      |    740
# (13 rows)


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

    name = 'bhm_spiders_agn_lsdr10'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        ls = Legacy_Survey_DR10.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()

        instrument = peewee.Value(self.instrument)

        fibertotflux_r_max = AB2nMgy(self.parameters['fibertotmag_r_min'])
        fibertotflux_r_min = AB2nMgy(self.parameters['fibertotmag_r_max'])
        fibertotflux_z_max = AB2nMgy(self.parameters['fibertotmag_z_min'])
        fibertotflux_z_min = AB2nMgy(self.parameters['fibertotmag_z_max'])

        fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence1'])
        fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence2'])
        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']

        # #########################################################################
        # prepare the spectroscopy catalogue
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # SDSS DR19p
        # downslect only 'good' spectra
        c2s19 = CatalogToSDSS_DR19p_Speclite.alias()
        ss19 = SDSS_DR19p_Speclite.alias()
        s19 = (
            ss19.select(
                ss19.pk.alias('s19_pk'),
            )
            .where(
                ss19.snmedian >= spec_sn_thresh,
                ss19.zwarning == 0,
                ss19.zerr <= spec_z_err_thresh,
                ss19.zerr > 0.0,
                ss19.scienceprimary > 0,
            )
            .alias('s19')
        )
        # #########################################################################

        # compute the abs(Galactic latitude):
        gal_lat = peewee.fn.abs(90.0 - peewee.fn.q3c_dist(north_gal_pole_ra,
                                                          north_gal_pole_dec,
                                                          c.ra, c.dec))

        # value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        value = peewee.Case(
            None,
            ((gal_lat > self.parameters['in_plane_lat_cut'],
              self.parameters.get('value', 1.0)),),
            0.0
        ).cast('float')

        # priority is determined from individual target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +dpriority_match_flags if target is a secondary cross-match (match_flag > 1)
        # add +dpriority_det_like if target has a low value of ero_det_like
        # add +dpriority_has_spec if target has existing good SDSS spectroscopy

        priority_1 = peewee.Case(
            None,
            ((x.xmatch_flags > 1, 1), ),
            0)
        priority_2 = peewee.Case(
            None,
            ((x.ero_det_like < self.parameters['det_like_for_priority'], 1), ),
            0)
        priority_3 = peewee.Case(
            None,
            ((s19.c.s19_pk.is_null(False), 1), ),
            0)
        priority_4 = peewee.Case(
            None,
            ((gal_lat > self.parameters['in_plane_lat_cut'], 0),),
            1
        )

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec'] +
            priority_4 * self.parameters['dpriority_in_plane']
        )

        # choose cadence based on fiber magnitude in r-band
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                (
                    ((ls.fibertotflux_r > fibertotflux_r_min_for_cadence1) |
                     (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1)) |
                     (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))),
                    cadence1),
                (ls.fibertotflux_r > fibertotflux_r_min_for_cadence2, cadence2),
                (ls.fibertotflux_r <= fibertotflux_r_min_for_cadence2, cadence3),
            ),
            cadence4)

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
        g0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g)))
        r0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r)))
        i0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_i)))
        z0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z)))
        g_r_p = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g) /
                                      peewee.fn.greatest(nMgy_min, ls.flux_r)))
        r_z_p = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r) /
                                      peewee.fn.greatest(nMgy_min, ls.flux_z)))

        # extended - start from ls8 fiberfluxes
        g0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g)))
        r0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r)))
        i0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_i)))
        z0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_z)))
        g_r_e = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g) /
                                      peewee.fn.greatest(nMgy_min, ls.fiberflux_r)))
        r_z_e = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r) /
                                      peewee.fn.greatest(nMgy_min, ls.fiberflux_z)))

        g_p = (g0_p + coeffs['g0_p'] + coeffs['g1_p'] * g_r_p + coeffs['g2_p'] * g_r_p * g_r_p)
        r_p = (r0_p + coeffs['r0_p'] + coeffs['r1_p'] * g_r_p + coeffs['r2_p'] * g_r_p * g_r_p)
        i_p = (r0_p + coeffs['i0_p'] + coeffs['i1_p'] * r_z_p + coeffs['i2_p'] * r_z_p * r_z_p)
        z_p = (z0_p + coeffs['z0_p'] + coeffs['z1_p'] * r_z_p + coeffs['z2_p'] * r_z_p * r_z_p)

        g_e = (g0_e + coeffs['g0_e'] + coeffs['g1_e'] * g_r_e + coeffs['g2_e'] * g_r_e * g_r_e)
        r_e = (r0_e + coeffs['r0_e'] + coeffs['r1_e'] * g_r_e + coeffs['r2_e'] * g_r_e * g_r_e)
        i_e = (r0_e + coeffs['i0_e'] + coeffs['i1_e'] * r_z_e + coeffs['i2_e'] * r_z_e * r_z_e)
        z_e = (z0_e + coeffs['z0_e'] + coeffs['z1_e'] * r_z_e + coeffs['z2_e'] * r_z_e * r_z_e)

        # validity checks - set limits semi-manually
        g_r_p_min = -0.25
        g_r_p_max = 1.75
        r_z_p_min = -0.5
        r_z_p_max = 2.5
        g_r_e_min = 0.0
        g_r_e_max = 1.75
        r_z_e_min = 0.2
        r_z_e_max = 1.6
        valid_p = (g0_p.between(0.1, 29.9) &
                   r0_p.between(0.1, 29.9) &
                   z0_p.between(0.1, 29.9) &
                   g_r_p.between(g_r_p_min, g_r_p_max) &
                   r_z_p.between(r_z_p_min, r_z_p_max))
        valid_e = (g0_e.between(0.1, 29.9) &
                   r0_e.between(0.1, 29.9) &
                   z0_e.between(0.1, 29.9) &
                   g_r_e.between(g_r_e_min, g_r_e_max) &
                   r_z_e.between(r_z_e_min, r_z_e_max))

        # We want to switch between psfmags and fibertotmags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            None,
            (
                ((ls.type == 'PSF') & valid_p, 'sdss_psfmag_from_lsdr10'),
                ((ls.type != 'PSF') & valid_e, 'sdss_fiber2mag_from_lsdr10'),
            ),
            'undefined')

        magnitude_g = peewee.Case(
            None,
            (
                ((ls.type == 'PSF') & valid_p, g_p.cast('float')),
                ((ls.type != 'PSF') & valid_e, g_e.cast('float')),
            ),
            'NaN')
        magnitude_r = peewee.Case(
            None,
            (
                ((ls.type == 'PSF') & valid_p, r_p.cast('float')),
                ((ls.type != 'PSF') & valid_e, r_e.cast('float')),
            ),
            'NaN')
        magnitude_i = peewee.Case(
            None,
            (
                ((ls.type == 'PSF') & valid_p, i_p.cast('float')),
                ((ls.type != 'PSF') & valid_e, i_e.cast('float')),
            ),
            'NaN')
        magnitude_z = peewee.Case(
            None,
            (
                ((ls.type == 'PSF') & valid_p, z_p.cast('float')),
                ((ls.type != 'PSF') & valid_e, z_e.cast('float')),
            ),
            'NaN')

        magnitude_gaia_g = peewee.Case(
            None,
            ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),),
            'NaN')
        magnitude_gaia_bp = peewee.Case(
            None,
            ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),),
            'NaN')
        magnitude_gaia_rp = peewee.Case(
            None,
            ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),),
            'NaN')

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(ls.ls_id).alias('ls_id'),  # extra
                fn.min(x.ero_detuid).alias('ero_detuid'),  # extra
                fn.min(s19.c.s19_pk).alias('sdss_dr19p_speclite_pk'), # extra
                fn.min(c.ra).alias('ra'),   # extra
                fn.min(c.dec).alias('dec'),   # extra
                priority.alias('priority'),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(opt_prov).alias('optical_prov'),
                fn.min(magnitude_g).alias('g'),
                fn.min(magnitude_r).alias('r'),
                fn.min(magnitude_i).alias('i'),
                fn.min(magnitude_z).alias('z'),
                fn.min(magnitude_gaia_g).alias('gaia_g'),
                fn.min(magnitude_gaia_bp).alias('bp'),
                fn.min(magnitude_gaia_rp).alias('rp'),
                fn.min(g0_p).alias('ls10_mag_g'),   # extra
                fn.min(r0_p).alias('ls10_mag_r'),  # extra
                fn.min(i0_p).alias('ls10_mag_i'),  # extra
                fn.min(z0_p).alias('ls10_mag_z'),  # extra
                fn.min(g0_e).alias('ls10_fibermag_g'),  # extra
                fn.min(r0_e).alias('ls10_fibermag_r'),  # extra
                fn.min(i0_e).alias('ls10_fibermag_i'),  # extra
                fn.min(z0_e).alias('ls10_fibermag_z'),  # extra
                fn.min(ls.type).alias('ls10_type'),  # extra
            )
            .join(c2ls)
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
                c2ls.best >> True,
            )
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s19, JOIN.LEFT_OUTER)
            .join(s19, JOIN.LEFT_OUTER,
                  on=(s19.c.s19_pk == c2s19.target_id))
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (x.xmatch_version == self.parameters['xmatch_version']),
                (x.opt_cat == self.parameters['opt_cat']),
                (x.xmatch_metric >= self.parameters['p_any_min']),
                (
                    (ls.fibertotflux_r.between(fibertotflux_r_min, fibertotflux_r_max)) |
                    (ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max))
                ),
                (x.ero_det_like > self.parameters['det_like_min']),
                (ls.maskbits.bin_and(2**2 + 2**3 + 2**4) == 0),  # avoid saturated sources
                (ls.nobs_r > 0),                        # always require r-band coverage
                ((ls.nobs_g > 0) | (ls.nobs_z > 0)),    # plus at least one other optical band
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters['gaia_g_mag_limit'])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters['gaia_rp_mag_limit'])),
            )
            .group_by(ls)   # avoid duplicates - we trust the legacy survey entries
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query

#
# END BhmSpidersAgnLsdr10Carton
# ##################################################################################


# class BhmSpidersAgnEfedsStragglersCarton(BaseCarton):
#
#     name = 'bhm_spiders_agn_efeds_stragglers'
#     category = 'science'
#     mapper = 'BHM'
#     program = 'bhm_spiders'
#     tile = False
#     instrument = 'BOSS'
#     cadence = None
#
#     def build_query(self, version_id, query_region=None):
#
#         c = Catalog.alias()
#         x = EROSITASupersetAGN.alias()
#         ls = Legacy_Survey_DR8.alias()
#         c2ls = CatalogToLegacy_Survey_DR8.alias()
#
#         instrument = peewee.Value(self.instrument)
#
#         fibertotflux_r_max = AB2nMgy(self.parameters['fibertotmag_r_min'])
#         fibertotflux_r_min = AB2nMgy(self.parameters['fibertotmag_r_max'])
#         fibertotflux_z_max = AB2nMgy(self.parameters['fibertotmag_z_min'])
#         fibertotflux_z_min = AB2nMgy(self.parameters['fibertotmag_z_max'])
#
#         fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence1'])
#         fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence2'])
#         gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
#         gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']
#
#         # flux30 = AB2nMgy(30.00)
#         value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
#
#         match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
#         spec_sn_thresh = self.parameters['spec_sn_thresh']
#         spec_z_err_thresh = self.parameters['spec_z_err_thresh']
#
#         #
#         ##########################################################################
#         # prepare the spectroscopy catalogues
#
#         # SDSS DR16
#         c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
#         ss16 = SDSS_DR16_SpecObj.alias()
#         s16 = (
#             ss16.select(
#                 ss16.specobjid.alias('specobjid'),
#             )
#             .where(
#                 ss16.snmedian >= spec_sn_thresh,
#                 ss16.zwarning == 0,
#                 ss16.zerr <= spec_z_err_thresh,
#                 ss16.zerr > 0.0,
#                 ss16.scienceprimary > 0,
#             )
#             .alias('s16')
#         )
#
#         # SDSS-IV/eFEDS March2020
#         c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
#         ss2020 = BHM_eFEDS_Veto.alias()
#         s2020 = (
#             ss2020.select(
#                 ss2020.pk.alias('pk'),
#             )
#             .where(
#                 ss2020.sn_median_all >= spec_sn_thresh,
#                 ss2020.zwarning == 0,
#                 ss2020.z_err <= spec_z_err_thresh,
#                 ss2020.z_err > 0.0,
#             )
#             .alias('s2020')
#         )
#
#         # SDSS-V spAll
#         ssV = SDSSV_BOSS_SPALL.alias()
#         sV = (
#             ssV.select(
#                 ssV.specobjid.alias('specobjid'),
#                 ssV.plug_ra.alias('plug_ra'),
#                 ssV.plug_dec.alias('plug_dec'),
#             )
#             .where(
#                 ssV.sn_median_all >= spec_sn_thresh,
#                 ssV.zwarning == 0,
#                 ssV.z_err <= spec_z_err_thresh,
#                 ssV.z_err > 0.0,
#                 ssV.specprimary > 0,
#             )
#             .alias('sV')
#         )
#
#         # priority is determined by target properties
#         # start with a priority floor value (per carton)
#         # then increment if any conditions are met:
#         # add +dpriority_match_flags if target is a secondary cross-match (match_flag > 1)
#         # add +dpriority_det_like if target has a low value of ero_det_like
#         # add +dpriority_ero_version if target is from the secondary eFEDS catalogue
#         # add +dpriority_has_spec if target has existing good SDSS spectroscopy
#
#         priority_1 = peewee.Case(
#             None,
#             ((x.xmatch_flags > 1, 1), ),
#             0)
#         priority_2 = peewee.Case(
#             None,
#             ((x.ero_det_like < self.parameters['det_like_for_priority'], 1), ),
#             0)
#         # add a step down in priority for anything only selected by the secondary xmatch_version
#         priority_3 = peewee.Case(
#             None,
#             ((x.ero_version == self.parameters['ero_version2'], 1), ),
#             0)
#
#         # no need for the step below because we just reject everything with a spectrum.
#         priority_4 = peewee.Case(
#             None,
#             (
#                 (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
#                 (s2020.c.pk.is_null(False), 1),
#                 (sV.c.specobjid.is_null(False), 1),
#                 # redundant # (sph.c.pkey.is_null(False), 1),
#             ),
#             0)
#
#         priority = (
#             self.parameters['priority_floor'] +
#             fn.min(priority_1) * self.parameters['dpriority_match_flags'] +
#             fn.min(priority_2) * self.parameters['dpriority_det_like'] +
#             fn.min(priority_3) * self.parameters['dpriority_ero_version'] +
#             fn.max(priority_4) * self.parameters['dpriority_has_spec']
#         )
#
#         # choose cadence based on fiber magnitude in r-band
#         cadence1 = self.parameters['cadence1']
#         cadence2 = self.parameters['cadence2']
#         cadence3 = self.parameters['cadence3']
#         cadence4 = 'unknown_cadence'
#         cadence = peewee.Case(
#             None,
#             (
#                 (
#                     ((ls.fibertotflux_r > fibertotflux_r_min_for_cadence1) |
#                      (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1)) |
#                      (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))),
#                     cadence1),
#                 (ls.fibertotflux_r > fibertotflux_r_min_for_cadence2, cadence2),
#                 (ls.fibertotflux_r <= fibertotflux_r_min_for_cadence2, cadence3),
#             ),
#             cadence4)
#
#         # compute transformed SDSS mags for pointlike and extended sources separately
#         # transform the legacysurvey grz into sdss psfmag griz
#
#         # extract coeffs from fit logs via:
#         # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_agn_lsdr8_*/lsdr8_*mag_to_sdss_*mag_?_results.log  # noqa
#         coeffs = {
#             "g2_e": -0.113816,
#             "g1_e": 0.317176,
#             "g0_e": 0.094145,
#             "i2_e": -0.415858,
#             "i1_e": 0.168922,
#             "i0_e": -0.010771,
#             "r2_e": 0.029398,
#             "r1_e": -0.019938,
#             "r0_e": 0.354042,
#             "z2_e": -0.111262,
#             "z1_e": 0.237656,
#             "z0_e": 0.148923,
#             "g2_p": 0.187193,
#             "g1_p": -0.184362,
#             "g0_p": 0.049492,
#             "i2_p": -0.098979,
#             "i1_p": -0.405518,
#             "i0_p": 0.009688,
#             "r2_p": -0.001935,
#             "r1_p": 0.098201,
#             "r0_p": 0.050321,
#             "z2_p": -0.034163,
#             "z1_p": 0.109878,
#             "z0_p": -0.030167,
#         }
#
#         nMgy_min = 1e-3  # equiv to AB=30
#         # pointlike - start from ls8 (psf)fluxes
#         g0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g)))
#         r0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r)))
#         z0_p = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z)))
#         g_r_p = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g) /
#                                       peewee.fn.greatest(nMgy_min, ls.flux_r)))
#         r_z_p = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r) /
#                                       peewee.fn.greatest(nMgy_min, ls.flux_z)))
#
#         # extended - start from ls8 fiberfluxes
#         g0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g)))
#         r0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r)))
#         z0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_z)))
#         g_r_e = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g) /
#                                       peewee.fn.greatest(nMgy_min, ls.fiberflux_r)))
#         r_z_e = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r) /
#                                       peewee.fn.greatest(nMgy_min, ls.fiberflux_z)))
#
#         g_p = (g0_p + coeffs['g0_p'] + coeffs['g1_p'] * g_r_p + coeffs['g2_p'] * g_r_p * g_r_p)
#         r_p = (r0_p + coeffs['r0_p'] + coeffs['r1_p'] * g_r_p + coeffs['r2_p'] * g_r_p * g_r_p)
#         i_p = (r0_p + coeffs['i0_p'] + coeffs['i1_p'] * r_z_p + coeffs['i2_p'] * r_z_p * r_z_p)
#         z_p = (z0_p + coeffs['z0_p'] + coeffs['z1_p'] * r_z_p + coeffs['z2_p'] * r_z_p * r_z_p)
#
#         g_e = (g0_e + coeffs['g0_e'] + coeffs['g1_e'] * g_r_e + coeffs['g2_e'] * g_r_e * g_r_e)
#         r_e = (r0_e + coeffs['r0_e'] + coeffs['r1_e'] * g_r_e + coeffs['r2_e'] * g_r_e * g_r_e)
#         i_e = (r0_e + coeffs['i0_e'] + coeffs['i1_e'] * r_z_e + coeffs['i2_e'] * r_z_e * r_z_e)
#         z_e = (z0_e + coeffs['z0_e'] + coeffs['z1_e'] * r_z_e + coeffs['z2_e'] * r_z_e * r_z_e)
#
#         # validity checks - set limits semi-manually
#         g_r_p_min = -0.25
#         g_r_p_max = 1.75
#         r_z_p_min = -0.5
#         r_z_p_max = 2.5
#         g_r_e_min = 0.0
#         g_r_e_max = 1.75
#         r_z_e_min = 0.2
#         r_z_e_max = 1.6
#         valid_p = (g0_p.between(0.1, 29.9) &
#                    r0_p.between(0.1, 29.9) &
#                    z0_p.between(0.1, 29.9) &
#                    g_r_p.between(g_r_p_min, g_r_p_max) &
#                    r_z_p.between(r_z_p_min, r_z_p_max))
#         valid_e = (g0_e.between(0.1, 29.9) &
#                    r0_e.between(0.1, 29.9) &
#                    z0_e.between(0.1, 29.9) &
#                    g_r_e.between(g_r_e_min, g_r_e_max) &
#                    r_z_e.between(r_z_e_min, r_z_e_max))
#
#         # We want to switch between psfmags and fibertotmags depending on
#         # ls.type parameter (PSF or extended)
#         # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
#         opt_prov = peewee.Case(
#             None,
#             (
#                 ((ls.type == 'PSF') & valid_p, 'sdss_psfmag_from_lsdr8'),
#                 ((ls.type != 'PSF') & valid_e, 'sdss_fiber2mag_from_lsdr8'),
#             ),
#             'undefined')
#
#         magnitude_g = peewee.Case(
#             None,
#             (
#                 ((ls.type == 'PSF') & valid_p, g_p.cast('float')),
#                 ((ls.type != 'PSF') & valid_e, g_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_r = peewee.Case(
#             None,
#             (
#                 ((ls.type == 'PSF') & valid_p, r_p.cast('float')),
#                 ((ls.type != 'PSF') & valid_e, r_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_i = peewee.Case(
#             None,
#             (
#                 ((ls.type == 'PSF') & valid_p, i_p.cast('float')),
#                 ((ls.type != 'PSF') & valid_e, i_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_z = peewee.Case(
#             None,
#             (
#                 ((ls.type == 'PSF') & valid_p, z_p.cast('float')),
#                 ((ls.type != 'PSF') & valid_e, z_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_gaia_g = peewee.Case(
#             None,
#             ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),),
#             'NaN')
#         magnitude_gaia_bp = peewee.Case(
#             None,
#             ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),),
#             'NaN')
#         magnitude_gaia_rp = peewee.Case(
#             None,
#             ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),),
#             'NaN')
#
#         query = (
#             c.select(
#                 fn.min(c.catalogid).alias('catalogid'),
#                 fn.min(ls.ls_id).alias('ls_id'),
#                 fn.min(x.ero_detuid).alias('ero_detuid'),
#                 fn.min(c.ra).alias('ra'),
#                 fn.min(c.dec).alias('dec'),
#                 priority.alias("priority"),
#                 fn.min(value).alias('value'),
#                 fn.min(cadence).alias('cadence'),
#                 fn.min(instrument).alias('instrument'),
#                 fn.min(opt_prov).alias('optical_prov'),
#                 fn.min(magnitude_g).alias('g'),
#                 fn.min(magnitude_r).alias('r'),
#                 fn.min(magnitude_i).alias('i'),
#                 fn.min(magnitude_z).alias('z'),
#                 fn.min(magnitude_gaia_g).alias('gaia_g'),
#                 fn.min(magnitude_gaia_bp).alias('bp'),
#                 fn.min(magnitude_gaia_rp).alias('rp'),
#                 fn.min(g0_p).alias("ls8_mag_g"),   # extra
#                 fn.min(r0_p).alias("ls8_mag_r"),  # extra
#                 fn.min(z0_p).alias("ls8_mag_z"),  # extra
#                 fn.min(g0_e).alias("ls8_fibermag_g"),  # extra
#                 fn.min(r0_e).alias("ls8_fibermag_r"),  # extra
#                 fn.min(z0_e).alias("ls8_fibermag_z"),  # extra
#                 fn.min(ls.type).alias("ls8_type"),  # extra
#             )
#             .join(c2ls)
#             .where(
#                 c.version_id == version_id,
#                 c2ls.version_id == version_id,
#                 c2ls.best >> True,
#             )
#             .join(ls)
#             .join(x, on=(ls.ls_id == x.ls_id))
#             # start joining the spectroscopy
#             .switch(c)
#             .join(c2s16, JOIN.LEFT_OUTER)
#             .join(
#                 s16, JOIN.LEFT_OUTER,
#                 on=(
#                     (c2s16.target_id == s16.c.specobjid) &
#                     (c2s16.version_id == version_id)
#                 )
#             )
#             .switch(c)
#             .join(c2s2020, JOIN.LEFT_OUTER)
#             .join(
#                 s2020, JOIN.LEFT_OUTER,
#                 on=(
#                     (c2s2020.target_id == s2020.c.pk) &
#                     (c2s2020.version_id == version_id)
#                 )
#             )
#             .join(
#                 sV, JOIN.LEFT_OUTER,
#                 on=(
#                     fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
#                                 c.ra, c.dec,
#                                 match_radius_spectro)
#                 )
#             )
#             # .join(
#             #     sph, JOIN.LEFT_OUTER,
#             #     on=(
#             #         fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
#             #                     c.ra, c.dec,
#             #                     match_radius_spectro)
#             #     )
#             # )
#             # finished joining the spectroscopy
#             .where(
#                 (
#                     (x.ero_version == self.parameters['ero_version1']) &
#                     (x.xmatch_method == self.parameters['xmatch_method1']) &
#                     (x.xmatch_version == self.parameters['xmatch_version1'])
#                 ) |
#                 (
#                     (x.ero_version == self.parameters['ero_version2']) &
#                     (x.xmatch_method == self.parameters['xmatch_method2']) &
#                     (x.xmatch_version == self.parameters['xmatch_version2'])
#                 ),
#                 x.opt_cat == self.parameters['opt_cat'],
#                 x.xmatch_metric >= self.parameters['p_any_min'],
#                 (ls.fibertotflux_r.between(fibertotflux_r_min, fibertotflux_r_max) |
#                  ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max)),
#                 x.ero_det_like > self.parameters['det_like_min'],
#                 ls.maskbits.bin_and(2**2 + 2**3 + 2**4) == 0,  # avoid saturated sources
#                 ls.nobs_r > 0,                        # always require r-band coverage
#                 (ls.nobs_g > 0) | (ls.nobs_z > 0),    # plus at least one other optical band
#                 # gaia safety checks to avoid bad ls photometry
#                 ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters['gaia_g_mag_limit'])),
#                 ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters['gaia_rp_mag_limit'])),
#             )
#             .group_by(x.ls_id)   # avoid duplicates - we trust the legacy survey entries
#         )
#
#         if query_region:
#             query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
#                                                            query_region[0],
#                                                            query_region[1],
#                                                            query_region[2]))
#
#         return query
# #
# # END BhmSpidersAgnEfedsStragglersCarton
# # ##################################################################################


class BhmSpidersAgnGaiadr3Carton(BaseCarton):

    name = 'bhm_spiders_agn_gaiadr3'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        g3 = Gaia_DR3.alias()
        c2g3 = CatalogToGaia_DR3.alias()

        instrument = peewee.Value(self.instrument)

        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']
        gaia_g_max_for_cadence2 = self.parameters['gaia_g_max_for_cadence2']
        gaia_rp_max_for_cadence2 = self.parameters['gaia_rp_max_for_cadence2']

        # these control matching to spectroscopy
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # #########################################################################
        # prepare the spectroscopy catalogues

        # SDSS DR19p
        # downslect only 'good' spectra
        c2s19 = CatalogToSDSS_DR19p_Speclite.alias()
        ss19 = SDSS_DR19p_Speclite.alias()
        s19 = (
            ss19.select(
                ss19.pk.alias('s19_pk'),
            )
            .where(
                ss19.snmedian >= spec_sn_thresh,
                ss19.zwarning == 0,
                ss19.zerr <= spec_z_err_thresh,
                ss19.zerr > 0.0,
                ss19.scienceprimary > 0,
            )
            .alias('s19')
        )
        # #########################################################################

        # compute the abs(Galactic latitude):
        gal_lat = peewee.fn.abs(90.0 - peewee.fn.q3c_dist(north_gal_pole_ra,
                                                          north_gal_pole_dec,
                                                          c.ra, c.dec))

        # value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        value = peewee.Case(
            None,
            ((gal_lat > self.parameters['in_plane_lat_cut'],
              self.parameters.get('value', 1.0)),),
            0.0
        ).cast('float')

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +dpriority_match_flags if target is a secondary cross-match (match_flag > 1)
        # add +dpriority_det_like if target has a low value of ero_det_like
        # add +dpriority_has_spec if target has existing good SDSS spectroscopy

        priority_1 = peewee.Case(
            None,
            ((x.xmatch_flags > 1, 1), ),
            0)
        priority_2 = peewee.Case(
            None,
            ((x.ero_det_like < self.parameters['det_like_for_priority'], 1), ),
            0)
        priority_3 = peewee.Case(
            None,
            ((s19.c.s19_pk.is_null(False), 1), ),
            0)
        priority_4 = peewee.Case(
            None,
            ((gal_lat > self.parameters['in_plane_lat_cut'], 0),),
            1
        )

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec'] +
            priority_4 * self.parameters['dpriority_in_plane']
        )

        # choose cadence based on magnitude in Gaia G and RP-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((g3.phot_g_mean_mag < gaia_g_max_for_cadence1) |
                 (g3.phot_rp_mean_mag < gaia_rp_max_for_cadence1), cadence1),
                ((g3.phot_g_mean_mag < gaia_g_max_for_cadence2) |
                 (g3.phot_rp_mean_mag < gaia_rp_max_for_cadence2), cadence2),
                ((g3.phot_g_mean_mag >= gaia_g_max_for_cadence2) &
                 (g3.phot_rp_mean_mag >= gaia_rp_max_for_cadence2), cadence3),
            ),
            cadence4)

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

        g_blue = (g3.gaiamag + coeffs['g0_blue'] + coeffs['g1_blue'] * g3.bp_rp +
                  coeffs['g2_blue'] * g3.bp_rp * g3.bp_rp +
                  coeffs['g3_blue'] * g3.bp_rp * g3.bp_rp * g3.bp_rp)
        r_blue = (g3.gaiamag + coeffs['r0_blue'] + coeffs['r1_blue'] * g3.bp_rp +
                  coeffs['r2_blue'] * g3.bp_rp * g3.bp_rp +
                  coeffs['r3_blue'] * g3.bp_rp * g3.bp_rp * g3.bp_rp)
        i_blue = (g3.gaiamag + coeffs['i0_blue'] + coeffs['i1_blue'] * g3.bp_rp +
                  coeffs['i2_blue'] * g3.bp_rp * g3.bp_rp +
                  coeffs['i3_blue'] * g3.bp_rp * g3.bp_rp * g3.bp_rp)
        z_blue = (g3.gaiamag + coeffs['z0_blue'] + coeffs['z1_blue'] * g3.bp_rp +
                  coeffs['z2_blue'] * g3.bp_rp * g3.bp_rp +
                  coeffs['z3_blue'] * g3.bp_rp * g3.bp_rp * g3.bp_rp)

        g_red = (g3.gaiamag + coeffs['g0_red'] + coeffs['g1_red'] * g3.bp_rp +
                 coeffs['g2_red'] * g3.bp_rp * g3.bp_rp)
        r_red = (g3.gaiamag + coeffs['r0_red'] + coeffs['r1_red'] * g3.bp_rp +
                 coeffs['r2_red'] * g3.bp_rp * g3.bp_rp)
        i_red = (g3.gaiamag + coeffs['i0_red'] + coeffs['i1_red'] * g3.bp_rp +
                 coeffs['i2_red'] * g3.bp_rp * g3.bp_rp)
        z_red = (g3.gaiamag + coeffs['z0_red'] + coeffs['z1_red'] * g3.bp_rp +
                 coeffs['z2_red'] * g3.bp_rp * g3.bp_rp)

        # validity checks - set limits semi-manually
        bp_rp_min = 0.0
        bp_rp_max = 3.0
        valid = (g3.gaiamag.between(0.1, 29.9) &
                 g3.gaiabp.between(0.1, 29.9) &
                 g3.gaiarp.between(0.1, 29.9) &
                 g3.bp_rp.between(bp_rp_min, bp_rp_max))
        opt_prov = peewee.Case(None, ((valid, 'sdss_psfmag_from_gdr2'),), 'undefined')
        magnitude_g = peewee.Case(None,
                                  (
                                      (valid & (g3.bp_rp < 1.8), g_blue),
                                      (valid & (g3.bp_rp > 1.8), g_red),
                                  ), 'NaN')
        magnitude_r = peewee.Case(None,
                                  (
                                      (valid & (g3.bp_rp < 1.8), r_blue),
                                      (valid & (g3.bp_rp > 1.8), r_red),
                                  ), 'NaN')
        magnitude_i = peewee.Case(None,
                                  (
                                      (valid & (g3.bp_rp < 1.8), i_blue),
                                      (valid & (g3.bp_rp > 1.8), i_red),
                                  ), 'NaN')
        magnitude_z = peewee.Case(None,
                                  (
                                      (valid & (g3.bp_rp < 1.8), z_blue),
                                      (valid & (g3.bp_rp > 1.8), z_red),
                                  ), 'NaN')

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(x.ero_detuid).alias('ero_detuid'),  # extra
                fn.min(g3.source_id).alias('gaia_dr3_source_id'),  # extra
                fn.min(s19.c.s19_pk).alias('sdss_dr19p_speclite_pk'), # extra
                fn.min(c.ra).alias('ra'),  # extra
                fn.min(c.dec).alias('dec'),  # extra
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(opt_prov).alias('optical_prov'),
                fn.min(magnitude_g).alias('g'),
                fn.min(magnitude_r).alias('r'),
                fn.min(magnitude_i).alias('i'),
                fn.min(magnitude_z).alias('z'),
                fn.min(g3.phot_g_mean_mag).alias('gaia_g'),
                fn.min(g3.phot_bp_mean_mag).alias('bp'),
                fn.min(g3.phot_rp_mean_mag).alias('rp'),
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
            .join(
                s19, JOIN.LEFT_OUTER,
                on=(s19.c.s19_pk == c2s19.target_id))
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (x.xmatch_version == self.parameters['xmatch_version']),
                (x.opt_cat == self.parameters['opt_cat']),
                (x.xmatch_metric >= self.parameters['p_any_min']),
                (g3.phot_g_mean_mag > self.parameters['gaia_g_mag_limit']),
                (g3.phot_rp_mean_mag > self.parameters['gaia_rp_mag_limit']),
                (x.ero_det_like > self.parameters['det_like_min']),
            )
            .group_by(g3.source_id)   # avoid duplicates - we trust the gaia ids
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
#
# END BhmSpidersAgnGaiadr3Carton
# ##################################################################################


# we can get away with just inheriting the selection code from
# the gaia dr3 hemisphere match and adjusting the parameters only
class BhmSpidersAgnSepCarton(BhmSpidersAgnGaiadr3Carton):
    name = 'bhm_spiders_agn_sep'


# class BhmSpidersAgnPs1dr2Carton(BaseCarton):
#
#     name = 'bhm_spiders_agn_ps1dr2'
#     category = 'science'
#     mapper = 'BHM'
#     program = 'bhm_spiders'
#     tile = False
#     instrument = 'BOSS'
#
#     def build_query(self, version_id, query_region=None):
#
#         c = Catalog.alias()
#         x = EROSITASupersetAGN.alias()
#         ps = Panstarrs1.alias()
#         c2ps = CatalogToPanstarrs1.alias()
#         g3 = Gaia_DR3.alias()
#         c2g3 = CatalogToGaia_DR3.alias()
#
#         instrument = peewee.Value(self.instrument)
#
#         g_psf_flux_max = AB2Jy(self.parameters['g_psf_mag_min'])
#         r_psf_flux_max = AB2Jy(self.parameters['r_psf_mag_min'])
#         i_psf_flux_max = AB2Jy(self.parameters['i_psf_mag_min'])
#         z_psf_flux_max = AB2Jy(self.parameters['z_psf_mag_min'])
#         g_psf_flux_min = AB2Jy(self.parameters['g_psf_mag_max'])
#         r_psf_flux_min = AB2Jy(self.parameters['r_psf_mag_max'])
#         i_psf_flux_min = AB2Jy(self.parameters['i_psf_mag_max'])
#         z_psf_flux_min = AB2Jy(self.parameters['z_psf_mag_max'])
#         g_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['g_psf_mag_max_for_cadence1'])
#         r_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence1'])
#         i_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence1'])
#         g_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['g_psf_mag_max_for_cadence2'])
#         r_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence2'])
#         i_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence2'])
#
#         # value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
#
#         # these control matching to spectroscopy
#         spec_sn_thresh = self.parameters['spec_sn_thresh']
#         spec_z_err_thresh = self.parameters['spec_z_err_thresh']
#
#         # #########################################################################
#         # prepare the spectroscopy catalogues
#
#         # SDSS DR19p
#         # downslect only 'good' spectra
#         c2s19 = CatalogToSDSS_DR19p_Speclite.alias()
#         ss19 = SDSS_DR19p_Speclite.alias()
#         s19 = (
#             ss19.select(
#                 ss19.pk.alias('s19_pk'),
#             )
#             .where(
#                 ss19.snmedian >= spec_sn_thresh,
#                 ss19.zwarning == 0,
#                 ss19.zerr <= spec_z_err_thresh,
#                 ss19.zerr > 0.0,
#                 ss19.scienceprimary > 0,
#             )
#             .alias('s19')
#         )
#         # #########################################################################
#
#         # compute the abs(Galactic latitude):
#         gal_lat = peewee.fn.abs(90.0 - peewee.fn.q3c_dist(north_gal_pole_ra,
#                                                           north_gal_pole_dec,
#                                                           c.ra, c.dec))
#         value = peewee.Case(
#             None,
#             ((gal_lat > self.parameters['in_plane_lat_cut'],
#               self.parameters.get('value', 1.0)),),
#             0.0
#         ).cast('float')
#
#         # priority is determined by target properties
#         # start with a priority floor value (per carton)
#         # then increment if any conditions are met:
#         # add +dpriority_match_flags if target is a secondary cross-match (match_flag > 1)
#         # add +dpriority_det_like if target has a low value of ero_det_like
#         # add +dpriority_has_spec if target has existing good SDSS spectroscopy
#         # add +dpriority_in_plane if target lies at |b| < in_plane_lat_cut
#
#         priority_1 = peewee.Case(
#             None,
#             ((x.xmatch_flags > 1, 1), ),
#             0)
#         priority_2 = peewee.Case(
#             None,
#             ((x.ero_det_like < self.parameters['det_like_for_priority'], 1), ),
#             0)
#         priority_3 = peewee.Case(
#             None,
#             ((s19.c.s19_pk.is_null(False), 1), ),
#             0)
#         priority_4 = peewee.Case(
#             None,
#             ((gal_lat > self.parameters['in_plane_lat_cut'], 0),),
#             1
#         )
#
#         priority = fn.max(
#             self.parameters['priority_floor'] +
#             priority_1 * self.parameters['dpriority_match_flags'] +
#             priority_2 * self.parameters['dpriority_det_like'] +
#             priority_3 * self.parameters['dpriority_has_spec'] +
#             priority_4 * self.parameters['dpriority_in_plane']
#         )
#
#         # choose cadence based on psf_flux magnitude in panstarrs1 g,r,i-bands
#         cadence1 = self.parameters['cadence1']
#         cadence2 = self.parameters['cadence2']
#         cadence3 = self.parameters['cadence3']
#         cadence4 = 'unknown_cadence'
#         cadence = peewee.Case(
#             None,
#             (
#                 ((ps.g_stk_psf_flux > g_psf_flux_min_for_cadence1) |
#                  (ps.r_stk_psf_flux > r_psf_flux_min_for_cadence1) |
#                  (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1), cadence1),
#                 ((ps.g_stk_psf_flux > g_psf_flux_min_for_cadence2) |
#                  (ps.r_stk_psf_flux > r_psf_flux_min_for_cadence2) |
#                  (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence2), cadence2),
#                 ((ps.g_stk_psf_flux <= g_psf_flux_min_for_cadence2) &
#                  (ps.r_stk_psf_flux <= r_psf_flux_min_for_cadence2) &
#                  (ps.i_stk_psf_flux <= i_psf_flux_min_for_cadence2), cadence3),
#             ),
#             cadence4)
#
#         # compute transformed SDSS mags for pointlike and extended sources separately
#         # transform the panstarrs1-dr2 griz into sdss psfmag griz
#
#         # extract coeffs from fit logs via:
#         # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_agn_ps1dr2_pointlike/ps1dr2_stk_psf_to_sdss_psfmag_?_results.log  bhm_spiders_agn_ps1dr2_extended/ps1dr2_stk_psf_to_sdss_fiber2mag_?_results.log  # noqa
#
#         coeffs = {
#             "g2_p": 0.275586,
#             "g1_p": -0.178727,
#             "g0_p": 0.024900,
#             "i2_p": -0.051817,
#             "i1_p": 0.098077,
#             "i0_p": -0.028243,
#             "r2_p": -0.031567,
#             "r1_p": 0.056499,
#             "r0_p": -0.013487,
#             "z2_p": -0.290196,
#             "z1_p": 0.156009,
#             "z0_p": -0.079393,
#             "g2_e": 0.084856,
#             "g1_e": -0.076550,
#             "g0_e": 0.841168,
#             "i2_e": 0.048106,
#             "i1_e": 0.025289,
#             "i0_e": 0.652371,
#             "r2_e": 0.066827,
#             "r1_e": -0.118807,
#             "r0_e": 0.752550,
#             "z2_e": 0.558727,
#             "z1_e": -0.006461,
#             "z0_e": 0.512403,
#         }
#
#         Jy_min = AB2Jy(30.00)
#
#         # pointlike and extended - both start from ps1dr2 stk psf fluxes
#         g0 = (8.9 - 2.5 * peewee.fn.log(peewee.fn.greatest(Jy_min, ps.g_stk_psf_flux)))
#         r0 = (8.9 - 2.5 * peewee.fn.log(peewee.fn.greatest(Jy_min, ps.r_stk_psf_flux)))
#         i0 = (8.9 - 2.5 * peewee.fn.log(peewee.fn.greatest(Jy_min, ps.i_stk_psf_flux)))
#         z0 = (8.9 - 2.5 * peewee.fn.log(peewee.fn.greatest(Jy_min, ps.z_stk_psf_flux)))
#         g_r = g0 - r0
#         r_i = r0 - i0
#         i_z = i0 - z0
#
#         # use different transform coeffs for pointlike and extended sources
#         g_p = (g0 + coeffs['g0_p'] + coeffs['g1_p'] * g_r + coeffs['g2_p'] * g_r * g_r)
#         r_p = (r0 + coeffs['r0_p'] + coeffs['r1_p'] * g_r + coeffs['r2_p'] * g_r * g_r)
#         i_p = (i0 + coeffs['i0_p'] + coeffs['i1_p'] * r_i + coeffs['i2_p'] * r_i * r_i)
#         z_p = (z0 + coeffs['z0_p'] + coeffs['z1_p'] * i_z + coeffs['z2_p'] * i_z * i_z)
#
#         g_e = (g0 + coeffs['g0_e'] + coeffs['g1_e'] * g_r + coeffs['g2_e'] * g_r * g_r)
#         r_e = (r0 + coeffs['r0_e'] + coeffs['r1_e'] * g_r + coeffs['r2_e'] * g_r * g_r)
#         i_e = (i0 + coeffs['i0_e'] + coeffs['i1_e'] * r_i + coeffs['i2_e'] * r_i * r_i)
#         z_e = (z0 + coeffs['z0_e'] + coeffs['z1_e'] * i_z + coeffs['z2_e'] * i_z * i_z)
#
#         # validity checks - set limits semi-manually
#         g_r_p_min = -0.5
#         g_r_p_max = 1.4
#         r_i_p_min = -0.5
#         r_i_p_max = 2.0
#         i_z_p_min = -0.5
#         i_z_p_max = 1.0
#         g_r_e_min = -0.5
#         g_r_e_max = 1.5
#         r_i_e_min = -0.2
#         r_i_e_max = 1.3
#         i_z_e_min = -0.3
#         i_z_e_max = 0.8
#
#         valid_p = (g0.between(0.1, 29.9) &
#                    r0.between(0.1, 29.9) &
#                    z0.between(0.1, 29.9) &
#                    g_r.between(g_r_p_min, g_r_p_max) &
#                    r_i.between(r_i_p_min, r_i_p_max) &
#                    i_z.between(i_z_p_min, i_z_p_max))
#
#         valid_e = (g0.between(0.1, 29.9) &
#                    r0.between(0.1, 29.9) &
#                    z0.between(0.1, 29.9) &
#                    g_r.between(g_r_e_min, g_r_e_max) &
#                    r_i.between(r_i_e_min, r_i_e_max) &
#                    i_z.between(i_z_e_min, i_z_e_max))
#
#         # We want to switch between psfmags and fibertotmags depending on
#         # ps.flags EXT+EXT_ALT (i.e. extended sources)
#         ext_flags = 8388608 + 16777216
#         good_stack_flag = 134217728
#
#         opt_prov = peewee.Case(
#             None,
#             (
#                 ((ps.flags.bin_and(ext_flags) == 0) & valid_p, 'sdss_psfmag_from_ps1dr2'),
#                 ((ps.flags.bin_and(ext_flags) > 0) & valid_e, 'sdss_fiber2mag_from_ps1dr2'),
#             ),
#             'undefined')
#         magnitude_g = peewee.Case(
#             None,
#             (
#                 ((ps.flags.bin_and(ext_flags) == 0) & valid_p, g_p.cast('float')),
#                 ((ps.flags.bin_and(ext_flags) > 0) & valid_e, g_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_r = peewee.Case(
#             None,
#             (
#                 ((ps.flags.bin_and(ext_flags) == 0) & valid_p, r_p.cast('float')),
#                 ((ps.flags.bin_and(ext_flags) > 0) & valid_e, r_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_i = peewee.Case(
#             None,
#             (
#                 ((ps.flags.bin_and(ext_flags) == 0) & valid_p, i_p.cast('float')),
#                 ((ps.flags.bin_and(ext_flags) > 0) & valid_e, i_e.cast('float')),
#             ),
#             'NaN')
#         magnitude_z = peewee.Case(
#             None,
#             (
#                 ((ps.flags.bin_and(ext_flags) == 0) & valid_p, z_p.cast('float')),
#                 ((ps.flags.bin_and(ext_flags) > 0) & valid_e, z_e.cast('float')),
#             ),
#             'NaN')
#
#         ##############################
#
#         query = (
#             c.select(
#                 fn.min(c.catalogid).alias('catalogid'),
#                 fn.min(ps.catid_objid).alias('ps1_catid_objid'),  # extra
#                 fn.min(g3.source_id).alias('gaia_dr3_source_id'),  # extra
#                 fn.min(x.ero_detuid).alias('ero_detuid'),  # extra
#                 fn.min(s19.c.s19_pk).alias('sdss_dr19p_speclite_pk'),  # extra
#                 fn.min(c.ra).alias('ra'),  # extra
#                 fn.min(c.dec).alias('dec'),  # extra
#                 priority.alias("priority"),
#                 fn.min(value).alias('value'),
#                 fn.min(cadence).alias('cadence'),
#                 fn.min(instrument).alias('instrument'),
#                 fn.min(magnitude_g).alias('g'),
#                 fn.min(magnitude_r).alias('r'),
#                 fn.min(magnitude_i).alias('i'),
#                 fn.min(magnitude_z).alias('z'),
#                 fn.min(g3.phot_g_mean_mag).alias('gaia_g'),
#                 fn.min(g3.phot_bp_mean_mag).alias('bp'),
#                 fn.min(g3.phot_rp_mean_mag).alias('rp'),
#                 fn.min(opt_prov).alias('optical_prov'),
#                 fn.min(g0).alias("ps1dr2_stk_psf_mag_g"),   # extra
#                 fn.min(r0).alias("ps1dr2_stk_psf_mag_r"),   # extra
#                 fn.min(i0).alias("ps1dr2_stk_psf_mag_i"),   # extra
#                 fn.min(z0).alias("ps1dr2_stk_psf_mag_z"),   # extra
#                 fn.min(ps.flags).alias("ps1dr2_flags"),  # extra
#             )
#             .join(c2ps)
#             .where(
#                 c.version_id == version_id,
#                 c2ps.version_id == version_id,
#                 fn.coalesce(c2s19.version_id, version_id) == version_id,
#                 fn.coalesce(c2g3.version_id, version_id) == version_id,
#                 c2ps.best >> True,
#             )
#             .join(ps)
#             .join(x, on=(ps.catid_objid == x.ps1_dr2_id))
#             .switch(c)
#             .join(c2g3, JOIN.LEFT_OUTER)  # do not reject c2g3.best=False entries
#             .join(g3, JOIN.LEFT_OUTER)
#             # start joining the spectroscopy
#             .switch(c)
#             .join(c2s19, JOIN.LEFT_OUTER)
#             .join(s19, JOIN.LEFT_OUTER,
#                   on=(s19.c.s19_pk == c2s19.target_id))
#             # finished joining the spectroscopy
#             .where(
#                 (x.ero_version == self.parameters['ero_version']),
#                 (x.xmatch_method == self.parameters['xmatch_method']),
#                 (x.xmatch_version == self.parameters['xmatch_version']),
#                 (x.opt_cat == self.parameters['opt_cat']),
#                 (x.xmatch_metric >= self.parameters['p_any_min']),
#                 (x.ero_det_like > self.parameters['det_like_min']),
#                 (ps.g_stk_psf_flux < g_psf_flux_max),
#                 (ps.r_stk_psf_flux < r_psf_flux_max),
#                 (ps.i_stk_psf_flux < i_psf_flux_max),
#                 (ps.z_stk_psf_flux < z_psf_flux_max),
#                 (ps.g_stk_psf_flux != 'NaN'),   # these NaNs are NOT represented as nulls
#                 (ps.r_stk_psf_flux != 'NaN'),   # these NaNs are NOT represented as nulls
#                 (ps.i_stk_psf_flux != 'NaN'),   # these NaNs are NOT represented as nulls
#                 (
#                     (ps.g_stk_psf_flux > g_psf_flux_min) |
#                     (ps.r_stk_psf_flux > r_psf_flux_min) |
#                     (ps.i_stk_psf_flux > i_psf_flux_min) |
#                     (ps.z_stk_psf_flux > z_psf_flux_min)
#                 ),
#                 (ps.flags.bin_and(good_stack_flag) > 0),
#             )
#             .group_by(ps.catid_objid)   # avoid duplicates - we trust the ps1 ids
#             .having(
#                 # each and every match to gaia must satisfy the bright star rejection criteria
#                 fn.sum(
#                     peewee.Case(
#                         None,
#                         (
#                             (g3.phot_g_mean_mag < self.parameters['gaia_g_mag_limit'], 1),
#                             (g3.phot_rp_mean_mag < self.parameters['gaia_rp_mag_limit'], 1),
#                         ),
#                         0)
#                 ) == 0
#             )
#         )
#
#         if query_region:
#             query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
#                                                            query_region[0],
#                                                            query_region[1],
#                                                            query_region[2]))
#
#         return query
# #
# # END BhmSpidersAgnPs1dr2Carton
# # ##################################################################################
#
