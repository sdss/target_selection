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

# imports of existing spectro catalogues
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToSDSS_DR16_SpecObj,
    SDSS_DR16_SpecObj,
    CatalogToBHM_eFEDS_Veto,
    BHM_eFEDS_Veto,
    SDSSV_BOSS_SPALL,
    SDSSV_BOSS_Conflist,
    SDSSV_Plateholes,
    SDSSV_Plateholes_Meta,
)

# additional imports required by bhm_spiders_agn_lsdr8
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToLegacy_Survey_DR8,
    Legacy_Survey_DR8,
)

# additional imports required by bhm_spiders_agn_gaiadr2
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToTIC_v8,
    TIC_v8,
)

# additional imports required by bhm_spiders_agn_ps1dr2
from sdssdb.peewee.sdss5db.catalogdb import (
    Panstarrs1,
    CatalogToPanstarrs1,
)

# additional imports required by bhm_spiders_agn_skymapperdr2
from sdssdb.peewee.sdss5db.catalogdb import (
    SkyMapper_DR2,
    CatalogToSkyMapper_DR2,
)


# additional imports required by bhm_spiders_agn_supercosmos
from sdssdb.peewee.sdss5db.catalogdb import (
    SuperCosmos,
    CatalogToSuperCosmos,
    CatalogToCatWISE2020,
)


from target_selection.mag_flux import AB2nMgy, AB2Jy


# ############################################
# ############################################
# ############################################
# ############################################
# This provides the following BHM SPIDERS AGN cartons in v0.5:
#  *  bhm_spiders_agn_lsdr8
#  *  bhm_spiders_agn_efeds_stragglers
#  *  bhm_spiders_agn_gaiadr2
#  *  bhm_spiders_agn_sep
#  *  bhm_spiders_agn_ps1dr2
#  *  bhm_spiders_agn_skymapperdr2
#     bhm_spiders_agn_supercosmos
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


class BhmSpidersAgnLsdr8Carton(BaseCarton):

    name = 'bhm_spiders_agn_lsdr8'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()

        instrument = peewee.Value(self.instrument)

        fibertotflux_r_max = AB2nMgy(self.parameters['fibertotmag_r_min'])
        fibertotflux_r_min = AB2nMgy(self.parameters['fibertotmag_r_max'])
        fibertotflux_z_max = AB2nMgy(self.parameters['fibertotmag_z_min'])
        fibertotflux_z_min = AB2nMgy(self.parameters['fibertotmag_z_max'])

        fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence1'])
        fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence2'])
        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']

        flux30 = AB2nMgy(30.00)
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # #########################################################################
        # prepare the spectroscopy catalogues

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

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
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                (sph.c.pkey.is_null(False), 1),
            ),
            0)

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec']
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

        # We want to switch between psfmags and fibertotmags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            ls.type,
            (('PSF', 'ls_psfmag'),),
            'ls_fibertotmag')

        magnitude_g = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_g))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_g))).cast('float'))

        magnitude_r = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_r))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_r))).cast('float'))

        magnitude_z = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_z))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_z))).cast('float'))

        magnitude_i = peewee.Case(
            ls.type,
            (('PSF',
              (22.5 - 2.5 * fn.log10(
                  fn.greatest(flux30, 0.5 * (ls.flux_r + ls.flux_z)))).cast('float')),),
            (22.5 - 2.5 * fn.log10(
                fn.greatest(flux30, 0.5 * (ls.fibertotflux_r + ls.fibertotflux_z)))).cast('float'))

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(ls.ls_id).alias('ls_id'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(magnitude_g).alias('g'),
                fn.min(magnitude_r).alias('r'),
                fn.min(magnitude_i).alias('i'),
                fn.min(magnitude_z).alias('z'),
                fn.min(opt_prov).alias('opt_prov'),
            )
            .join(c2ls)
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True,
            )
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
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
# END BhmSpidersAgnLsdr8Carton
# ##################################################################################


class BhmSpidersAgnEfedsStragglersCarton(BaseCarton):

    name = 'bhm_spiders_agn_efeds_stragglers'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'
    cadence = None

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()

        instrument = peewee.Value(self.instrument)

        fibertotflux_r_max = AB2nMgy(self.parameters['fibertotmag_r_min'])
        fibertotflux_r_min = AB2nMgy(self.parameters['fibertotmag_r_max'])
        fibertotflux_z_max = AB2nMgy(self.parameters['fibertotmag_z_min'])
        fibertotflux_z_min = AB2nMgy(self.parameters['fibertotmag_z_max'])

        fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence1'])
        fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence2'])
        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']

        flux30 = AB2nMgy(30.00)
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        #
        ##########################################################################
        # prepare the spectroscopy catalogues

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # All eFEDS plates have been observed so ignore plateholes now
        # # SDSS-V plateholes - only consider plateholes that
        # # were drilled+shipped but that were not yet observed
        # ssph = SDSSV_Plateholes.alias()
        # ssphm = SDSSV_Plateholes_Meta.alias()
        # ssconf = SDSSV_BOSS_Conflist.alias()
        # sph = (
        #     ssph.select(
        #         ssph.pkey.alias('pkey'),
        #         ssph.target_ra.alias('target_ra'),
        #         ssph.target_dec.alias('target_dec'),
        #     )
        #     .join(
        #         ssphm,
        #         on=(ssph.yanny_uid == ssphm.yanny_uid)
        #     )
        #     .join(
        #         ssconf, JOIN.LEFT_OUTER,
        #         on=(ssphm.plateid == ssconf.plate)
        #     )
        #     .where(
        #         (ssph.holetype == 'BOSS_SHARED'),
        #         (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
        #         ssphm.isvalid > 0,
        #         ssconf.plate.is_null(),
        #     )
        #     .alias('sph')
        # )

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target is from the secondary eFEDS catalogue
        # add +80 if target has existing good SDSS spectroscopy

        priority_1 = peewee.Case(
            None,
            ((x.xmatch_flags > 1, 1), ),
            0)
        priority_2 = peewee.Case(
            None,
            ((x.ero_det_like < self.parameters['det_like_for_priority'], 1), ),
            0)
        # add a step down in priority for anything only selected by the secondary xmatch_version
        priority_3 = peewee.Case(
            None,
            ((x.ero_version == self.parameters['ero_version2'], 1), ),
            0)

        # no need for the step below because we just reject everything with a spectrum.
        priority_4 = peewee.Case(
            None,
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                # (sph.c.pkey.is_null(False), 1),
            ),
            0)

        priority = (
            self.parameters['priority_floor'] +
            fn.min(priority_1) * self.parameters['dpriority_match_flags'] +
            fn.min(priority_2) * self.parameters['dpriority_det_like'] +
            fn.min(priority_3) * self.parameters['dpriority_ero_version'] +
            fn.max(priority_4) * self.parameters['dpriority_has_spec']
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

        # We want to switch between psfmags and fibertotmags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            ls.type,
            (('PSF', 'ls_psfmag'),),
            'ls_fibertotmag')

        magnitude_g = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_g))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_g))).cast('float'))

        magnitude_r = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_r))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_r))).cast('float'))

        magnitude_z = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_z))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_z))).cast('float'))

        magnitude_i = peewee.Case(
            ls.type,
            (('PSF',
              (22.5 - 2.5 * fn.log10(
                  fn.greatest(flux30, 0.5 * (ls.flux_r + ls.flux_z)))).cast('float')),),
            (22.5 - 2.5 * fn.log10(
                fn.greatest(flux30, 0.5 * (ls.fibertotflux_r + ls.fibertotflux_z)))).cast('float'))

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(ls.ls_id).alias('ls_id'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(magnitude_g).alias('g'),
                fn.min(magnitude_r).alias('r'),
                fn.min(magnitude_i).alias('i'),
                fn.min(magnitude_z).alias('z'),
                fn.min(ls.gaia_phot_g_mean_mag).alias('gaia_g'),   # extra
                fn.min(ls.gaia_phot_g_mean_mag).alias('gaia_rp'),  # extra
                fn.min(opt_prov).alias('opt_prov'),
            )
            .join(c2ls)
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True,
            )
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # .join(
            #     sph, JOIN.LEFT_OUTER,
            #     on=(
            #         fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
            #                     c.ra, c.dec,
            #                     match_radius_spectro)
            #     )
            # )
            # finished joining the spectroscopy
            .where(
                (
                    (x.ero_version == self.parameters['ero_version1']) &
                    (x.xmatch_method == self.parameters['xmatch_method1']) &
                    (x.xmatch_version == self.parameters['xmatch_version1'])
                ) |
                (
                    (x.ero_version == self.parameters['ero_version2']) &
                    (x.xmatch_method == self.parameters['xmatch_method2']) &
                    (x.xmatch_version == self.parameters['xmatch_version2'])
                ),
                x.opt_cat == self.parameters['opt_cat'],
                x.xmatch_metric >= self.parameters['p_any_min'],
                (ls.fibertotflux_r.between(fibertotflux_r_min, fibertotflux_r_max) |
                 ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max)),
                x.ero_det_like > self.parameters['det_like_min'],
                ls.maskbits.bin_and(2**2 + 2**3 + 2**4) == 0,  # avoid saturated sources
                ls.nobs_r > 0,                        # always require r-band coverage
                (ls.nobs_g > 0) | (ls.nobs_z > 0),    # plus at least one other optical band
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters['gaia_g_mag_limit'])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters['gaia_rp_mag_limit'])),
            )
            .group_by(x.ls_id)   # avoid duplicates - we trust the legacy survey entries
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
#
# END BhmSpidersAgnEfedsStragglersCarton
# ##################################################################################


class BhmSpidersAgnGaiadr2Carton(BaseCarton):

    name = 'bhm_spiders_agn_gaiadr2'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        instrument = peewee.Value(self.instrument)

        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']
        gaia_g_max_for_cadence2 = self.parameters['gaia_g_max_for_cadence2']
        gaia_rp_max_for_cadence2 = self.parameters['gaia_rp_max_for_cadence2']

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # #########################################################################
        # prepare the spectroscopy catalogues

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )
        # #########################################################################

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

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
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                (sph.c.pkey.is_null(False), 1),
            ),
            0)

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec']
        )

        # choose cadence based on magnitude in Gaia G and RP-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((tic.gaiamag < gaia_g_max_for_cadence1) |
                 (tic.gaiarp < gaia_rp_max_for_cadence1), cadence1),
                ((tic.gaiamag < gaia_g_max_for_cadence2) |
                 (tic.gaiarp < gaia_rp_max_for_cadence2), cadence2),
                ((tic.gaiamag >= gaia_g_max_for_cadence2) &
                 (tic.gaiarp >= gaia_rp_max_for_cadence2), cadence3),
            ),
            cadence4)

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(tic.gaia_int).alias('gaia_source'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(tic.gaiamag).alias('gaia_g'),
                fn.min(tic.gaiabp).alias('bp'),
                fn.min(tic.gaiarp).alias('rp'),
            )
            .join(c2tic)
            .where(
                c.version_id == version_id,
                c2tic.version_id == version_id,
                c2tic.best >> True,
            )
            .join(tic)
            .join(x, on=(tic.gaia_int == x.gaia_dr2_id))
            .switch(c)
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (x.xmatch_version == self.parameters['xmatch_version']),
                (x.opt_cat == self.parameters['opt_cat']),
                (x.xmatch_metric >= self.parameters['p_any_min']),
                (tic.gaiamag > self.parameters['gaia_g_mag_limit']),
                (tic.gaiarp > self.parameters['gaia_rp_mag_limit']),
                (x.ero_det_like > self.parameters['det_like_min']),
            )
            .group_by(tic.gaia_int)   # avoid duplicates - we trust the gaia ids
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
#
# END BhmSpidersAgnGaiadr2Carton
# ##################################################################################


class BhmSpidersAgnSepCarton(BaseCarton):

    name = 'bhm_spiders_agn_sep'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        instrument = peewee.Value(self.instrument)

        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']
        gaia_g_max_for_cadence2 = self.parameters['gaia_g_max_for_cadence2']
        gaia_rp_max_for_cadence2 = self.parameters['gaia_rp_max_for_cadence2']

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy - N/A here

        priority_1 = peewee.Case(
            None,
            ((x.xmatch_flags > 1, 1), ),
            0)
        priority_2 = peewee.Case(
            None,
            ((x.ero_det_like < self.parameters['det_like_for_priority'], 1), ),
            0)
        # add a step down in priority for anything selected by the secondary xmatch_version
        priority_3 = peewee.Case(
            None,
            ((x.xmatch_method == self.parameters['xmatch_version2'], 1), ),
            0)

        # choose the maximum priority option for all combinations of this target
        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_match_method']
        )

        # choose cadence based on magnitude in Gaia G and RP-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((tic.gaiamag < gaia_g_max_for_cadence1) |
                 (tic.gaiarp < gaia_rp_max_for_cadence1), cadence1),
                ((tic.gaiamag < gaia_g_max_for_cadence2) |
                 (tic.gaiarp < gaia_rp_max_for_cadence2), cadence2),
                ((tic.gaiamag >= gaia_g_max_for_cadence2) &
                 (tic.gaiarp >= gaia_rp_max_for_cadence2), cadence3),
            ),
            cadence4)

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(tic.gaia_int).alias('gaia_source'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(tic.gaiamag).alias('gaia_g'),
                fn.min(tic.gaiabp).alias('bp'),
                fn.min(tic.gaiarp).alias('rp'),
            )
            .join(c2tic)
            .where(
                c.version_id == version_id,
                c2tic.version_id == version_id,
                c2tic.best >> True,
            )
            .join(tic)
            .join(x, on=(tic.gaia_int == x.gaia_dr2_id))
            .switch(c)
            .where(
                x.ero_version == self.parameters['ero_version'],
                x.xmatch_method == self.parameters['xmatch_method'],
                (
                    (x.xmatch_version == self.parameters['xmatch_version1']) |
                    (x.xmatch_version == self.parameters['xmatch_version2'])
                ),
                x.opt_cat == self.parameters['opt_cat'],
                x.xmatch_metric >= self.parameters['p_any_min'],
                tic.gaiamag > self.parameters['gaia_g_mag_limit'],
                tic.gaiarp > self.parameters['gaia_rp_mag_limit'],
                x.ero_det_like > self.parameters['det_like_min'],
            )
            .group_by(tic.gaia_int)   # avoid duplicates - we trust the gaia ids
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query

#
# END BhmSpidersAgnSepCarton
# ##################################################################################


class BhmSpidersAgnPs1dr2Carton(BaseCarton):

    name = 'bhm_spiders_agn_ps1dr2'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        ps = Panstarrs1.alias()
        c2ps = CatalogToPanstarrs1.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        instrument = peewee.Value(self.instrument)

        g_psf_flux_max = AB2Jy(self.parameters['g_psf_mag_min'])
        r_psf_flux_max = AB2Jy(self.parameters['r_psf_mag_min'])
        i_psf_flux_max = AB2Jy(self.parameters['i_psf_mag_min'])
        z_psf_flux_max = AB2Jy(self.parameters['z_psf_mag_min'])
        g_psf_flux_min = AB2Jy(self.parameters['g_psf_mag_max'])
        r_psf_flux_min = AB2Jy(self.parameters['r_psf_mag_max'])
        i_psf_flux_min = AB2Jy(self.parameters['i_psf_mag_max'])
        z_psf_flux_min = AB2Jy(self.parameters['z_psf_mag_max'])
        g_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['g_psf_mag_max_for_cadence1'])
        r_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence1'])
        i_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence1'])
        g_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['g_psf_mag_max_for_cadence2'])
        r_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence2'])
        i_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence2'])

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # #########################################################################
        # prepare the spectroscopy catalogues

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )
        # #########################################################################

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

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
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                (sph.c.pkey.is_null(False), 1),
            ),
            0)

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec']
        )

        # choose cadence based on psf_flux magnitude in panstarrs1 g,r,i-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((ps.g_stk_psf_flux > g_psf_flux_min_for_cadence1) |
                 (ps.r_stk_psf_flux > r_psf_flux_min_for_cadence1) |
                 (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1), cadence1),
                ((ps.g_stk_psf_flux > g_psf_flux_min_for_cadence2) |
                 (ps.r_stk_psf_flux > r_psf_flux_min_for_cadence2) |
                 (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence2), cadence2),
                ((ps.g_stk_psf_flux <= g_psf_flux_min_for_cadence2) &
                 (ps.r_stk_psf_flux <= r_psf_flux_min_for_cadence2) &
                 (ps.i_stk_psf_flux <= i_psf_flux_min_for_cadence2), cadence3),
            ),
            cadence4)
        # We want to switch between psfmags and fibertotmags depending on
        # ps.flags EXT+EXT_ALT (i.e. extended sources)
        # For non-extended targets, we use psfmags, but for extended sources use apermag
        flux30 = AB2Jy(30.00)
        ext_flags = 8388608 + 16777216
        good_stack_flag = 134217728
        opt_prov = peewee.Case(
            ps.flags.bin_and(ext_flags),
            ((0, 'ps_psfmag'),),
            'ps_apermag')

        magnitude_g = peewee.Case(
            ps.flags.bin_and(ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.g_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.g_stk_aper_flux))).cast('float'))

        magnitude_r = peewee.Case(
            ps.flags.bin_and(ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.r_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.r_stk_aper_flux))).cast('float'))

        magnitude_i = peewee.Case(
            ps.flags.bin_and(ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.i_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.i_stk_aper_flux))).cast('float'))

        magnitude_z = peewee.Case(
            ps.flags.bin_and(ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.z_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.z_stk_aper_flux))).cast('float'))

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(ps.catid_objid).alias('ps1_catid_objid'),
                fn.min(tic.gaia_int).alias('gaia_source'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(magnitude_g).alias('g'),
                fn.min(magnitude_r).alias('r'),
                fn.min(magnitude_i).alias('i'),
                fn.min(magnitude_z).alias('z'),
                fn.min(tic.gaiamag).alias('gaia_g'),
                fn.min(tic.gaiabp).alias('bp'),
                fn.min(tic.gaiarp).alias('rp'),
                fn.min(opt_prov).alias('opt_prov'),
            )
            .join(c2ps)
            .where(
                c.version_id == version_id,
                c2ps.version_id == version_id,
                c2ps.best >> True,
            )
            .join(ps)
            .join(x, on=(ps.catid_objid == x.ps1_dr2_id))
            .switch(c)
            .join(
                c2tic, JOIN.LEFT_OUTER,
                on=(
                    (c.catalogid == c2tic.catalogid) &
                    (c2tic.version_id == version_id)
                    # explicitly do not want to take only c2tic.best entries
                )
            )
            .join(tic, JOIN.LEFT_OUTER)
            .switch(c)
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (x.xmatch_version == self.parameters['xmatch_version']),
                (x.opt_cat == self.parameters['opt_cat']),
                (x.xmatch_metric >= self.parameters['p_any_min']),
                (x.ero_det_like > self.parameters['det_like_min']),
                (ps.g_stk_psf_flux < g_psf_flux_max),
                (ps.r_stk_psf_flux < r_psf_flux_max),
                (ps.i_stk_psf_flux < i_psf_flux_max),
                (ps.z_stk_psf_flux < z_psf_flux_max),
                (ps.g_stk_psf_flux != 'NaN'),   # these NaNs are NOT represented as nulls
                (ps.r_stk_psf_flux != 'NaN'),   # these NaNs are NOT represented as nulls
                (ps.i_stk_psf_flux != 'NaN'),   # these NaNs are NOT represented as nulls
                (
                    (ps.g_stk_psf_flux > g_psf_flux_min) |
                    (ps.r_stk_psf_flux > r_psf_flux_min) |
                    (ps.i_stk_psf_flux > i_psf_flux_min) |
                    (ps.z_stk_psf_flux > z_psf_flux_min)
                ),
                (ps.flags.bin_and(good_stack_flag) > 0),
            )
            .group_by(ps.catid_objid)   # avoid duplicates - we trust the ps1 ids
            .having(
                # each and every match to the tic must satisfy the bright star rejection criteria
                fn.sum(
                    peewee.Case(
                        None,
                        (
                            (tic.gaiamag < self.parameters['gaia_g_mag_limit'], 1),
                            (tic.gaiarp < self.parameters['gaia_rp_mag_limit'], 1),
                            (tic.tmag < self.parameters['tic_t_mag_limit'], 1),
                        ),
                        0)
                ) == 0
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
#
# END BhmSpidersAgnPs1dr2Carton
# ##################################################################################


class BhmSpidersAgnSkyMapperDr2Carton(BaseCarton):

    name = 'bhm_spiders_agn_skymapperdr2'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        sm = SkyMapper_DR2.alias()
        c2sm = CatalogToSkyMapper_DR2.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        instrument = peewee.Value(self.instrument)

        g_psf_mag_min = self.parameters['g_psf_mag_min']
        r_psf_mag_min = self.parameters['r_psf_mag_min']
        i_psf_mag_min = self.parameters['i_psf_mag_min']
        z_psf_mag_min = self.parameters['z_psf_mag_min']
        g_psf_mag_max = self.parameters['g_psf_mag_max']
        r_psf_mag_max = self.parameters['r_psf_mag_max']
        i_psf_mag_max = self.parameters['i_psf_mag_max']
        z_psf_mag_max = self.parameters['z_psf_mag_max']
        g_psf_mag_max_for_cadence1 = self.parameters['g_psf_mag_max_for_cadence1']
        r_psf_mag_max_for_cadence1 = self.parameters['r_psf_mag_max_for_cadence1']
        i_psf_mag_max_for_cadence1 = self.parameters['i_psf_mag_max_for_cadence1']
        g_psf_mag_max_for_cadence2 = self.parameters['g_psf_mag_max_for_cadence2']
        r_psf_mag_max_for_cadence2 = self.parameters['r_psf_mag_max_for_cadence2']
        i_psf_mag_max_for_cadence2 = self.parameters['i_psf_mag_max_for_cadence2']

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # #########################################################################
        # prepare the spectroscopy catalogues

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )
        # #########################################################################

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

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
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                (sph.c.pkey.is_null(False), 1),
            ),
            0)

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec']
        )

        # choose cadence based on psf_flux magnitude in skymapper g,r,i-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((sm.g_psf < g_psf_mag_max_for_cadence1) |
                 (sm.r_psf < r_psf_mag_max_for_cadence1) |
                 (sm.i_psf < i_psf_mag_max_for_cadence1), cadence1),
                ((sm.g_psf < g_psf_mag_max_for_cadence2) |
                 (sm.r_psf < r_psf_mag_max_for_cadence2) |
                 (sm.i_psf < i_psf_mag_max_for_cadence2), cadence2),
                ((sm.g_psf >= g_psf_mag_max_for_cadence2) &
                 (sm.r_psf >= r_psf_mag_max_for_cadence2) &
                 (sm.i_psf >= i_psf_mag_max_for_cadence2), cadence3),
            ),
            cadence4)

        # We want to only use psfmags
        opt_prov = peewee.Value('sm_psfmag')

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(sm.object_id).alias('sm2_object_id'),
                fn.min(tic.gaia_int).alias('gaia_source'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(sm.g_psf).alias('g'),
                fn.min(sm.r_psf).alias('r'),
                fn.min(sm.i_psf).alias('i'),
                fn.min(sm.z_psf).alias('z'),
                fn.min(tic.gaiamag).alias('gaia_g'),
                fn.min(tic.gaiabp).alias('bp'),
                fn.min(tic.gaiarp).alias('rp'),
                fn.min(opt_prov).alias('opt_prov'),
            )
            .join(c2sm)
            .where(
                c.version_id == version_id,
                c2sm.version_id == version_id,
                c2sm.best >> True,
            )
            .join(sm)
            .join(x, on=(sm.object_id == x.skymapper_dr2_id))
            .switch(c)
            .join(
                c2tic, JOIN.LEFT_OUTER,
                on=(
                    (c.catalogid == c2tic.catalogid) &
                    (c2tic.version_id == version_id)
                    # & (c2tic.best >> True)   # No, we want all possible matches
                )
            )
            .join(tic, JOIN.LEFT_OUTER)
            .switch(c)
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (x.xmatch_version == self.parameters['xmatch_version']),
                (x.opt_cat == self.parameters['opt_cat']),
                (x.xmatch_metric >= self.parameters['p_any_min']),
                (x.ero_det_like > self.parameters['det_like_min']),
                (sm.g_psf > g_psf_mag_min),
                (sm.r_psf > r_psf_mag_min),
                (sm.i_psf > i_psf_mag_min),
                (sm.z_psf > z_psf_mag_min),
                (sm.g_psf.is_null(False) | sm.r_psf.is_null(False) | sm.i_psf.is_null(False)),
                (
                    (sm.g_psf < g_psf_mag_max) |
                    (sm.r_psf < r_psf_mag_max) |
                    (sm.i_psf < i_psf_mag_max) |
                    (sm.z_psf < z_psf_mag_max)
                ),
                # see description of skymapper flags here:
                # http://skymapper.anu.edu.au/table-browser/dr2/
                # "flags - Bitwise OR of Source Extractor flags across all observations"
                # Sextractor flags described here: https://sextractor.readthedocs.io/en/latest/Flagging.html  # noqa
                # everything beyond flag>=4 sounds pretty bad
                (sm.flags < 4),
            )
            .group_by(sm.object_id)   # avoid duplicates - we trust the sm2 ids
            .having(
                # each and every match to the tic must satisfy the bright star rejection criteria
                fn.sum(
                    peewee.Case(
                        None,
                        (
                            (tic.gaiamag < self.parameters['gaia_g_mag_limit'], 1),
                            (tic.gaiarp < self.parameters['gaia_rp_mag_limit'], 1),
                            (tic.tmag < self.parameters['tic_t_mag_limit'], 1),
                        ),
                        0)
                ) == 0
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
#
# END BhmSpidersAgnSkyMapperDr2Carton
# ##################################################################################


class BhmSpidersAgnSuperCosmosCarton(BaseCarton):

    name = 'bhm_spiders_agn_supercosmos'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = EROSITASupersetAGN.alias()
        sc = SuperCosmos.alias()
        c2sc = CatalogToSuperCosmos.alias()
        c2cw = CatalogToCatWISE2020.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        instrument = peewee.Value(self.instrument)

        b_psf_mag_min = self.parameters['b_psf_mag_min']
        r_psf_mag_min = self.parameters['r_psf_mag_min']
        i_psf_mag_min = self.parameters['i_psf_mag_min']
        b_psf_mag_max = self.parameters['b_psf_mag_max']
        r_psf_mag_max = self.parameters['r_psf_mag_max']
        i_psf_mag_max = self.parameters['i_psf_mag_max']
        b_psf_mag_max_for_cadence1 = self.parameters['b_psf_mag_max_for_cadence1']
        r_psf_mag_max_for_cadence1 = self.parameters['r_psf_mag_max_for_cadence1']
        i_psf_mag_max_for_cadence1 = self.parameters['i_psf_mag_max_for_cadence1']
        b_psf_mag_max_for_cadence2 = self.parameters['b_psf_mag_max_for_cadence2']
        r_psf_mag_max_for_cadence2 = self.parameters['r_psf_mag_max_for_cadence2']
        i_psf_mag_max_for_cadence2 = self.parameters['i_psf_mag_max_for_cadence2']

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # #########################################################################
        # prepare the spectroscopy catalogues

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )
        # #########################################################################

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

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
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                (sph.c.pkey.is_null(False), 1),
            ),
            0)

        priority = fn.max(
            self.parameters['priority_floor'] +
            priority_1 * self.parameters['dpriority_match_flags'] +
            priority_2 * self.parameters['dpriority_det_like'] +
            priority_3 * self.parameters['dpriority_has_spec']
        )

        # choose cadence based on psf-like magnitude in supercosmos b,r2,i-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence = peewee.Case(
            None,
            (
                ((sc.scormagb.between(b_psf_mag_min, b_psf_mag_max_for_cadence1)) |
                 (sc.scormagr2.between(r_psf_mag_min, r_psf_mag_max_for_cadence1)) |
                 (sc.scormagi.between(i_psf_mag_min, i_psf_mag_max_for_cadence1)),
                 cadence1),
                ((sc.scormagb.between(b_psf_mag_max_for_cadence1, b_psf_mag_max_for_cadence2)) |
                 (sc.scormagr2.between(r_psf_mag_max_for_cadence1, r_psf_mag_max_for_cadence2)) |
                 (sc.scormagi.between(i_psf_mag_max_for_cadence1, i_psf_mag_max_for_cadence2)),
                 cadence2),
            ),
            cadence3)

        # We only use psfmags
        opt_prov = peewee.Value('sc_psfmag')
        # - transform the photographic B,R,I -> to griz
        # some very crude by-eye fits to SPIDERS AGN targets matched to SC and SDSSdr9 (via stilts)
        # completely ignore differences between psfmags and total mags
        magnitude_g = (sc.scormagb + 0.1 + (sc.scormagb - sc.scormagr2) * -0.23).cast('float')
        magnitude_r = (sc.scormagr2 + 0.25 + (sc.scormagb - sc.scormagr2) * 0.12).cast('float')
        magnitude_i = (fn.greatest(sc.scormagr2, sc.scormagi) + 0.1).cast('float')
        magnitude_z = (
            fn.greatest(sc.scormagr2, sc.scormagi) +
            0.2 + (sc.scormagb - sc.scormagr2) * -0.2 +
            -0.28 * (sc.scormagb - sc.scormagr2) * (sc.scormagb - sc.scormagr2)
        ).cast('float')

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                # fn.array_agg(sc.objid, coerce=False).alias('sc_objids'),
                x.catwise2020_id.alias('cw2020_source_id'),
                fn.min(tic.gaia_int).alias('gaia_source'),
                fn.min(x.ero_detuid).alias('ero_detuid'),
                fn.min(c.ra).alias('ra'),
                fn.min(c.dec).alias('dec'),
                priority.alias("priority"),
                fn.min(value).alias('value'),
                fn.min(cadence).alias('cadence'),
                fn.min(instrument).alias('instrument'),
                fn.min(magnitude_g).alias('g'),
                fn.min(magnitude_r).alias('r'),
                fn.min(magnitude_i).alias('i'),
                fn.min(magnitude_z).alias('z'),
                fn.min(sc.scormagb).alias('scormagb'),
                fn.min(sc.scormagr2).alias('scormagr2'),
                fn.min(sc.scormagi).alias('scormagi'),
                fn.min(sc.classmagb).alias('classmagb'),
                fn.min(sc.classmagr2).alias('classmagr2'),
                fn.min(sc.classmagi).alias('classmagi'),
                fn.min(sc.classr1).alias('classr1'),
                fn.min(tic.gaiamag).alias('gaia_g'),
                fn.min(tic.gaiabp).alias('bp'),
                fn.min(tic.gaiarp).alias('rp'),
                fn.min(opt_prov).alias('opt_prov'),
            )
            .join(c2cw)
            .join(x, on=(c2cw.target_id == x.catwise2020_id))
            .switch(c)
            .join(c2sc)
            .join(sc)
            .where(
                c.version_id == version_id,
                c2cw.version_id == version_id,
                c2sc.version_id == version_id,
                c2cw.best >> True,
                c2sc.best >> True,
            )
            .switch(c)
            .join(
                c2tic, JOIN.LEFT_OUTER,
                on=(
                    (c.catalogid == c2tic.catalogid) &
                    (c2tic.version_id == version_id)
                )
            )
            .join(tic, JOIN.LEFT_OUTER)
            .switch(c)
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # finished joining the spectroscopy
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (x.xmatch_version == self.parameters['xmatch_version']),
                (x.opt_cat == self.parameters['opt_cat']),
                (x.xmatch_metric >= self.parameters['p_any_min']),
                (x.ero_det_like > self.parameters['det_like_min']),
                (  # include anything that falls in suitable mag range
                    sc.scormagb.between(b_psf_mag_min, b_psf_mag_max) |
                    sc.scormagr2.between(r_psf_mag_min, r_psf_mag_max) |
                    sc.scormagi.between(i_psf_mag_min, i_psf_mag_max)
                ),
                # exclude anything that strays too bright:
                ~sc.scormagb.between(-99.0, b_psf_mag_min),
                ~sc.scormagr2.between(-99.0, r_psf_mag_min),
                ~sc.scormagi.between(-99.0, i_psf_mag_min),
            )
            .group_by(x.catwise2020_id)   # avoid duplicates - we trust the CatWISE2020 ids
            .having(
                # each and every match to the tic must satisfy the bright star rejection criteria
                fn.sum(
                    peewee.Case(
                        None,
                        (
                            (tic.gaiamag < self.parameters['gaia_g_mag_limit'], 1),
                            (tic.gaiarp < self.parameters['gaia_rp_mag_limit'], 1),
                            (tic.tmag < self.parameters['tic_t_mag_limit'], 1),
                        ),
                        0)
                ) == 0
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
#
# END BhmSpidersAgnSuperCosmosCarton
# ##################################################################################
