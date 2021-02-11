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
    BHM_eFEDS_Veto,
    SDSSV_BOSS_SPALL,
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
    # CatalogToPanstarrs1,    # only exists after v0.5 cross-match
)


from target_selection.mag_flux import AB2nMgy, AB2Jy
# from target_selection.mag_flux import psfmag_minus_fiber2mag


# ############################################
# ############################################
# ############################################
# ############################################
# This provides the following BHM SPIDERS AGN cartons in v0.5:
#  *  bhm_spiders_agn_lsdr8
#  *  bhm_spiders_agn_efeds_stragglers
#  *  bhm_spiders_agn_gaiadr2
#  *  bhm_spiders_agn_ps1dr2
#     bhm_spiders_agn_skymapperdr2
#     bhm_spiders_agn_supercosmos
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

    '''

    DROP TABLE sandbox.temp_td_bhm_spiders_agn_lsdr8;

    SELECT
         MAX(c.catalogid) AS catalogid,
         ls.ls_id AS ls_id,
         MAX(x.ero_detuid) AS ero_detuid,
         MAX(c.ra) AS ra,
         MAX(c.dec) AS dec,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.fiberflux_g))) as fibermag_g,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.fiberflux_r))) as fibermag_r,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.fiberflux_z))) as fibermag_z,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.fibertotflux_g))) as fibertotmag_g,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.fibertotflux_r))) as fibertotmag_r,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.fibertotflux_z))) as fibertotmag_z,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.flux_g))) as mag_g,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.flux_r))) as mag_r,
         MIN(22.5-2.5*log10(greatest(1e-3,ls.flux_z))) as mag_z,
         (CASE WHEN MAX(s16.specobjid) IS NOT NULL THEN 1
               WHEN MAX(s2020.pk) IS NOT NULL THEN 1
               WHEN MAX(sV.specobjid) IS NOT NULL THEN 1
               ELSE 0 END) AS has_spec,
         (1520 +
          (CASE WHEN MAX(s16.specobjid) IS NOT NULL THEN 4
                WHEN MAX(s2020.pk)      IS NOT NULL THEN 4
                WHEN MAX(sV.specobjid)  IS NOT NULL THEN 4
                ELSE 0 END) +
          (CASE WHEN MAX(x.xmatch_flags) > 1 THEN 1
                ELSE 0 END) +
          (CASE WHEN MAX(x.ero_det_like) < 8.0 THEN 2
                ELSE 0 END)
         ) AS priority,
         (CASE WHEN (     MAX(ls.fibertotflux_r) > 251.189
                       OR MAX(ls.gaia_phot_g_mean_mag) BETWEEN 0.1 AND 16.0
                       OR MAX(ls.gaia_phot_rp_mean_mag) BETWEEN 0.1 AND 16.0
                     ) THEN 'bright_2x1'
               WHEN MAX(ls.fibertotflux_r) > 39.8107 THEN 'dark_1x2'
               WHEN MAX(ls.fibertotflux_r) <= 39.8107 THEN 'dark_1x4'
               ELSE 'unknown_cadence' END) AS cadence,
         MAX(ls.gaia_phot_g_mean_mag) as gaia_g,
         MAX(ls.gaia_phot_rp_mean_mag) as gaia_rp,
         MAX(ls.ref_cat) AS ref_cat,
         MAX(ls.maskbits) AS ls_maskbits,
         MAX(CASE WHEN ls.type = 'PSF' THEN 'ls_psfmag' ELSE 'ls_fibertotmag' END) AS opt_prov,
         MAX(x.xmatch_flags) AS xmatch_flags,
         MAX(x.ero_det_like) AS ero_det_like

    INTO sandbox.temp_td_bhm_spiders_agn_lsdr8
    FROM catalogdb.catalog AS c
    JOIN catalog_to_legacy_survey_dr8 as c2ls
         ON c.catalogid = c2ls.catalogid
    JOIN legacy_survey_dr8 AS ls
         ON c2ls.target_id = ls.ls_id
    JOIN erosita_superset_agn AS x
         ON x.ls_id = ls.ls_id
    LEFT OUTER JOIN catalog_to_sdss_dr16_specobj AS c2s16
         ON c.catalogid = c2s16.catalogid
    LEFT OUTER JOIN sdss_dr16_specobj AS s16
         ON ( c2s16.target_id = s16.specobjid
              AND s16.zwarning = 0
              AND s16.snmedian > 1.6
              AND s16.zerr < 0.002
              AND s16.zerr > 0.0
              AND s16.scienceprimary > 0 )
    LEFT OUTER JOIN bhm_efeds_veto AS s2020
          ON ( q3c_join(s2020.plug_ra,s2020.plug_dec,c.ra,c.dec,1.0/3600.)
               AND s2020.zwarning = 0
               AND s2020.sn_median_all > 1.6
               AND s2020.z_err < 0.002
               AND s2020.z_err > 0.0
             )
    LEFT OUTER JOIN sdssv_boss_spall AS sV
          ON ( q3c_join(sV.plug_ra,sV.plug_dec,c.ra,c.dec,1.0/3600.)
               AND sV.zwarning = 0
               AND sV.sn_median_all > 1.6
               AND sV.z_err < 0.002
               AND sV.z > 0.0
             )
    WHERE
            x.ero_version = 'em01_c946_201008_poscorr'
        AND x.xmatch_method = 'XPS/NWAY'
        AND x.xmatch_version = 'JWMS_v_40'
        AND x.opt_cat = 'lsdr8'
        AND x.xmatch_metric > 0.1
        AND x.ero_det_like > 6.0
        AND ((ls.fibertotflux_r BETWEEN 1.0 AND 3981.07 ) OR
             (ls.fibertotflux_z BETWEEN 3.98107 AND 3981.07))
        AND (ls.gaia_phot_g_mean_mag NOT BETWEEN 0.1 AND 13.5 )
        AND (ls.gaia_phot_rp_mean_mag NOT BETWEEN 0.1 AND 13.5 )
        AND ls.nobs_r > 0
        AND ((ls.nobs_g > 0) OR (ls.nobs_z > 0))
        AND (ls.maskbits & 28) = 0
        AND c.version_id = 21
        AND c2ls.version_id = 21
        AND c2ls.best IS TRUE
    GROUP BY ls.ls_id
    ;


    select cadence,count(*) from sandbox.temp_td_bhm_spiders_agn_lsdr8 group by cadence;
    select priority,count(*) from sandbox.temp_td_bhm_spiders_agn_lsdr8 group by priority;
    select count(*) from sandbox.temp_td_bhm_spiders_agn_lsdr8 where has_spec > 0;

    select cadence,count(*) from sandbox.temp_bhm_spiders_agn_lsdr8 group by cadence;
    select priority,count(*) from sandbox.temp_bhm_spiders_agn_lsdr8 group by priority;

    \copy (SELECT * FROM sandbox.temp_td_bhm_spiders_agn_lsdr8) TO '/home/dwelly/scratch/bhm_spiders_agn_lsdr8.csv' with csv header

    stilts tpipe \
    in=/home/dwelly/scratch/bhm_spiders_agn_lsdr8.csv \
    out=/home/dwelly/scratch/bhm_spiders_agn_lsdr8.fits ofmt=fits-basic

    '''

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
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        s16 = SDSS_DR16_SpecObj.alias()
        s2020 = BHM_eFEDS_Veto.alias()
        sV = SDSSV_BOSS_SPALL.alias()

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

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

        p_f = self.parameters['priority_floor']
        priority_1 = peewee.Case(
            None, ((fn.max(x.xmatch_flags) > 1, 1), ), 0)
        priority_2 = peewee.Case(
            None, ((fn.max(x.ero_det_like) < self.parameters['det_like_for_priority'], 2), ), 0)
        priority_3 = peewee.Case(
            None,
            (
                (fn.max(s16.specobjid).is_null(False), 4),  # any of these can be satisfied
                (fn.max(s2020.pk).is_null(False), 4),
                (fn.max(sV.specobjid).is_null(False), 4),
            ),
            0)
        priority = (p_f + priority_1 + priority_2 + priority_3)

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

        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

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
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(  # rely on catalogdb.catalog cross-matches to keep processing time down
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.specobjid) &
                    (s16.snmedian >= spec_sn_thresh) &
                    (s16.zwarning == 0) &
                    (s16.zerr <= spec_z_err_thresh) &
                    (s16.zerr > 0.0) &
                    (s16.scienceprimary > 0)
                )
            )
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(s2020.plug_ra, s2020.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (s2020.sn_median_all >= spec_sn_thresh) &
                    (s2020.zwarning == 0) &
                    (s2020.z_err <= spec_z_err_thresh) &
                    (s2020.z_err > 0.0)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.plug_ra, sV.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (sV.sn_median_all >= spec_sn_thresh) &
                    (sV.zwarning == 0) &
                    (sV.z_err <= spec_z_err_thresh) &
                    (sV.z_err > 0.0)
                )
            )
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True
            )
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

    '''

    DROP TABLE sandbox.temp_td_bhm_spiders_agn_efeds_stragglers;

    SELECT
         MAX(c.catalogid) AS catalogid,
         ls.ls_id AS ls_id,
         MAX(x.ero_detuid) AS ero_detuid,
         MAX(c.ra) AS ra,
         MAX(c.dec) AS dec,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.fiberflux_g)))) as fibermag_g,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.fiberflux_r)))) as fibermag_r,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.fiberflux_z)))) as fibermag_z,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.fibertotflux_g)))) as fibertotmag_g,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.fibertotflux_r)))) as fibertotmag_r,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.fibertotflux_z)))) as fibertotmag_z,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.flux_g)))) as mag_g,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.flux_r)))) as mag_r,
         (22.5-2.5*log10(greatest(1e-3,AVG(ls.flux_z)))) as mag_z,
         (CASE WHEN MAX(s16.specobjid) IS NOT NULL THEN 1
               WHEN MAX(s2020.pk)      IS NOT NULL THEN 1
               WHEN MAX(sV.specobjid)  IS NOT NULL THEN 1
               WHEN     (MAX(ph.pkey)       IS NOT NULL
                    AND  MAX(phm.yanny_uid) IS NOT NULL) THEN 1
               ELSE 0 END) AS has_spec,
         (1510 +
           (CASE WHEN MAX(s16.specobjid)  IS NOT NULL  THEN 10000
                 WHEN MAX(s2020.pk)       IS NOT NULL  THEN 10000
                 WHEN MAX(sV.specobjid)   IS NOT NULL  THEN 10000
                 WHEN (MAX(ph.pkey)       IS NOT NULL
                   AND MAX(phm.yanny_uid) IS NOT NULL) THEN 10000
                 ELSE 0 END) +
           (CASE WHEN MAX(x.xmatch_flags) > 1 THEN 1
                 ELSE 0 END) +
           (CASE WHEN MAX(x.ero_det_like) < 8.0 THEN 2
                 ELSE 0 END) +
           MIN(CASE WHEN x.ero_version = 'eFEDS_c001_V18C_V3_ext' THEN 4
                 ELSE 0 END)
         ) AS priority,
         (CASE WHEN (     MAX(ls.fibertotflux_r) > 251.189
                       OR MAX(ls.gaia_phot_g_mean_mag) BETWEEN 0.1 AND 16.0
                       OR MAX(ls.gaia_phot_rp_mean_mag) BETWEEN 0.1 AND 16.0
                     ) THEN 'bright_2x1'
               WHEN MAX(ls.fibertotflux_r) > 39.8107 THEN 'dark_1x2'
               WHEN MAX(ls.fibertotflux_r) <= 39.8107 THEN 'dark_1x4'
               ELSE 'unknown_cadence' END) AS cadence,
         MAX(ls.gaia_phot_g_mean_mag) as gaia_g,
         MAX(ls.gaia_phot_rp_mean_mag) as gaia_rp,
         MAX(ls.ref_cat) AS ref_cat,
         MAX(ls.maskbits) AS ls_maskbits,
         MAX(CASE WHEN ls.type = 'PSF' THEN 'ls_psfmag' ELSE 'ls_fibertotmag' END) AS opt_prov,
         MAX(x.xmatch_flags) AS xmatch_flags,
         MAX(x.ero_det_like) AS ero_det_like,
         array_agg(x.ero_version) AS ero_version

    INTO sandbox.temp_td_bhm_spiders_agn_efeds_stragglers
    FROM catalogdb.catalog AS c
    JOIN catalog_to_legacy_survey_dr8 as c2ls
         ON c.catalogid = c2ls.catalogid
    JOIN legacy_survey_dr8 AS ls
         ON c2ls.target_id = ls.ls_id
    JOIN erosita_superset_agn AS x
         ON x.ls_id = ls.ls_id
    LEFT OUTER JOIN catalog_to_sdss_dr16_specobj AS c2s16
         ON c.catalogid = c2s16.catalogid
    LEFT OUTER JOIN sdss_dr16_specobj AS s16
         ON ( c2s16.target_id = s16.specobjid
              AND s16.zwarning = 0
              AND s16.snmedian > 1.6
              AND s16.zerr < 0.002
              AND s16.zerr > 0.0
              AND s16.scienceprimary > 0 )
    LEFT OUTER JOIN bhm_efeds_veto AS s2020
          ON ( q3c_join(s2020.plug_ra,s2020.plug_dec,c.ra,c.dec,1.0/3600.)
               AND s2020.zwarning = 0
               AND s2020.sn_median_all > 1.6
               AND s2020.z_err < 0.002
               AND s2020.z_err > 0.0
             )
    LEFT OUTER JOIN sdssv_boss_spall AS sV
          ON ( q3c_join(sV.plug_ra,sV.plug_dec,c.ra,c.dec,1.0/3600.)
               AND sV.zwarning = 0
               AND sV.sn_median_all > 1.6
               AND sV.z_err < 0.002
               AND sV.z > 0.0
             )
    LEFT OUTER JOIN sdssv_plateholes AS ph
          ON ( q3c_join(ph.target_ra,ph.target_dec,c.ra,c.dec,1.0/3600.)
               AND ph.holetype = 'BOSS_SHARED'
               AND (ph.sourcetype = 'SCI' OR ph.sourcetype = 'STA' )
             )
    LEFT OUTER JOIN sdssv_plateholes_meta AS phm
          ON (     ph.yanny_uid = phm.yanny_uid )
    WHERE
            ((     x.ero_version = 'eFEDS_c001_V18C_V3_main'
               AND x.xmatch_method = 'XPS/NWAY'
               AND x.xmatch_version = 'Merged_03DEC2020' ) OR
             (     x.ero_version = 'eFEDS_c001_V18C_V3_ext'
               AND x.xmatch_method = 'XPS/NWAY'
               AND x.xmatch_version = 'Merged_03DEC2020' ))
        AND x.opt_cat = 'lsdr8'
        AND x.xmatch_metric > 0.1
        AND x.ero_det_like > 6.0
        AND ((ls.fibertotflux_r BETWEEN 1.0 AND 3981.07 ) OR
             (ls.fibertotflux_z BETWEEN 3.98107 AND 3981.07))
        AND (ls.gaia_phot_g_mean_mag NOT BETWEEN 0.1 AND 13.5 )
        AND (ls.gaia_phot_rp_mean_mag NOT BETWEEN 0.1 AND 13.5 )
        AND ls.nobs_r > 0
        AND ( (ls.nobs_g > 0) OR (ls.nobs_z > 0) )
        AND ( ls.maskbits & 28 )  = 0
        AND c.version_id = 21
        AND c2ls.version_id = 21
        AND c2ls.best IS TRUE
    GROUP BY ls.ls_id
    ;


    #########################
    ##### TODO add this back in when isvalid column is added to sdssv_plateholes_meta
               AND phm.isvalid IS TRUE )
    ###########################

    select cadence,count(*) from sandbox.temp_td_bhm_spiders_agn_efeds_stragglers group by cadence;
    select priority,count(*) from sandbox.temp_td_bhm_spiders_agn_efeds_stragglers group by priority ORDER BY priority;
    select count(*) from sandbox.temp_td_bhm_spiders_agn_efeds_stragglers where has_spec > 0;

    select cadence,count(*) from sandbox.temp_bhm_spiders_agn_efeds_stragglers group by cadence;
    select priority,count(*) from sandbox.temp_bhm_spiders_agn_efeds_stragglers group by priority ORDER BY priority;

    \copy (SELECT * FROM sandbox.temp_td_bhm_spiders_agn_efeds_stragglers) TO '/home/dwelly/scratch/bhm_spiders_agn_efeds_stragglers.csv' with csv header

    stilts tpipe \
    in=/home/dwelly/scratch/bhm_spiders_agn_efeds_stragglers.csv \
    out=/home/dwelly/scratch/bhm_spiders_agn_efeds_stragglers.fits ofmt=fits-basic

    '''

    name = 'bhm_spiders_agn_efeds_stragglers'
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
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        s16 = SDSS_DR16_SpecObj.alias()
        s2020 = BHM_eFEDS_Veto.alias()
        sV = SDSSV_BOSS_SPALL.alias()
        ph = SDSSV_Plateholes.alias()
        phm = SDSSV_Plateholes_Meta.alias()

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

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

        p_f = self.parameters['priority_floor']
        priority_1 = peewee.Case(
            None, ((fn.max(x.xmatch_flags) > 1, 1), ), 0)
        priority_2 = peewee.Case(
            None, ((fn.max(x.ero_det_like) < self.parameters['det_like_for_priority'], 2), ), 0)
        priority_3 = peewee.Case(
            None,
            (
                (fn.max(s16.specobjid).is_null(False), 10000),  # any of these can be satisfied
                (fn.max(s2020.pk).is_null(False), 10000),
                (fn.max(sV.specobjid).is_null(False), 10000),
                (fn.max(ph.pkey).is_null(False) & fn.max(phm.yanny_uid ).is_null(False), 10000),
            ),
            0)
        # add a step down in priority for anything selected by the secondary xmatch_version
        priority_4 = fn.min(peewee.Case(
            None, ((x.ero_version == self.parameters['ero_version2'], 8), ), 0))

        priority = (p_f + priority_1 + priority_2 + priority_3 + priority_4)

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

        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

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
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(  # rely on catalogdb.catalog cross-matches to keep processing time down
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.specobjid) &
                    (s16.snmedian >= spec_sn_thresh) &
                    (s16.zwarning == 0) &
                    (s16.zerr <= spec_z_err_thresh) &
                    (s16.zerr > 0.0) &
                    (s16.scienceprimary > 0)
                )
            )
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(s2020.plug_ra, s2020.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (s2020.sn_median_all >= spec_sn_thresh) &
                    (s2020.zwarning == 0) &
                    (s2020.z_err <= spec_z_err_thresh) &
                    (s2020.z_err > 0.0)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.plug_ra, sV.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (sV.sn_median_all >= spec_sn_thresh) &
                    (sV.zwarning == 0) &
                    (sV.z_err <= spec_z_err_thresh) &
                    (sV.z_err > 0.0)
                )
            )
            .join(
                ph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(ph.target_ra, ph.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (ph.holetype == 'BOSS_SHARED') &
                    (
                        (ph.sourcetype == 'SCI') |
                        (ph.sourcetype == 'STA')
                    )
                )
            )
            .join(
                phm, JOIN.LEFT_OUTER,
                on=(
                    (ph.yanny_uid == phm.yanny_uid) # &
                    ## TODO add this back in when isvalid column is added to sdssv_plateholes_meta
                    ## (phm.isvalid >> True)
                )
            )
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True
            )
            .where(
                (
                    (
                        (x.ero_version == self.parameters['ero_version1']) &
                        (x.xmatch_method == self.parameters['xmatch_method1']) &
                        (x.xmatch_version == self.parameters['xmatch_version1'])
                    ) |
                    (
                        (x.ero_version == self.parameters['ero_version2']) &
                        (x.xmatch_method == self.parameters['xmatch_method2']) &
                        (x.xmatch_version == self.parameters['xmatch_version2'])
                    )
                ),
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
# END BhmSpidersAgnEfedsStragglersCarton
# ##################################################################################


class BhmSpidersAgnGaiadr2Carton(BaseCarton):

    '''

    DROP TABLE sandbox.temp_td_bhm_spiders_agn_gaiadr2;

    SELECT
         MAX(c.catalogid) AS catalogid,
         tic.gaia_int AS gaia_source,
         MAX(x.ero_detuid) AS ero_detuid,
         MAX(c.ra) AS ra,
         MAX(c.dec) AS dec,
         MIN(tic.gaiamag) as gaia_g,
         MIN(tic.gaiabp) as bp,
         MIN(tic.gaiarp) as rp,
         (1530 +
          (CASE WHEN MAX(s16.specobjid) IS NOT NULL THEN 4
                WHEN MAX(s2020.pk) IS NOT NULL THEN 4
                WHEN MAX(sV.specobjid) IS NOT NULL THEN 4
                ELSE 0 END) +
          (CASE WHEN MAX(x.xmatch_flags) > 1 THEN 1
                ELSE 0 END) +
          (CASE WHEN MAX(x.ero_det_like) < 8.0 THEN 2
                ELSE 0 END)
         ) AS priority,
         (CASE WHEN (MAX(tic.gaiamag) < 16.0 OR MAX(tic.gaiarp) < 16.0) THEN 'bright_2x1'
               WHEN (MAX(tic.gaiamag) < 19.0 OR MAX(tic.gaiarp) < 19.0) THEN 'dark_1x2'
               WHEN (MAX(tic.gaiamag) >= 19.0 AND MAX(tic.gaiarp) >= 19.0) THEN 'dark_1x4'
               ELSE 'unknown_cadence' END) AS cadence,
         MAX(x.xmatch_flags) AS xmatch_flags,
         MAX(x.ero_det_like) AS ero_det_like

    INTO sandbox.temp_td_bhm_spiders_agn_gaiadr2
    FROM catalogdb.catalog AS c
    JOIN catalog_to_TIC_v8 as c2tic
         ON c.catalogid = c2tic.catalogid
    JOIN tic_v8 AS tic
         ON c2tic.target_id = tic.id
    JOIN erosita_superset_agn AS x
         ON x.gaia_dr2_id = tic.gaia_int
    LEFT OUTER JOIN catalog_to_sdss_dr16_specobj AS c2s16
         ON c.catalogid = c2s16.catalogid
    LEFT OUTER JOIN sdss_dr16_specobj AS s16
         ON ( c2s16.target_id = s16.specobjid
              AND s16.zwarning = 0
              AND s16.snmedian > 1.6
              AND s16.zerr < 0.002
              AND s16.zerr > 0.0
              AND s16.scienceprimary > 0 )
    LEFT OUTER JOIN bhm_efeds_veto AS s2020
          ON ( q3c_join(s2020.plug_ra,s2020.plug_dec,c.ra,c.dec,1.0/3600.)
               AND s2020.zwarning = 0
               AND s2020.sn_median_all > 1.6
               AND s2020.z_err < 0.002
               AND s2020.z_err > 0.0
             )
    LEFT OUTER JOIN sdssv_boss_spall AS sV
          ON ( q3c_join(sV.plug_ra,sV.plug_dec,c.ra,c.dec,1.0/3600.)
               AND sV.zwarning = 0
               AND sV.sn_median_all > 1.6
               AND sV.z_err < 0.002
               AND sV.z > 0.0
             )
    WHERE
            x.ero_version = 'em01_c946_201008_poscorr'
        AND x.xmatch_method = 'XPS/NWAY'
        AND x.xmatch_version = 'JWMS_CW2_v_03_TDopt'
        AND x.opt_cat = 'gaiadr2'
        AND x.xmatch_metric > 0.1
        AND x.ero_det_like > 6.0
        AND tic.gaiamag > 13.5
        AND tic.gaiarp > 13.5
        AND c.version_id = 21
        AND c2tic.version_id = 21
        AND c2tic.best IS TRUE
    GROUP BY tic.gaia_int
    ;


    select cadence,count(*) from sandbox.temp_td_bhm_spiders_agn_gaiadr2 group by cadence;
    select priority,count(*) from sandbox.temp_td_bhm_spiders_agn_gaiadr2 group by priority order by priority;

    select cadence,count(*) from sandbox.temp_bhm_spiders_agn_gaiadr2 group by cadence;
    select priority,count(*) from sandbox.temp_bhm_spiders_agn_gaiadr2 group by priority order by priorit;

    \copy (SELECT * FROM sandbox.temp_bhm_spiders_agn_gaiadr2) TO '/home/dwelly/scratch/bhm_spiders_agn_gaiadr2.csv' with csv header

    stilts tpipe \
    in=/home/dwelly/scratch/bhm_spiders_agn_gaiadr2.csv \
    out=/home/dwelly/scratch/bhm_spiders_agn_gaiadr2.fits ofmt=fits-basic

    '''

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
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        s16 = SDSS_DR16_SpecObj.alias()
        s2020 = BHM_eFEDS_Veto.alias()
        sV = SDSSV_BOSS_SPALL.alias()

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
        # add +4 if target has existing good SDSS spectroscopy

        p_f = self.parameters['priority_floor']
        priority_1 = peewee.Case(
            None, ((fn.max(x.xmatch_flags) > 1, 1), ), 0)
        priority_2 = peewee.Case(
            None, ((fn.max(x.ero_det_like) < self.parameters['det_like_for_priority'], 2), ), 0)
        priority_3 = peewee.Case(
            None,
            (
                (fn.max(s16.specobjid).is_null(False), 4),  # any of these can be satisfied
                (fn.max(s2020.pk).is_null(False), 4),
                (fn.max(sV.specobjid).is_null(False), 4),
            ),
            0)
        priority = (p_f + priority_1 + priority_2 + priority_3)

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

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

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
            .join(tic)
            .join(x, on=(tic.gaia_int == x.gaia_dr2_id))
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(  # rely on catalogdb.catalog cross-matches to keep processing time down
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.specobjid) &
                    (s16.snmedian >= spec_sn_thresh) &
                    (s16.zwarning == 0) &
                    (s16.zerr <= spec_z_err_thresh) &
                    (s16.zerr > 0.0) &
                    (s16.scienceprimary > 0)
                )
            )
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(s2020.plug_ra, s2020.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (s2020.sn_median_all >= spec_sn_thresh) &
                    (s2020.zwarning == 0) &
                    (s2020.z_err <= spec_z_err_thresh) &
                    (s2020.z_err > 0.0)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.plug_ra, sV.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (sV.sn_median_all >= spec_sn_thresh) &
                    (sV.zwarning == 0) &
                    (sV.z_err <= spec_z_err_thresh) &
                    (sV.z_err > 0.0)
                )
            )
            .where(
                c.version_id == version_id,
                c2tic.version_id == version_id,
                c2tic.best >> True
            )
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
    '''

    DROP TABLE sandbox.temp_td_bhm_spiders_agn_sep;

    SELECT
         MAX(c.catalogid) AS catalogid,
         tic.gaia_int AS gaia_source,
         MAX(x.ero_detuid) AS ero_detuid,
         MAX(c.ra) AS ra,
         MAX(c.dec) AS dec,
         MAX(tic.gaiamag) as gaia_g,
         MAX(tic.gaiabp) as bp,
         MAX(tic.gaiarp) as rp,
         MIN(
           1590 +
           (CASE WHEN x.xmatch_flags > 1   THEN 1 ELSE 0 END) +
           (CASE WHEN x.ero_det_like < 8.0 THEN 2 ELSE 0 END) +
           (CASE WHEN x.xmatch_version = 'SEP_CW2_NOV2020_MSopt' THEN 4 ELSE 0 END)
         ) AS priority,
         (CASE WHEN (MAX(tic.gaiamag) < 16.0 OR MAX(tic.gaiarp) < 16.0) THEN 'bright_2x1'
               WHEN (MAX(tic.gaiamag) < 19.0 OR MAX(tic.gaiarp) < 19.0) THEN 'dark_1x2'
               WHEN (MAX(tic.gaiamag) >= 19.0 AND MAX(tic.gaiarp) >= 19.0) THEN 'dark_1x4'
               ELSE 'unknown_cadence' END) AS cadence,
         MAX(x.ero_det_like) AS ero_det_like,
         array_agg(x.xmatch_version) AS xmatch_version
    INTO sandbox.temp_td_bhm_spiders_agn_sep
    FROM catalogdb.catalog AS c
    JOIN catalog_to_TIC_v8 as c2tic
         ON c.catalogid = c2tic.catalogid
    JOIN tic_v8 AS tic
         ON c2tic.target_id = tic.id
    JOIN erosita_superset_agn AS x
         ON x.gaia_dr2_id = tic.gaia_int
    WHERE
            x.ero_version = 'em01_SEP_c946'
        AND x.xmatch_method = 'XPS/NWAY'
        AND ( x.xmatch_version = 'SEP_CW2_07DEC2020_TDopt'
           OR x.xmatch_version = 'SEP_CW2_NOV2020_MSopt' )
        AND x.opt_cat = 'gaiadr2'
        AND x.xmatch_metric > 0.1
        AND x.ero_det_like > 6.0
        AND tic.gaiamag > 13.5
        AND tic.gaiarp > 13.5
        AND c.version_id = 21
        AND c2tic.version_id = 21
        AND c2tic.best IS TRUE
    GROUP BY tic.gaia_int
    ;


    select cadence,count(*) from sandbox.temp_td_bhm_spiders_agn_sep group by cadence;
    select priority,count(*) from sandbox.temp_td_bhm_spiders_agn_sep group by priority order by priority;

    select cadence,count(*) from sandbox.temp_bhm_spiders_agn_sep group by cadence;
    select priority,count(*) from sandbox.temp_bhm_spiders_agn_sep group by priority order by priority;

    \copy (SELECT * FROM sandbox.temp_bhm_spiders_agn_sep) TO '/home/dwelly/scratch/bhm_spiders_agn_sep.csv' with csv header

    stilts tpipe \
    in=/home/dwelly/scratch/bhm_spiders_agn_sep.csv \
    out=/home/dwelly/scratch/bhm_spiders_agn_sep.fits ofmt=fits-basic

    '''

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
        # add +4 if target has existing good SDSS spectroscopy

        p_f = self.parameters['priority_floor']
        priority_1 = peewee.Case(
            None, ((x.xmatch_flags > 1, 1), ), 0)
        priority_2 = peewee.Case(
            None, ((x.ero_det_like < self.parameters['det_like_for_priority'], 2), ), 0)
        # add a step down in priority for anything selected by the secondary xmatch_version
        priority_3 = peewee.Case(
            None, ((x.xmatch_method == self.parameters['xmatch_version2'], 4), ), 0)
        # choose the minimum priority option for all versions of this target
        priority = fn.min(p_f + priority_1 + priority_2 + priority_3)

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
            .join(tic)
            .join(x, on=(tic.gaia_int == x.gaia_dr2_id))
            .switch(c)
            .where(
                c.version_id == version_id,
                c2tic.version_id == version_id,
                c2tic.best >> True
            )
            .where(
                (x.ero_version == self.parameters['ero_version']),
                (x.xmatch_method == self.parameters['xmatch_method']),
                (
                    (x.xmatch_version == self.parameters['xmatch_version1']) |
                    (x.xmatch_version == self.parameters['xmatch_version2'])
                )
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
# END BhmSpidersAgnSepCarton
# ##################################################################################


class BhmSpidersAgnPs1dr2Carton(BaseCarton):

    '''

    DROP TABLE sandbox.temp_td_bhm_spiders_agn_ps1dr2;

    SELECT
         MAX(c.catalogid) AS catalogid,
         ps.catid_objid AS ps1dr2_id,
         MAX(x.ero_detuid) AS ero_detuid,
         MAX(c.ra) AS ra,
         MAX(c.dec) AS dec,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.g_stk_psf_flux)))  AS psf_g,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.r_stk_psf_flux)))  AS psf_r,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.i_stk_psf_flux)))  AS psf_i,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.z_stk_psf_flux)))  AS psf_z,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.g_stk_aper_flux))) AS aper_g,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.r_stk_aper_flux))) AS aper_r,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.i_stk_aper_flux))) AS aper_i,
         MAX(8.9-2.5*log10(greatest(3.631e-9,ps.z_stk_aper_flux))) AS aper_z,
         MIN(tic.gaiamag) as gaia_g,
         MIN(tic.gaiabp) as bp,
         MIN(tic.gaiarp) as rp,
         MAX(CASE WHEN (ps.flags & (8388608 + 16777216)) = 0 THEN 'ps_psfmag'
                  ELSE 'ps_apermag' END) AS opt_prov,
         (1530 +
          (CASE WHEN MAX(s16.specobjid) IS NOT NULL THEN 4
                WHEN MAX(s2020.pk)      IS NOT NULL THEN 4
                WHEN MAX(sV.specobjid)  IS NOT NULL THEN 4
                ELSE 0 END) +
          (CASE WHEN MAX(x.xmatch_flags) > 1 THEN 1
                ELSE 0 END) +
          (CASE WHEN MAX(x.ero_det_like) < 8.0 THEN 2
                ELSE 0 END)
         ) AS priority,
         (CASE WHEN (   MAX(ps.g_stk_psf_flux) > 9.120e-4
                     OR MAX(ps.r_stk_psf_flux) > 9.120e-4
                     OR MAX(ps.i_stk_psf_flux) > 9.120e-4) THEN 'bright_2x1'
               WHEN (   MAX(ps.g_stk_psf_flux) > 1.445e-4
                     OR MAX(ps.r_stk_psf_flux) > 1.445e-4
                     OR MAX(ps.i_stk_psf_flux) > 1.445e-4) THEN 'dark_1x2'
               WHEN (    MAX(ps.g_stk_psf_flux) < 1.445e-4
                     AND MAX(ps.r_stk_psf_flux) < 1.445e-4
                     AND MAX(ps.i_stk_psf_flux) < 1.445e-4) THEN 'dark_1x4'
               ELSE 'unknown_cadence' END) AS cadence,
         MAX(x.xmatch_flags) AS xmatch_flags,
         MAX(x.ero_det_like) AS ero_det_like

    INTO sandbox.temp_td_bhm_spiders_agn_ps1dr2
    FROM catalogdb.catalog AS c
    JOIN catalog_to_panstarrs1 as c2ps
         ON c.catalogid = c2ps.catalogid
    JOIN panstarrs1 AS ps
         ON c2ps.target_id = ps.catid_objid
    JOIN erosita_superset_agn AS x
         ON x.ps1_dr2_id = ps.catid_objid
    LEFT OUTER JOIN catalog_to_TIC_v8 as c2tic
         ON (c.catalogid = c2tic.catalogid AND c2tic.version_id = 21 )
    LEFT OUTER JOIN tic_v8 AS tic
         ON c2tic.target_id = tic.id
    LEFT OUTER JOIN catalog_to_sdss_dr16_specobj AS c2s16
         ON c.catalogid = c2s16.catalogid
    LEFT OUTER JOIN sdss_dr16_specobj AS s16
         ON ( c2s16.target_id = s16.specobjid
              AND s16.zwarning = 0
              AND s16.snmedian > 1.6
              AND s16.zerr < 0.002
              AND s16.zerr > 0.0
              AND s16.scienceprimary > 0 )
    LEFT OUTER JOIN bhm_efeds_veto AS s2020
          ON ( q3c_join(s2020.plug_ra,s2020.plug_dec,c.ra,c.dec,1.0/3600.)
               AND s2020.zwarning = 0
               AND s2020.sn_median_all > 1.6
               AND s2020.z_err < 0.002
               AND s2020.z_err > 0.0
             )
    LEFT OUTER JOIN sdssv_boss_spall AS sV
          ON ( q3c_join(sV.plug_ra,sV.plug_dec,c.ra,c.dec,1.0/3600.)
               AND sV.zwarning = 0
               AND sV.sn_median_all > 1.6
               AND sV.z_err < 0.002
               AND sV.z > 0.0
             )
    WHERE
            x.ero_version = 'em01_c946_201008_poscorr'
        AND x.xmatch_method = 'XPS/NWAY'
        AND x.xmatch_version = 'JWMS_CW2_v_03_TDopt'
        AND x.opt_cat = 'ps1dr2'
        AND x.xmatch_metric > 0.1
        AND x.ero_det_like > 6.0
        AND ps.g_stk_psf_flux < 1.445e-2
        AND ps.r_stk_psf_flux < 1.445e-2
        AND ps.i_stk_psf_flux < 1.445e-2
        AND ps.g_stk_psf_flux != 'NaN'
        AND ps.r_stk_psf_flux != 'NaN'
        AND ps.i_stk_psf_flux != 'NaN'
        AND (ps.flags & 134217728 ) > 0
        AND c.version_id = 21
        AND c2ps.version_id = 21
        AND c2ps.best IS TRUE
    GROUP BY ps.catid_objid
    HAVING
      SUM(CASE WHERE
           ( tic.gaiamag < 13.5
          OR tic.gaiarp < 13.5
          OR tic.tmag < 13.0
           ) THEN 1 ELSE 0 END
         ) = 0
    ;

    #    AND q3c_radial_query(c.ra,c.dec,135.0,+1.0,1.0)

    select cadence,count(*) from sandbox.temp_td_bhm_spiders_agn_ps1dr2 group by cadence;
    select priority,count(*) from sandbox.temp_td_bhm_spiders_agn_ps1dr2 group by priority order by priority;

    select cadence,count(*) from sandbox.temp_bhm_spiders_agn_ps1dr2 group by cadence;
    select priority,count(*) from sandbox.temp_bhm_spiders_agn_ps1dr2 group by priority order by priority;

    \copy (SELECT * FROM sandbox.temp_bhm_spiders_agn_ps1dr2) TO '/home/dwelly/scratch/bhm_spiders_agn_ps1dr2.csv' with csv header

    stilts tpipe \
    in=/home/dwelly/scratch/bhm_spiders_agn_ps1dr2.csv \
    out=/home/dwelly/scratch/bhm_spiders_agn_ps1dr2.fits ofmt=fits-basic

    '''

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
        c2ps = CatalogToPanstarrs1.alias()   # only exists after v0.5 cross-match
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        s16 = SDSS_DR16_SpecObj.alias()
        s2020 = BHM_eFEDS_Veto.alias()
        sV = SDSSV_BOSS_SPALL.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        instrument = peewee.Value(self.instrument)

        g_psf_flux_max = AB2Jy(self.parameters['g_psf_mag_min'])
        r_psf_flux_max = AB2Jy(self.parameters['r_psf_mag_min'])
        i_psf_flux_max = AB2Jy(self.parameters['i_psf_mag_min'])
        g_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['g_psf_mag_max_for_cadence1'])
        r_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence1'])
        i_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence1'])
        g_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['g_psf_mag_max_for_cadence2'])
        r_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence2'])
        i_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence2'])

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')

        # priority is determined by target properties
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:
        # add +1 if target is a secondary cross-match (match_flag > 1)
        # add +2 if target has a low value of ero_det_like
        # add +4 if target has existing good SDSS spectroscopy

        p_f = self.parameters['priority_floor']
        priority_1 = peewee.Case(
            None, ((fn.max(x.xmatch_flags) > 1, 1), ), 0)
        priority_2 = peewee.Case(
            None, ((fn.max(x.ero_det_like) < self.parameters['det_like_for_priority'], 2), ), 0)
        priority_3 = peewee.Case(
            None,
            (
                (fn.max(s16.specobjid).is_null(False), 4),  # any of these can be satisfied
                (fn.max(s2020.pk).is_null(False), 4),
                (fn.max(sV.specobjid).is_null(False), 4),
            ),
            0)
        priority = (p_f + priority_1 + priority_2 + priority_3)

        # choose cadence based on psf_flux magnitude in panstarrs1 g,r,i-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((ps.g_stk_psf_flux < g_psf_flux_min_for_cadence1) |
                 (ps.r_stk_psf_flux < r_psf_flux_min_for_cadence1) |
                 (ps.i_stk_psf_flux < i_psf_flux_min_for_cadence1), cadence1),
                ((ps.g_stk_psf_flux < g_psf_flux_min_for_cadence2) |
                 (ps.r_stk_psf_flux < r_psf_flux_min_for_cadence2) |
                 (ps.i_stk_psf_flux < i_psf_flux_min_for_cadence2), cadence2),
                ((ps.g_stk_psf_flux >= g_psf_flux_min_for_cadence2) &
                 (ps.r_stk_psf_flux >= r_psf_flux_min_for_cadence2) &
                 (ps.i_stk_psf_flux >= i_psf_flux_min_for_cadence2), cadence3),
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

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        query = (
            c.select(
                fn.min(c.catalogid).alias('catalogid'),
                fn.min(ps.catid_objid).alias('ps1_catid_objid'),
                fn.min(tic.ps1_int).alias('gaia_source'),
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
                fn.min(tic.gaiamag).alias('ps1_g'),
                fn.min(tic.gaiabp).alias('bp'),
                fn.min(tic.gaiarp).alias('rp'),
                fn.min(opt_prov).alias('opt_prov'),
            )
            .join(c2ps)
            .join(ps)
            .join(x, on=(ps.catid_objid == x.ps1_dr2_id))
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
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(  # rely on catalogdb.catalog cross-matches to keep processing time down
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.specobjid) &
                    (s16.snmedian >= spec_sn_thresh) &
                    (s16.zwarning == 0) &
                    (s16.zerr <= spec_z_err_thresh) &
                    (s16.zerr > 0.0) &
                    (s16.scienceprimary > 0)
                )
            )
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(s2020.plug_ra, s2020.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (s2020.sn_median_all >= spec_sn_thresh) &
                    (s2020.zwarning == 0) &
                    (s2020.z_err <= spec_z_err_thresh) &
                    (s2020.z_err > 0.0)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.plug_ra, sV.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (sV.sn_median_all >= spec_sn_thresh) &
                    (sV.zwarning == 0) &
                    (sV.z_err <= spec_z_err_thresh) &
                    (sV.z_err > 0.0)
                )
            )
            .where(
                c.version_id == version_id,
                c2ps.version_id == version_id,
                c2ps.best >> True
            )
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
                (ps.g_stk_psf_flux != 'NaN'),   #  TODO check this is correct test via peewee
                (ps.r_stk_psf_flux != 'NaN'),
                (ps.i_stk_psf_flux != 'NaN'),
                (ps.flags.bin_and(good_stack_flag) > 0),
            )
            .group_by(ps.catid_objid)   # avoid duplicates - we trust the ps1 ids
            .having(
                # each and every match to the tic must satisfy the bright star rejection criteria
                fn.sum(
                    peewee.Case
                    (None,
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


#
#
#
#
#
# # #This provides the following BHM cartons:
# #    bhm_spiders_agn-efeds
# # and will eventually provide:
# #    bhm_spiders_agn-wide
# #    bhm_spiders_agn-deep
# # maybe:
# #    erosita-pointlike-bright-boss
#
#
#
#
#
# class BhmSpidersAgnEfedsCarton(BaseCarton):
#
#     '''
#     SELECT * from bhm_spiders_agn_superset AS x
#     INNER JOIN legacy_survey_dr8 AS ls ON  x.ls_id = ls.ls_id
#     LEFT JOIN bhm_efeds_veto AS v ON q3c_join(x.opt_ra,x.opt_dec,v.plug_ra,v.plug_dec,1.0)
#     WHERE x.ero_version = "efeds_c940_V2T"
#     AND WHERE x.ero_det_like > X.X
#     AND WHERE (ls.fiberflux_r > D.D OR
#                ls.fiberflux_z > E.E)                           # faint limits
#     AND WHERE ls.fibertotflux_r < CCC.C                        # bright limit - use total flux
#     AND WHERE (v.plug_ra = NULL OR
#                v.sn_median_all < 1.x OR
#                v.zwarning > 0 OR
#                v.z_err > 2e-3 OR
#                v.z_err <= 0.0)
#     '''
#
#     category = 'science'
#     mapper = 'BHM'
#     program = 'bhm_spiders'
#     tile = False
#
#
#     name = 'bhm_spiders_agn-efeds'
#     cadence = 'bhm_spiders_1x8'
#
#
#
#     def build_query(self, version_id, query_region=None):
#
#         c = Catalog.alias()
#         x = BHM_Spiders_AGN_Superset.alias()
#         ls = Legacy_Survey_DR8.alias()
#         c2ls = CatalogToLegacy_Survey_DR8.alias()
#         v = BHM_eFEDS_Veto.alias()
# #TODO        c2s = CatalogToSDSS_DR16_SpecObj.alias()
#         s = SDSS_DR16_SpecObj.alias()
#
#         fiberflux_r_max = AB2nMgy(self.parameters['fibermag_r_min'])
#         fiberflux_r_min = AB2nMgy(self.parameters['fibermag_r_max'])
#         fiberflux_z_min = AB2nMgy(self.parameters['fibermag_z_max'])
#
#         flux30 = AB2nMgy(30.00)
#
#         value = peewee.Value(self.parameters.get('value', 1.0)).cast('float').alias('value')
#         match_radius_spectro = self.parameters['spec_join_radius']/3600.0
#
#         p_f = self.parameters['priority_floor']
#         priority = peewee.Case(None,
#                                (
#                                    ((x.ero_det_like < self.parameters['det_like_for_priority']) & (s.specobjid.is_null(True)), p_f+6),  # noqa: E501
#                                    ((x.ero_det_like < self.parameters['det_like_for_priority']) & (s.specobjid.is_null(False)), p_f+7),  # noqa: E501
#                                    ((x.xmatch_flags == 1) & (s.specobjid.is_null(True)), p_f+0),  # noqa: E501
#                                    ((x.xmatch_flags == 0) & (s.specobjid.is_null(True)), p_f+1),  # noqa: E501
#                                    ((x.xmatch_flags > 1) & (s.specobjid.is_null(True)), p_f+2),  # noqa: E501
#                                    ((x.xmatch_flags == 1) & (s.specobjid.is_null(False)), p_f+3),  # noqa: E501
#                                    ((x.xmatch_flags == 0) & (s.specobjid.is_null(False)), p_f+4),  # noqa: E501
#                                    ((x.xmatch_flags > 1) & (s.specobjid.is_null(False)), p_f+5),  # noqa: E501
#                                ),
#                                p_f+9) ## should never get here
#
#
#         # Notes on convertion from ls_fibermag to sdss_fiber2mag:
#         # https://wiki.mpe.mpg.de/eRosita/EroAGN_eFEDS/SDSSIVSpecialPlates#Estimating_SDSS_fiber2mag_.2A_from_legacysurvey_photometry # noqa: E501
#
#         # Notes on converting from sdss_fiber2mag to sdss_psfmag
#         # https://wiki.sdss.org/display/OPS/Contents+of+targetdb.magnitude#Contentsoftargetdb.magnitude-WhatmagnitudestoputintotheplPlugMapfilesforBOSSplatetargets? # noqa: E501
#
#         # A flux ratio of 0.6 (roughly what is seen in all three ls bands)
#         # is a magnitude difference of
#         # fiber2mag(SDSS)-fibermag(LS) = 0.55mags
#         flux_ratio = {'g' : 0.60, 'r' : 0.60, 'i' : 0.60, 'z' : 0.60 }
#         # Then, also add the correction from sdss_fiber2mag to sdss_psfmag: psfmag_minus_fiber2mag('filter')  # noqa: E501
#
#         # legacysurvey fibermags - derived from fiberfluxes - with limits to avoid divide by zero errors  # noqa: E501
#         magnitude_g = (psfmag_minus_fiber2mag('g') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_g*flux_ratio['g']))).cast('float')
#         magnitude_r = (psfmag_minus_fiber2mag('r') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_r*flux_ratio['r']))).cast('float')
#         magnitude_z = (psfmag_minus_fiber2mag('z') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_z*flux_ratio['z']))).cast('float')
#         # the simplest possible interpolation
#         # TODO - we can do this better
#         magnitude_i = (psfmag_minus_fiber2mag('i') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,
#                                                      0.5*(ls.fiberflux_r+
#                                                           ls.fiberflux_z)*flux_ratio['i']))).cast('float')
#
#
#         query = (
#             c
#             .select(c.catalogid,
#                     priority.alias('priority'),
#                     value,
#                     magnitude_g.alias("g"),
#                     magnitude_r.alias("r"),
#                     magnitude_i.alias("i"),
#                     magnitude_z.alias("z"),
#             )
#             .join(c2ls)
#             .join(ls)
#             .join(x)
#             .join(v, JOIN.LEFT_OUTER,
#                   on=(fn.q3c_join(c.ra,c.dec,
#                                  v.plug_ra,v.plug_dec,
#                                  match_radius_spectro) &
#                       (v.sn_median_all >= self.parameters['spec_sn_thresh']) &
#                       (v.zwarning == 0) &
#                       (v.z_err <= self.parameters['spec_z_err_thresh']) &
#                       (v.z_err > 0.0)
#                       )
#                   )
#             .join(s, JOIN.LEFT_OUTER,
#                   on=(fn.q3c_join(c.ra,c.dec,
#                                   s.ra,s.dec,
#                                   match_radius_spectro) &
#                       (s.snmedian >= self.parameters['spec_sn_thresh']) &
#                       (s.zwarning == 0) &
#                       (s.zerr <= self.parameters['spec_z_err_thresh']) &
#                       (s.zerr > 0.0) &
#                       (s.scienceprimary > 0)
#                   )
#             )
#             .where(c.version_id == version_id,
#                    c2ls.version_id == version_id,
#                    c2ls.best is True)
#             .where(
#                 (x.ero_version == self.parameters['ero_version'] ),
#                 (v.pk.is_null()),
#                 (ls.fibertotflux_r < fiberflux_r_max),
#                 (
#                     (ls.fiberflux_r >= fiberflux_r_min) |
#                     (ls.fiberflux_z >= fiberflux_z_min)
#                 ),
#                 (x.ero_det_like > self.parameters['det_like_min']),
#                 (
#                     (
#                         (x.xmatch_method == 'XPS-ML/NWAY') &
#                         (x.xmatch_metric >= self.parameters['p_any_min'])
#                     ) |
#                     (
#                         (x.xmatch_method == 'XPS-LR') &
#                         (x.xmatch_metric >= self.parameters['lr_min'])
#                     )
#                 )
#             )
#             .distinct([ls.ls_id])   # avoid duplicates - trust the ls_id
#         )
#
#         return query
#
#
#
# ### some useful snippets:
# '''
# # example to get the listing of fields from a PeeWee model
# print(catalogdb.ErositaAGNMock._meta.fields)
#
# x = BHM_Spiders_AGN_Superset.alias()
# ls = Legacy_Survey_DR8.alias()
# #q14 = SDSS_DR14_QSO.alias()
# c = Catalog.alias()
#
#
# for f in c._meta.fields:
#     print (f)
# for f in x._meta.fields:
#     print (f)
# for f in ls._meta.fields:
#     print (f)
# #for f in q14._meta.fields:
# #    print (f)
#
# '''
# ###

# deferred#
# deferred#
# deferred#
# deferred#
# deferred# class BhmSpidersBaseCarton(BaseCarton):
# deferred#     ''' Parent class that provides the mask selections for any SPIDERs catalogue'''
# deferred#
# deferred#     name = 'bhm_spiders'
# deferred#     category = 'science'
# deferred#     mapper = 'BHM'
# deferred#     program = 'BHM-SPIDERS'
# deferred#     tile = False
# deferred#
# deferred#     # list of skymasks - move to the config file?
# deferred#     skymasks = [
# deferred#     ]
# deferred#
# deferred#
# deferred#     def post_process(self, model, **kwargs):
# deferred#         # this is where we select on mask location
# deferred#         # get the ids and coords of all of the objects in the temp table
# deferred#         cat_id = np.array(model.catalog_id[:])
# deferred#         ra = np.array(model.ra[:])
# deferred#         dec = np.array(model.dec[:])
# deferred#
# deferred#         flags = np.ones(len(cat_id), np.bool)
# deferred#
# deferred#         for sm in self.skymasks:
# deferred#             sm.apply(lon=ra, lat=dec, flags=flags)
# deferred#
# deferred#         # not sure what to return at this point - a list of tuples?
# deferred#         result = [(i,flag) for i,flag in zip(cat_id,flags)]
# deferred#
# deferred#         return result
# deferred#
# deferred#
# deferred#
# deferred# class BhmSpidersWideBaseCarton(BhmSpidersBaseCarton):
# deferred#     ''' Parent class that provides the mask slections for any SPIDER-wide catalogue'''
# deferred#
# deferred#     name = 'bhm_spiders_wide'
# deferred#     category = 'science'
# deferred#     mapper = 'BHM'
# deferred#     program = 'BHM-SPIDERS'
# deferred#     tile = False
# deferred#
# deferred#     # list of skymasks
# deferred#     skymasks = [
# deferred#         SkyMask(filename=pkg_resources.resource_filename(
# deferred#             __name__,
# deferred#             'masks/eROSITA-DE_exgal_lsdr8_or_psdr2_proc.ply'),
# deferred#                 name="spiders_wide",
# deferred#                 masktype="mangle",
# deferred#                 sense="include",
# deferred#         ),
# deferred#         SkyMask(filename=pkg_resources.resource_filename(
# deferred#             __name__,
# deferred#             'masks/rsFields-annotated-lco-deep_proc.ply'),
# deferred#                 name="spiders_deep",
# deferred#                 masktype="mangle",
# deferred#                 sense="exclude",
# deferred#         ),
# deferred#     ]
# deferred#
# deferred#
# deferred# class BhmSpidersAgnWideLsCarton(BhmSpidersWideBaseCarton):
# deferred#
# deferred#     '''
# deferred#     spiders_agn_wide_ls:
# deferred#
# deferred#     SELECT * from bhm_spiders_agn_superset AS x
# deferred#     INNER JOIN legacy_survey_dr8 AS ls ON x.ls_id = ls.ls_id
# deferred#     WHERE x.ero_version = "version_code_TBD"
# deferred#     AND WHERE x.ero_det_like > X.X
# deferred#     AND WHERE x.xmatch_metric > 0.x
# deferred#     AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
# deferred#     AND WHERE ls.fibertotflux_r < CCC.C
# deferred#     '''
# deferred#
# deferred#     name = 'bhm_spiders_agn_wide_ls'
# deferred#     cadence = 'bhm_spiders_1x4'
# deferred#
# deferred#     def build_query(self, version_id, query_region=None):
# deferred#         c = Catalog.alias()
# deferred#         x = BHM_Spiders_AGN_Superset.alias()
# deferred#         ls = Legacy_Survey_DR8.alias()
# deferred#         c2ls = CatalogToLegacy_Survey_DR8.alias()
# deferred#
# deferred#         flux_r_max =  AB2nMgy(self.parameters['mag_r_min'])
# deferred#         flux_r_min =  AB2nMgy(self.parameters['mag_r_max'])
# deferred#         flux_z_min =  AB2nMgy(self.parameters['mag_z_max'])
# deferred#         target_value = peewee.Value(self.parameters.get('value', 1.0)).alias('value')
# deferred#
# deferred#         query = (
# deferred#             c
# deferred#             .select(c.catalogid,
# deferred#                     c.ra,
# deferred#                     c.dec,
# deferred#                     c.pmra,
# deferred#                     c.pmdec,
# deferred#                     x.target_priority.alias('priority'),
# deferred#                     ls.fiberflux_g.alias('lsfiberflux_g'),
# deferred#                     ls.fiberflux_r.alias('lsfiberflux_r'),
# deferred#                     ls.fiberflux_z.alias('lsfiberflux_z'),
# deferred#             )
# deferred#             .join(c2ls)
# deferred#             .join(ls)
# deferred#             .join(x)
# deferred#             .where(c.version_id == version_id,
# deferred#                    c2ls.version_id == version_id)
# deferred#             .where(
# deferred#                 (x.ero_version == self.parameters['ero_version'] ),
# deferred#                 (ls.fibertotflux_r < flux_r_max),
# deferred#                 ((ls.fiberflux_r   > flux_r_min) |
# deferred#                  (ls.fiberflux_z > flux_z_min) ),
# deferred#                 (x.ero_det_like > self.parameters['det_like_min']),
# deferred#                 (x.xmatch_metric > self.parameters['p_any_min']),
# deferred#             )
# deferred#         )
# deferred#
# deferred#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")  # noqa: E501
# deferred#
# deferred#         return query
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred# #wait_for_psdr2# class BhmSpidersAgnWidePsCarton(BhmSpidersWideBaseCarton):
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#     '''
# deferred# #wait_for_psdr2#     spiders_agn_wide_ps:
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#     SELECT * from bhm_spiders_agn_superset AS x
# deferred# #wait_for_psdr2#     INNER JOIN panstarrs_dr2 AS ps ON x.ps1_dr2_objid = ls.ls_id
# deferred# #wait_for_psdr2#     WHERE x.ero_version = "version_code_TBD"
# deferred# #wait_for_psdr2#     AND WHERE x.ero_det_like > X.X
# deferred# #wait_for_psdr2#     AND WHERE x.xmatch_metric > 0.x
# deferred# #wait_for_psdr2#     AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
# deferred# #wait_for_psdr2#     AND WHERE ls.fiberflux_r < CCC.C
# deferred# #wait_for_psdr2#     '''
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#     name = 'bhm_spiders_agn_wide_ls'
# deferred# #wait_for_psdr2#     cadence = 'bhm_spiders_1x4'
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#     def build_query(self, version_id, query_region=None):
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#         c = Catalog.alias()
# deferred# #wait_for_psdr2#         x = BHM_Spiders_AGN_Superset.alias()
# deferred# #wait_for_psdr2#         ps = PanStarrsDr2.alias()
# deferred# #wait_for_psdr2#         c2ps = CatalogToPanStarrsDr2.alias()
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#         flux_r_max =  AB2Jy(self.parameters['r_mag_min'])
# deferred# #wait_for_psdr2#         flux_r_min =  AB2Jy(self.parameters['r_mag_max'])
# deferred# #wait_for_psdr2#         flux_z_min =  AB2Jy(self.parameters['z_mag_max'])
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#         query = (
# deferred# #wait_for_psdr2#             c
# deferred# #wait_for_psdr2#             .select(c.catalogid,
# deferred# #wait_for_psdr2#                     c.ra,
# deferred# #wait_for_psdr2#                     c.dec,
# deferred# #wait_for_psdr2#                     c.pmra,
# deferred# #wait_for_psdr2#                     c.pmdec,
# deferred# #wait_for_psdr2#                     x.target_priority.alias('priority'),
# deferred# #wait_for_psdr2#                     ps.g_stk_aper_flux.alias('psaperflux_g'),
# deferred# #wait_for_psdr2#                     ps.r_stk_aper_flux.alias('psaperflux_r'),
# deferred# #wait_for_psdr2#                     ps.i_stk_aper_flux.alias('psaperflux_i'),
# deferred# #wait_for_psdr2#                     ps.z_stk_aper_flux.alias('psaperflux_z'))
# deferred# #wait_for_psdr2#             .join(ps)
# deferred# #wait_for_psdr2#             .join(c2ps)
# deferred# #wait_for_psdr2#             .where(
# deferred# #wait_for_psdr2#                 (ps.r_stk_aper_flux < flux_r_max) &
# deferred# #wait_for_psdr2#                 ((ps.r_stk_aper_flux > flux_r_min) |
# deferred# #wait_for_psdr2#                  (ps.z_stk_aper_flux > flux_z_min) ) &
# deferred# #wait_for_psdr2#                 (tab.ero_det_like > self.parameters['det_like_min']) &
# deferred# #wait_for_psdr2#                 (tab.xmatch_metric > self.parameters['p_any_min'])
# deferred# #wait_for_psdr2#             )
# deferred# #wait_for_psdr2#         )
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#         print(f"This query will return nrows={query.count()} "
# deferred# #wait_for_psdr2#             f"(c.f. req_ntargets={self.parameters['req_ntargets']})")
# deferred# #wait_for_psdr2#
# deferred# #wait_for_psdr2#         return query
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred#
# deferred# #
# deferred# # ##############################################################
# deferred# # ##############################################################
# deferred# # ##############################################################
# deferred# # ## MOCK TARGETS ##############################################
# deferred# # class BhmSpidersAgnWideMockCarton(BaseCarton):
# deferred# #
# deferred# #     name = 'bhm_spiders_agn_wide_mock'
# deferred# #     category = 'science'
# deferred# #     survey = 'BHM'
# deferred# #     cadence = 'bhm_spiders_1x4'
# deferred# #     tile = False
# deferred# #
# deferred# #     # list of masks - remove to config file?
# deferred# #     masks = [
# deferred# #         {"name": "wide",
# deferred# #          "type": "mangle",
# deferred# #          "polarity": "include",
# deferred# #          "filename": "eROSITA-DE_exgal_lsdr8_or_psdr2_proc.ply",
# deferred# #         },
# deferred# #         {"name": "deep",
# deferred# #          "type": "mangle",
# deferred# #          "polarity": "exclude",
# deferred# #          "filename" : "rsFields-annotated-lco-deep_proc.ply",
# deferred# #         },
# deferred# #     ]
# deferred# #
# deferred# #     def build_query(self, version_id, query_region=None):
# deferred# #         '''
# deferred# #         Pure database level query - generates a super-set of potential targets
# deferred# #         '''
# deferred# #         # get the table name from the config
# deferred# #         #   - maybe replace this with a list of options
# deferred# #         exec(f'tab = catalogdb.{params["catalogdb_table"]}')
# deferred# #         assert tab is not None, 'Failed to locate catalogdb_table'
# deferred# #
# deferred# #         query = (tab.select(tab.gaia_dr2_source_id.alias('catalog_id'),
# deferred# #                             tab.target_ra.alias('ra'),
# deferred# #                             tab.target_dec.alias('dec'),
# deferred# #                             tab.target_pmra.alias('pmra'),
# deferred# #                             tab.target_pmdec.alias('pmdec'),
# deferred# #                             tab.target_epoch.alias('epoch'),
# deferred# #                             tab.target_mag_r.alias('magnitude_g'),
# deferred# #                             tab.target_mag_r.alias('magnitude_r'),
# deferred# #                             tab.target_mag_r.alias('magnitude_z'))
# deferred# #                  .where((tab.target_mag_r > self.parameters['r_mag_min']) &
# deferred# #                         (tab.target_mag_r < self.parameters['r_mag_max']) &
# deferred# #                         (tab.ero_det_like_0 > self.parameters['det_like_0_min'])))
# deferred# #
# deferred# #         print(f"This query will return nrows={query.count()} "
# deferred# #                " (c.f. req_ntargets={self.parameters['req_ntargets']})")
# deferred# #
# deferred# #         return query
# deferred# #
# deferred# #
# deferred# #     def post_process(self, model, **kwargs):
# deferred# #         # this is where we select on mask location
# deferred# #         # get the coords of all of the objects oin the temp table
# deferred# #         ra = np.array(model.ra[:])
# deferred# #         dec = np.array(model.dec[:])
# deferred# #
# deferred# #         return True
# deferred# #
# deferred# # ## MOCK TARGETS ##############################################
# deferred# # ##############################################################
# deferred# # ##############################################################
# deferred# # ##############################################################
# deferred#
# deferred#
# deferred#
# deferred#
# deferred# def _test_xmatch_stuff():
# deferred#     # check numbers of targets in a test patch
# deferred#
# deferred#     version_id = 11
# deferred#     search_radius_deg = 0.01
# deferred#     ra0 = 135.0
# deferred#     dec0 = 1.0
# deferred#
# deferred#     c = Catalog.alias()
# deferred# #    x = BHM_Spiders_AGN_Superset.alias()
# deferred#     ls = Legacy_Survey_DR8.alias()
# deferred#     c2ls = CatalogToLegacy_Survey_DR8.alias()
# deferred#
# deferred#
# deferred#     query = (
# deferred#         c
# deferred#         .select(c.catalogid,
# deferred#                 c.ra,
# deferred#                 c.dec,
# deferred#                 c.lead,
# deferred#                 c.version,
# deferred#                 ls.ls_id,
# deferred#                 ls.ra,
# deferred#                 ls.dec,
# deferred#         )
# deferred#         .join(c2ls)
# deferred#         .join(ls)
# deferred#         .where(c.version_id == version_id,
# deferred#                c2ls.version_id == version_id)
# deferred#         .where(
# deferred#             peewee.fn.q3c_radial_query(c.ra,c.dec,
# deferred#                                        ra0, dec0,
# deferred#                                        search_radius_deg)
# deferred#         )
# deferred#     )
# deferred#
# deferred#
# deferred#     query.select().limit(1000).count()
# deferred#


# '''
# Exporting from the temp table
#
# \copy (SELECT * FROM sandbox.temp_bhm_spiders_agn_efeds)
#  TO 'bhm_spiders_agn_efeds.csv' with csv header
# stilts tpipe in=bhm_spiders_agn_efeds.csv out=bhm_spiders_agn_efeds.fits ifmt=csv ofmt=fits-basic
#
#
# '''
#
