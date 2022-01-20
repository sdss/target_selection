#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# ## flake8: noqa
# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn
from sdssdb.peewee.sdss5db.catalogdb import (
    # Catalog,
    BHM_CSC_v2,
    Panstarrs1,
    CatalogToPanstarrs1,
    CatalogToTIC_v8,
    TIC_v8
)

# imports of existing spectro catalogues
# from sdssdb.peewee.sdss5db.catalogdb import (
#     CatalogToSDSS_DR16_SpecObj,
#     SDSS_DR16_SpecObj,
#     CatalogToBHM_eFEDS_Veto,
#     BHM_eFEDS_Veto,
#     SDSSV_BOSS_SPALL,
#     SDSSV_BOSS_Conflist,
#     SDSSV_Plateholes,
#     SDSSV_Plateholes_Meta,
# )

from target_selection.mag_flux import AB2Jy

from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ChandraSourceCatalogue(CSC)-inprogress  # noqa: E501

# This module provides the following BHM cartons:
#   bhm_csc_boss_dark
#   bhm_csc_boss_bright
#   bhm_csc_apogee


class BhmCscBossBrightCarton(BaseCarton):
    '''
    SELECT x.cxoid,x.ocat,
           COALESCE(c2t.catalogid,c2p.catalogid) as catalogid,
           x.ora,x.odec,
           x.omag,
           8.9 - 2.5*log10(p.i_stk_psf_flux) as i_mag,
           p.i_chp_psf,
           tic.gaiamag
        FROM catalogdb.bhm_csc_v2 as x
        LEFT OUTER JOIN panstarrs1 as p
            on x.idps = p.extid_hi_lo
        LEFT OUTER JOIN tic_v8 as tic
            on tic.gaia_int = x.idg2
        LEFT OUTER JOIN catalog_to_panstarrs1 as c2p
            on p.catid_objid = c2p.target_id
        LEFT OUTER JOIN catalog_to_tic_v8 as c2t
            on c2t.target_id = tic.id
    WHERE (   (     x.ocat = 'P'
                AND p.i_stk_psf_flux BETWEEN 0.000575 AND 0.01445
                AND p.i_stk_psf_flux != 'NaN')
           OR (x.ocat = 'G'
               AND (NOT tic.gaiamag IS NULL)
               AND tic.gaiamag BETWEEN 13. AND 17.0)
           )
           AND ( COALESCE(c2p.version_id,c2t.version_id) = 25 )
    LIMIT 20;

    '''
    name = 'bhm_csc_boss_bright'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_csc'
    instrument = 'BOSS'
    tile = False

    def build_query(self, version_id, query_region=None):
        x = BHM_CSC_v2.alias()
        ps = Panstarrs1.alias()
        c2ps = CatalogToPanstarrs1.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        g_psf_flux_max = AB2Jy(self.parameters['g_psf_mag_min'])
        r_psf_flux_max = AB2Jy(self.parameters['r_psf_mag_min'])
        i_psf_flux_max = AB2Jy(self.parameters['i_psf_mag_min'])
        z_psf_flux_max = AB2Jy(self.parameters['z_psf_mag_min'])
        g_psf_flux_min = AB2Jy(self.parameters['g_psf_mag_max'])
        r_psf_flux_min = AB2Jy(self.parameters['r_psf_mag_max'])
        i_psf_flux_min = AB2Jy(self.parameters['i_psf_mag_max'])
        z_psf_flux_min = AB2Jy(self.parameters['z_psf_mag_max'])

        gaia_g_mag_min = self.parameters['gaia_g_mag_min']
        gaia_g_mag_max = self.parameters['gaia_g_mag_max']
        gaia_rp_mag_min = self.parameters['gaia_rp_mag_min']
        gaia_rp_mag_max = self.parameters['gaia_rp_mag_max']

        value = peewee.Value(self.parameters['value'])
        priority = peewee.Value(self.parameters['priority_floor'])

        query = (
            x.select(
                fn.coalesce(c2ps.catalogid, c2tic.catalogid).alias('catalogid'),
                x.cxoid.alias('csc_cxoid'),
                x.ora.alias('csc_ora'),
                x.odec.alias('csc_odec'),
                x.omag.alias('csc_omag'),
                x.ocat.alias('csc_ocat'),
                x.idps.alias('csc_idps'),
                x.idg2.alias('csc_idg2'),
                priority.alias('priority'),
                value.alias('value'),
            )
            .join(ps, join=JOIN.LEFT_OUTER,
                  on=(x.idps == ps.extid_hi_lo))
            .join(tic, join=JOIN.LEFT_OUTER,
                  on=(x.idg2 == tic.gaia_int))
            .join(c2ps, join=JOIN.LEFT_OUTER,
                  on=(ps.catid_objid == c2ps.target_id))
            .join(c2tic, join=JOIN.LEFT_OUTER,
                  on=(tic.id == c2tic.target_id))
            .where(
                fn.coalesce(c2ps.version_id, c2tic.version_id) == version_id,
                fn.coalesce(c2ps.best, c2tic.best) >> True,
            )
            .where(
                (
                    x.ocat == 'P',
                    ps.g_stk_psf_flux < g_psf_flux_max,
                    ps.r_stk_psf_flux < r_psf_flux_max,
                    ps.i_stk_psf_flux < i_psf_flux_max,
                    ps.z_stk_psf_flux < z_psf_flux_max,
                    (
                        (ps.g_stk_psf_flux > g_psf_flux_min & ps.g_stk_psf_flux != 'NaN') |
                        (ps.r_stk_psf_flux > r_psf_flux_min & ps.r_stk_psf_flux != 'NaN') |
                        (ps.i_stk_psf_flux > i_psf_flux_min & ps.i_stk_psf_flux != 'NaN') |
                        (ps.z_stk_psf_flux > z_psf_flux_min & ps.z_stk_psf_flux != 'NaN')
                    ),
                ) |
                (
                    x.ocat == 'G',
                    tic.gaiamag > gaia_g_mag_min,
                    tic.gaiarp > gaia_rp_mag_min,
                    ((tic.gaiamag < gaia_g_mag_max) |
                     (tic.gaiarp < gaia_rp_mag_max)),
                )
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(x.ora, x.odec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query



# class BhmCscBaseCarton(BaseCarton):
#     '''
#     Parent class that provides the basic selections for all CSC cartons
#     To be sub-classed, not typically to be called directly.
#
#     SELECT DISTINCT ON (t.pk)
#       c.catalogid, c.lead, c.ra, c.dec,
#       q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
#       t.mag_g, t.mag_r, t.mag_i
#     FROM catalog AS c
#     JOIN catalog_to_bhm_csc AS c2t
#       ON c.catalogid = c2t.catalogid
#     JOIN bhm_csc AS t
#       ON c2t.target_id = t.pk
#     WHERE
#           c.version_id = 21
#       AND c2t.version_id = 21
#       AND c2t.best IS TRUE
#       AND t.spectrograph = 'BOSS'
#       AND t.mag_i BETWEEN 16.0 AND 22.0
#     ORDER BY t.pk
#     LIMIT 20;
#
#
#     # check how many csc targets have positions defined by the  bhm_csc.oir_ra,oir_dec fields
#     # these may actually be Chandra X-ray coords and not optical positions!
#
#     SELECT
#         count(*) as count,
#         AVG(sep) as mean_sep,
#         MIN(sep) as min_sep,
#         MAX(sep) as max_sep,
#         AVG(mag_r) as mean_mag_r,
#         MIN(mag_r) as min_mag_r,
#         MAX(mag_r) as max_mag_r
#     FROM (
#         SELECT DISTINCT ON (t.pk)
#            c.catalogid, c.lead, c.ra, c.dec,
#           q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
#           t.mag_g, t.mag_r, t.mag_i
#         FROM catalog AS c
#         JOIN catalog_to_bhm_csc AS c2t
#           ON c.catalogid = c2t.catalogid
#         JOIN bhm_csc AS t
#           ON c2t.target_id = t.pk
#         WHERE
#               c.version_id = 21
#           AND c2t.version_id = 21
#           AND c2t.best IS TRUE
#           AND t.spectrograph = 'BOSS'
#           AND t.mag_i BETWEEN 16.0 AND 22.0
#         ORDER BY t.pk
#         ) AS subq
#       WHERE subq.lead = 'bhm_csc';
#
#     #
#     SELECT DISTINCT ON (t.pk)
#           c.catalogid, c.lead, c.ra, c.dec,
#           q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
#           t.mag_g, t.mag_r, t.mag_i
#         INTO sandbox.temp_td_csc_subset
#         FROM catalog AS c
#         JOIN catalog_to_bhm_csc AS c2t
#           ON c.catalogid = c2t.catalogid
#         JOIN bhm_csc AS t
#           ON c2t.target_id = t.pk
#         WHERE
#               c.version_id = 21
#           AND c2t.version_id = 21
#           AND c2t.best IS TRUE
#           AND t.spectrograph = 'BOSS'
#           AND t.mag_i BETWEEN 16.0 AND 22.0
#           AND c.lead = 'bhm_csc'
#         ORDER BY t.pk;
#
#
#
#     '''
#
#     name = 'bhm_csc_base'
#     base_name = 'bhm_csc_base'
#     category = 'science'
#     mapper = 'BHM'
#     program = 'bhm_csc'
#     tile = False
#     priority = None
#     alias_t = None
#     instrument = None
#     this_cadence = None
#
#     def build_query(self, version_id, query_region=None):
#         c = Catalog.alias()
#         c2t = CatalogToBHM_CSC.alias()
#         t = BHM_CSC.alias()
#         self.alias_t = t
#         # c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
#         # s16 = SDSS_DR16_SpecObj.alias()
#         # s2020 = BHM_eFEDS_Veto.alias()
#         # sV = SDSSV_BOSS_SPALL.alias()
#         # ph = SDSSV_Plateholes.alias()
#         # phm = SDSSV_Plateholes_Meta.alias()
#
#         # set the Carton priority+values here - read from yaml
#         value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
#         instrument = peewee.Value(self.instrument)
#         cadence = peewee.Value(self.this_cadence)
#         # opt_prov = peewee.Value('ps1_psfmag')
#
#         if (self.instrument == 'BOSS'):
#
#             # #########################################################################
#             # prepare the spectroscopy catalogues
#             match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
#             spec_sn_thresh = self.parameters['spec_sn_thresh']
#             spec_z_err_thresh = self.parameters['spec_z_err_thresh']
#             dpriority_has_spec = self.parameters['dpriority_has_spec']
#
#             # SDSS DR16
#             c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
#             ss16 = SDSS_DR16_SpecObj.alias()
#             s16 = (
#                 ss16.select(
#                     ss16.specobjid.alias('specobjid'),
#                 )
#                 .where(
#                     ss16.snmedian >= spec_sn_thresh,
#                     ss16.zwarning == 0,
#                     ss16.zerr <= spec_z_err_thresh,
#                     ss16.zerr > 0.0,
#                     ss16.scienceprimary > 0,
#                 )
#                 .alias('s16')
#             )
#
#             # SDSS-IV/eFEDS March2020
#             c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
#             ss2020 = BHM_eFEDS_Veto.alias()
#             s2020 = (
#                 ss2020.select(
#                     ss2020.pk.alias('pk'),
#                 )
#                 .where(
#                     ss2020.sn_median_all >= spec_sn_thresh,
#                     ss2020.zwarning == 0,
#                     ss2020.z_err <= spec_z_err_thresh,
#                     ss2020.z_err > 0.0,
#                 )
#                 .alias('s2020')
#             )
#
#             # SDSS-V spAll
#             ssV = SDSSV_BOSS_SPALL.alias()
#             sV = (
#                 ssV.select(
#                     ssV.specobjid.alias('specobjid'),
#                     ssV.plug_ra.alias('plug_ra'),
#                     ssV.plug_dec.alias('plug_dec'),
#                 )
#                 .where(
#                     ssV.sn_median_all >= spec_sn_thresh,
#                     ssV.zwarning == 0,
#                     ssV.z_err <= spec_z_err_thresh,
#                     ssV.z_err > 0.0,
#                     ssV.specprimary > 0,
#                 )
#                 .alias('sV')
#             )
#
#             # SDSS-V plateholes - only consider plateholes that
#             # were drilled+shipped but that were not yet observed
#             ssph = SDSSV_Plateholes.alias()
#             ssphm = SDSSV_Plateholes_Meta.alias()
#             ssconf = SDSSV_BOSS_Conflist.alias()
#             sph = (
#                 ssph.select(
#                     ssph.pkey.alias('pkey'),
#                     ssph.target_ra.alias('target_ra'),
#                     ssph.target_dec.alias('target_dec'),
#                 )
#                 .join(
#                     ssphm,
#                     on=(ssph.yanny_uid == ssphm.yanny_uid)
#                 )
#                 .join(
#                     ssconf, JOIN.LEFT_OUTER,
#                     on=(ssphm.plateid == ssconf.plate)
#                 )
#                 .where(
#                     (ssph.holetype == 'BOSS_SHARED'),
#                     (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
#                     ssphm.isvalid > 0,
#                     ssconf.plate.is_null(),
#                 )
#                 .alias('sph')
#             )
#
#             # adjust priority if target aleady has an SDSS spectrum
#             priority_1 = peewee.Case(
#                 None,
#                 (
#                     (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
#                     (s2020.c.pk.is_null(False), 1),
#                     (sV.c.specobjid.is_null(False), 1),
#                     (sph.c.pkey.is_null(False), 1),
#                 ),
#                 0)
#             #
#             # Compute net priority
#             priority = (
#                 peewee.Value(self.parameters['priority_floor']) +
#                 priority_1 * dpriority_has_spec
#             )
#         else:
#             priority = peewee.Value(self.parameters['priority_floor'])
#
#         # compute transformed SDSS mags for pointlike and extended sources separately
#         # transform the csc (panstarrs1-dr1) griz into sdss psfmag griz
#
#         # extract coeffs from fit logs via:
#         # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_csc_boss/ts_mag_to_sdss_psfmag_?_results.log  # noqa
#
#         coeffs = {
#             "g2_p": 0.087878,
#             "g1_p": 0.063329,
#             "g0_p": 0.021488,
#             "i2_p": -0.011220,
#             "i1_p": 0.020782,
#             "i0_p": 0.000154,
#             "r2_p": -0.093371,
#             "r1_p": 0.136032,
#             "r0_p": -0.011477,
#             "z2_p": -0.180526,
#             "z1_p": 0.007284,
#             "z0_p": -0.037933,
#         }
#         # Note that the corrections for r and i are very small,
#         # however g+z both have non-negligible colour terms
#
#         g0 = peewee.Case(None, ((t.mag_g <= 0.0, None),), t.mag_g)
#         r0 = peewee.Case(None, ((t.mag_r <= 0.0, None),), t.mag_r)
#         i0 = peewee.Case(None, ((t.mag_i <= 0.0, None),), t.mag_i)
#         z0 = peewee.Case(None, ((t.mag_z <= 0.0, None),), t.mag_z)
#         g_r = g0 - r0
#         r_i = r0 - i0
#         i_z = i0 - z0
#
#         # use single set of transforms because we do not have any info in csc parent table to
#         # differentiate between pointlike and extended sources)
#         g = (g0 + coeffs['g0_p'] + coeffs['g1_p'] * g_r + coeffs['g2_p'] * g_r * g_r)
#         r = (r0 + coeffs['r0_p'] + coeffs['r1_p'] * g_r + coeffs['r2_p'] * g_r * g_r)
#         i = (i0 + coeffs['i0_p'] + coeffs['i1_p'] * r_i + coeffs['i2_p'] * r_i * r_i)
#         z = (z0 + coeffs['z0_p'] + coeffs['z1_p'] * i_z + coeffs['z2_p'] * i_z * i_z)
#
#         # validity checks (only griz) - set limits semi-manually
#         g_r_min = -0.3
#         g_r_max = 1.7
#         r_i_min = -0.5
#         r_i_max = 2.5
#         i_z_min = -0.3
#         i_z_max = 1.25
#         valid = (g0.between(0.1, 29.9) &
#                  r0.between(0.1, 29.9) &
#                  i0.between(0.1, 29.9) &
#                  z0.between(0.1, 29.9) &
#                  g_r.between(g_r_min, g_r_max) &
#                  r_i.between(r_i_min, r_i_max) &
#                  i_z.between(i_z_min, i_z_max))
#
#         opt_prov = peewee.Case(None, ((valid, 'sdss_psfmag_from_csc'),), 'undefined')
#         magnitude_g = peewee.Case(None, ((valid, g),), 'NaN')
#         magnitude_r = peewee.Case(None, ((valid, r),), 'NaN')
#         magnitude_i = peewee.Case(None, ((valid, i),), 'NaN')
#         magnitude_z = peewee.Case(None, ((valid, z),), 'NaN')
#         magnitude_h = peewee.Case(None, ((t.mag_h <= 0.0, None),), t.mag_h).cast('float')
#
#         # # Process the bhm_csc.[g,r,i,z,h] magnitudes to deal with zeros
#         # magnitude_g = peewee.Case(None, ((t.mag_g <= 0.0, None),), t.mag_g).cast('float')
#         # magnitude_r = peewee.Case(None, ((t.mag_r <= 0.0, None),), t.mag_r).cast('float')
#         # magnitude_i = peewee.Case(None, ((t.mag_i <= 0.0, None),), t.mag_i).cast('float')
#         # magnitude_z = peewee.Case(None, ((t.mag_z <= 0.0, None),), t.mag_z).cast('float')
#         # magnitude_h = peewee.Case(None, ((t.mag_h <= 0.0, None),), t.mag_h).cast('float')
#
#         # Create a subquery that will calculate the minimum catalog_to_bhm_csc.distance for each
#         # csc candidate target
#         subq = (
#             c2t
#             .select(
#                 c2t.target_id,
#                 fn.MIN(c2t.distance).alias('min_distance'))
#             .where(
#                 c2t.version_id == version_id,
#                 c2t.best >> True
#             )
#             .group_by(c2t.target_id)
#             .alias('min_dist_subq')
#         )
#
#         query = (
#             c.select(
#                 c.catalogid,
#                 t.cxo_name,   # extra
#                 t.pk.alias('csc_pk'),   # extra
#                 c.ra,  # extra
#                 c.dec,  # extra
#                 priority.alias('priority'),
#                 value.alias('value'),
#                 cadence.alias('cadence'),
#                 instrument.alias('instrument'),
#                 opt_prov.alias('optical_prov'),
#                 magnitude_g.alias('g'),
#                 magnitude_r.alias('r'),
#                 magnitude_i.alias('i'),
#                 magnitude_z.alias('z'),
#                 magnitude_h.alias('h'),
#                 t.mag_g.alias('csc_mag_g'),   # extra
#                 t.mag_r.alias('csc_mag_r'),   # extra
#                 t.mag_i.alias('csc_mag_i'),   # extra
#                 t.mag_z.alias('csc_mag_z'),   # extra
#                 t.oir_ra.alias('csc_ra'),   # extra
#                 t.oir_dec.alias('csc_dec'),   # extra
#             )
#             .join(c2t)
#             .join(t)
#             .join(
#                 subq,
#                 on=(
#                     (c2t.target_id == subq.c.target_id) &
#                     (
#                         (c2t.distance == subq.c.min_distance) |
#                         (c2t.distance.is_null() & subq.c.min_distance.is_null())
#                     )
#                 ),
#             )
#             .where(
#                 c.version_id == version_id,
#                 c2t.version_id == version_id,
#                 c2t.best >> True
#             )
#             # .distinct([c2t.target_id])  # avoid duplicates - trust the CSC parent sample,
#             # .distinct([c.catalogid])  # avoid duplicates - trust the catalogid,
#             # avoid duplicates - trust uniquness in both CSC name and catalogid
#             .distinct([c.catalogid])
#             # .distinct([t.cxo_name])
#             .where
#             (
#                 t.spectrograph == self.instrument
#             )
#         )
#
#         if (self.instrument == 'BOSS'):
#             # Append the spectro query
#             query = (
#                 query
#                 .switch(c)
#                 .join(c2s16, JOIN.LEFT_OUTER)
#                 .join(
#                     s16, JOIN.LEFT_OUTER,
#                     on=(
#                         (c2s16.target_id == s16.c.specobjid) &
#                         (c2s16.version_id == version_id)
#                     )
#                 )
#                 .switch(c)
#                 .join(c2s2020, JOIN.LEFT_OUTER)
#                 .join(
#                     s2020, JOIN.LEFT_OUTER,
#                     on=(
#                         (c2s2020.target_id == s2020.c.pk) &
#                         (c2s2020.version_id == version_id)
#                     )
#                 )
#                 .join(
#                     sV, JOIN.LEFT_OUTER,
#                     on=(
#                         fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
#                                     c.ra, c.dec,
#                                     match_radius_spectro)
#                     )
#                 )
#                 .join(
#                     sph, JOIN.LEFT_OUTER,
#                     on=(
#                         fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
#                                     c.ra, c.dec,
#                                     match_radius_spectro)
#                     )
#                 )
#             )
#
#         if query_region:
#             query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
#                                                            query_region[0],
#                                                            query_region[1],
#                                                            query_region[2]))
#
#         return query
#
#
# class BhmCscBossDarkCarton(BhmCscBaseCarton):
#     '''
#     SELECT * from bhm_csc AS c
#     WHERE c.spectrograph = "BOSS"
#        ND c.mag_i BETWEEN 17.x AND 21.x
#     '''
#     name = 'bhm_csc_boss_dark'
#     this_cadence = 'dark_1x4'
#     instrument = 'BOSS'
#
#     def build_query(self, version_id, query_region=None):
#         query = super().build_query(version_id, query_region)
#         t = self.alias_t
#         query = query.where(
#             (t.mag_i >= self.parameters['mag_i_min']) &
#             (t.mag_i < self.parameters['mag_i_max'])
#         )
#         return query
#
#
# class BhmCscBossBrightCarton(BhmCscBaseCarton):
#     '''
#     SELECT * from bhm_csc AS c
#     WHERE c.spectrograph = "BOSS"
#       AND c.mag_i BETWEEN 1x.0 AND 18.x
#     '''
#     name = 'bhm_csc_boss_bright'
#     this_cadence = 'bright_1x1'
#     instrument = 'BOSS'
#
#     def build_query(self, version_id, query_region=None):
#         query = super().build_query(version_id, query_region)
#         t = self.alias_t
#         query = query.where(
#             (t.mag_i >= self.parameters['mag_i_min']) &
#             (t.mag_i < self.parameters['mag_i_max'])
#         )
#         return query
#
#
# class BhmCscApogeeCarton(BhmCscBaseCarton):
#     '''
#     SELECT * from bhm_csc AS c
#     WHERE c.spectrograph = "APOGEE"
#       AND c.mag_H BETWEEN 1x.x AND 1x.x
#     '''
#     name = 'bhm_csc_apogee'
#     this_cadence = 'bright_3x1'
#     instrument = 'APOGEE'
#
#     def build_query(self, version_id, query_region=None):
#         query = super().build_query(version_id, query_region)
#         t = self.alias_t
#         query = query.where(
#             (t.mag_h >= self.parameters['mag_h_min']) &
#             (t.mag_h < self.parameters['mag_h_max'])
#         )
#         return query
#
# #
# # '''
# # Exporting from the temp table
# #
# # \copy (SELECT * FROM sandbox.temp_bhm_csc_boss_dark)  TO 'bhm_csc_boss_dark.csv' with csv header
# # \copy (SELECT * FROM sandbox.temp_bhm_csc_boss_bright)
# #      TO 'bhm_csc_boss_bright.csv' with csv header
# # \copy (SELECT * FROM sandbox.temp_bhm_csc_apogee)  TO 'bhm_csc_apogee.csv' with csv header
# # \copy (SELECT * FROM sandbox.temp_td_csc_subset)  TO 'temp_td_csc_subset.csv' with csv header
# # for F in bhm_csc_*.csv; do
# #    stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic;
# # done
# #
# # '''
#
