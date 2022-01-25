#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import JOIN, fn

from sdssdb.peewee.sdss5db.catalogdb import (SDSSV_BOSS_SPALL,
                                             BHM_CSC_v2, BHM_eFEDS_Veto,
                                             CatalogToBHM_eFEDS_Veto,
                                             CatalogToPanstarrs1,
                                             CatalogToSDSS_DR16_SpecObj,
                                             CatalogToTIC_v8, Panstarrs1,
                                             SDSS_DR16_SpecObj,
                                             SDSSV_BOSS_Conflist,
                                             SDSSV_Plateholes,
                                             SDSSV_Plateholes_Meta, TIC_v8)

from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import AB2Jy


# Details: Start here
# https://wiki.sdss.org/display/OPS/Cartons+for+v0.5#Cartonsforv0.5-BHMCSC # noqa: E501

# This module provides the following BHM cartons:
#   bhm_csc_boss
#   bhm_csc_apogee


class BhmCscBossCarton(BaseCarton):
    '''
    SELECT *
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
                AND p.i_stk_psf_flux BETWEEN xxx AND xxxx
                AND p.i_stk_psf_flux != 'NaN')
           OR (     x.ocat = 'G'
               AND (NOT tic.gaiamag IS NULL)
               AND tic.gaiamag BETWEEN 13. AND xxxx)
           )
           AND ( COALESCE(c2p.version_id,c2t.version_id) = 25 )
           AND ( COALESCE(c2p.best,c2t.best) = True )
    ;
    '''

    name = 'bhm_csc_boss'
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

        gaia_g_mag_min = self.parameters['gaia_g_mag_min']
        gaia_rp_mag_min = self.parameters['gaia_rp_mag_min']

        i_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence1'])
        gaia_g_mag_max_for_cadence1 = self.parameters['gaia_g_mag_max_for_cadence1']
        i_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence2'])
        gaia_g_mag_max_for_cadence2 = self.parameters['gaia_g_mag_max_for_cadence2']

        value = peewee.Value(self.parameters['value'])

        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence = peewee.Case(
            None,
            (
                (((ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1) |
                  (tic.gaiamag < gaia_g_mag_max_for_cadence1)),
                 cadence1),
                (((ps.i_stk_psf_flux > i_psf_flux_min_for_cadence2) |
                  (tic.gaiamag < gaia_g_mag_max_for_cadence2)),
                 cadence2),
            ),
            cadence3)

        # #########################################################################
        # prepare the spectroscopy catalogues
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']
        dpriority_has_spec = self.parameters['dpriority_has_spec']

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

        # adjust priority if target aleady has an SDSS spectrum
        priority_1 = peewee.Case(
            None,
            (
                (s16.c.specobjid.is_null(False), 1),  # any of these can be satisfied
                (s2020.c.pk.is_null(False), 1),
                (sV.c.specobjid.is_null(False), 1),
                (sph.c.pkey.is_null(False), 1),
            ),
            0)
        #
        # Compute net priority
        priority_floor_dark = peewee.Value(self.parameters['priority_floor_dark'])
        priority_floor_bright = peewee.Value(self.parameters['priority_floor_bright'])
        priority = peewee.Case(
            None,
            (
                (((ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1) |
                  (tic.gaiamag < gaia_g_mag_max_for_cadence1)),
                 priority_floor_bright + (priority_1 * dpriority_has_spec) + x.pri - 1),
            ),
            priority_floor_dark + (priority_1 * dpriority_has_spec) + x.pri - 1)

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
                ps.i_stk_psf_flux.alias('ps_i_stk_psf_flux'),
                tic.gaiamag.alias("tic_gaiamag"),
                priority.alias('priority'),
                cadence.alias('cadence'),
                value.alias('value'),
            )
            .join(ps, join_type=JOIN.LEFT_OUTER,
                  on=(x.idps == ps.extid_hi_lo))
            .join(tic, join_type=JOIN.LEFT_OUTER,
                  on=(x.idg2 == tic.gaia_int))
            .join(c2ps, join_type=JOIN.LEFT_OUTER,
                  on=(ps.catid_objid == c2ps.target_id))
            .join(c2tic, join_type=JOIN.LEFT_OUTER,
                  on=(tic.id == c2tic.target_id))
            .where(
                fn.coalesce(c2ps.version_id, c2tic.version_id) == version_id,
                fn.coalesce(c2ps.best, c2tic.best) >> True,
            )
            .where(
                (
                    (x.ocat == 'P') &
                    (ps.g_stk_psf_flux < g_psf_flux_max) &
                    (ps.r_stk_psf_flux < r_psf_flux_max) &
                    (ps.i_stk_psf_flux < i_psf_flux_max) &
                    (ps.z_stk_psf_flux < z_psf_flux_max) &
                    (ps.i_stk_psf_flux != 'NaN')
                ) |
                (
                    (x.ocat == 'G') &
                    (tic.gaiamag > gaia_g_mag_min) &
                    (tic.gaiarp > gaia_rp_mag_min)
                )
            )
            .distinct(fn.coalesce(c2ps.catalogid, c2tic.catalogid))
        )

        # Append the spectro query
        query = (
            query
            .join(c2s16, JOIN.LEFT_OUTER,
                  on=(c2s16.catalogid == fn.coalesce(c2ps.catalogid, c2tic.catalogid)))
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .join(c2s2020, JOIN.LEFT_OUTER,
                  on=(c2s2020.catalogid == fn.coalesce(c2ps.catalogid, c2tic.catalogid)))
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
                                x.ora, x.odec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                x.ora, x.odec,
                                match_radius_spectro)
                )
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(x.ora, x.odec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query


# #######
class BhmCscApogeeCarton(BaseCarton):
    '''
    SELECT * from bhm_csc_v2 AS x
    WHERE
      AND x.hmag BETWEEN 7 AND 14
    '''
    name = 'bhm_csc_apogee'
    mapper = 'BHM'
    program = 'bhm_csc'
    category = 'science'
    this_cadence = 'bright_3x1'
    instrument = 'APOGEE'
    tile = False

    def build_query(self, version_id, query_region=None):
        x = BHM_CSC_v2.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        hmag_max_for_cadence1 = self.parameters['hmag_max_for_cadence1']

        value = peewee.Value(self.parameters['value'])

        cadence = peewee.Case(
            None,
            (
                ((x.hmag < hmag_max_for_cadence1), self.parameters['cadence1']),
            ),
            self.parameters['cadence2'])

        # Compute net priority
        priority = peewee.Value(self.parameters['priority_floor']) + x.pri - 1

        query = (
            x.select(
                c2tic.catalogid.alias('catalogid'),
                x.cxoid.alias('csc_cxoid'),
                x.ra2m.alias('csc_ra2m'),
                x.dec2m.alias('csc_dec2m'),
                x.hmag.alias('csc_hmag'),
                x.designation2m.alias('csc_designation2m'),
                priority.alias('priority'),
                cadence.alias('cadence'),
                value.alias('value'),
            )
            .join(tic, on=(x.designation2m == tic.twomass_psc))
            .join(c2tic, on=(tic.id == c2tic.target_id))
            .where(
                c2tic.version_id == version_id,
                c2tic.best >> True,
                x.hmag >= self.parameters['hmag_min'],
                x.hmag < self.parameters['hmag_max'],
                x.hmag != 'NaN',
            )
            .distinct(c2tic.catalogid)
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(x.ra2m, x.dec2m,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
