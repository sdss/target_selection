#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import JOIN, fn

from sdssdb.peewee.sdss5db.catalogdb import (
    BHM_CSC_v2,
    CatalogToPanstarrs1,
    CatalogToSDSS_DR19p_Speclite,
    CatalogToGaia_DR2,
    Panstarrs1,
    SDSS_DR19p_Speclite,
    Gaia_DR2,
    TwoMassPSC,
    CatalogToTwoMassPSC,
#    TIC_v8,
#    CatalogToTIC_v8,
)

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
        g2 = Gaia_DR2.alias()
        c2g2 = CatalogToGaia_DR2.alias()

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

        value = peewee.Value(self.parameters['value']).cast('real')

        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence = peewee.Case(
            None,
            (
                (((ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1) |
                  (g2.phot_g_mean_mag < gaia_g_mag_max_for_cadence1)),
                 cadence1),
                (((ps.i_stk_psf_flux > i_psf_flux_min_for_cadence2) |
                  (g2.phot_g_mean_mag < gaia_g_mag_max_for_cadence2)),
                 cadence2),
            ),
            cadence3)

        # #########################################################################
        # prepare the existing spectroscopy catalogues - these affect only priority
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']
        dpriority_has_spec = self.parameters['dpriority_has_spec']

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

        # adjust priority if target aleady has a good SDSS spectrum
        priority_1 = peewee.Case(
            None,
            (
                (s19.c.s19_pk.is_null(False), 1),
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
                  (g2.phot_g_mean_mag < gaia_g_mag_max_for_cadence1)),
                 priority_floor_bright + (priority_1 * dpriority_has_spec) + x.pri - 1),
            ),
            priority_floor_dark + (priority_1 * dpriority_has_spec) + x.pri - 1)

        query = (
            x.select(
                fn.coalesce(fn.max(c2ps.catalogid), fn.max(c2g2.catalogid)).alias('catalogid'),
                fn.max(x.cxoid).alias('csc_cxoid'),
                fn.max(x.ora).alias('csc_ora'),
                fn.max(x.odec).alias('csc_odec'),
                fn.max(x.omag).alias('csc_omag'),
                fn.max(x.ocat).alias('csc_ocat'),
                fn.max(x.idps).alias('csc_idps'),
                fn.max(x.idg2).alias('csc_idg2'),
                fn.max(ps.i_stk_psf_flux).alias('ps_i_stk_psf_flux'),
                fn.max(g2.phot_g_mean_mag).alias("g2_g_mag"),
                fn.max(priority).alias('priority'),
                fn.max(cadence).alias('cadence'),
                fn.max(value).alias('value'),
            )
            .join(ps, join_type=JOIN.LEFT_OUTER,
                  on=(x.idps == ps.extid_hi_lo))
            .join(g2, join_type=JOIN.LEFT_OUTER,
                  on=(x.idg2 == g2.source_id))
            .join(c2ps, join_type=JOIN.LEFT_OUTER,
                  on=(ps.catid_objid == c2ps.target_id))
            .join(c2g2, join_type=JOIN.LEFT_OUTER,
                  on=(g2.source_id == c2g2.target_id))
            .join(c2s19, join_type=JOIN.LEFT_OUTER,
                  on=(fn.coalesce(fn.max(c2ps.catalogid),
                                  fn.max(c2g2.catalogid)) == c2s19.catalogid))
            .join(s19, join_type=JOIN.LEFT_OUTER,
                  on=(s19.c.s19_pk == c2s19.target_id))
            .where(
                fn.coalesce(c2ps.version_id, c2g2.version_id) == version_id,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
                fn.coalesce(c2ps.best, c2g2.best) >> True,
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
                    (g2.phot_g_mean_mag > gaia_g_mag_min) &
                    (g2.phot_rp_mean_mag > gaia_rp_mag_min)
                )
            )
            .distinct(fn.coalesce(fn.max(c2ps.catalogid), fn.max(c2g2.catalogid)))
            .group_by(x.cxoid)
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
    this_cadence = 'bright_3x1'  # TODO Check this is still a valid choice
    instrument = 'APOGEE'
    tile = False

    def build_query(self, version_id, query_region=None):
        x = BHM_CSC_v2.alias()
        # tic = TIC_v8.alias()
        # c2tic = CatalogToTIC_v8.alias()
        # change to now rely directly on twomass_psc table instead of TIC_v8
        tm = TwoMassPSC.alias()
        c2tm = CatalogToTwoMassPSC.alias()

        hmag_max_for_cadence1 = self.parameters['hmag_max_for_cadence1']

        value = peewee.Value(self.parameters['value']).cast('real')

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
                # c2tic.catalogid.alias('catalogid'),
                c2tm.catalogid.alias('catalogid'),
                x.cxoid.alias('csc_cxoid'),
                x.ra2m.alias('csc_ra2m'),
                x.dec2m.alias('csc_dec2m'),
                x.hmag.alias('csc_hmag'),
                x.designation2m.alias('csc_designation2m'),
                priority.alias('priority'),
                cadence.alias('cadence'),
                value.alias('value'),
            )
            # .join(tic, on=(x.designation2m == tic.twomass_psc))
            # .join(c2tic, on=(tic.id == c2tic.target_id))
            .join(tm, on=(x.designation2m == tm.designation))
            .join(c2tm, on=(tm.designation == c2tm.target_id))
            .where(
                c2tm.version_id == version_id,
                c2tm.best >> True,
                # c2tic.version_id == version_id,
                #  c2tic.best >> True,
                x.hmag >= self.parameters['hmag_min'],
                x.hmag < self.parameters['hmag_max'],
                x.hmag != 'NaN',
            )
            .distinct(x.cxoid)
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(x.ra2m, x.dec2m,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
