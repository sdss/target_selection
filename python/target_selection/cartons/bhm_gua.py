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
    CatalogToTIC_v8,
    TIC_v8,
    # Gaia_DR2,
    # CatalogToGaia_unWISE_AGN,  # <-- this is old and does not work
    CatalogToSDSS_DR16_SpecObj,
    SDSS_DR16_SpecObj,
    Gaia_unWISE_AGN,
    CatalogToBHM_eFEDS_Veto,
    BHM_eFEDS_Veto,
    SDSSV_BOSS_SPALL,
    SDSSV_BOSS_Conflist,
    SDSSV_Plateholes,
    SDSSV_Plateholes_Meta,
)
from target_selection.cartons.base import BaseCarton

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms
#
# This module provides the following BHM cartons:
# bhm_gua_dark
# bhm_gua_bright
#
# Updated by TD on 27/01/2021 to add v0.5 improvements
# switch to a q3c join to sdss_dr16_specobj

'''
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
'''


class BhmGuaBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for both Gaia UnWISE AGN cartons
    To be sub-classed, not to be called directly.

    To get from Catalog to GUA we join tables via :
    Catalog -> CatalogToTIC_v8 -> Gaia_DR2 -> Gaia_unWISE_AGN

    Cadence+priority per target depends on brightness
    Keep the two bhm_gua_* cartons separate to enable overlap in brightness ranges
    '''

    name = 'bhm_gua_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_filler'
    tile = False
    priority = None
    cadence = None
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        # ## c2t = CatalogToGaia_unWISE_AGN.alias() - deprecated - but leave this as a reminder
        c2tic = CatalogToTIC_v8.alias()
        tic = TIC_v8.alias()
        # s2020 = BHM_eFEDS_Veto.alias()
        # sV = SDSSV_BOSS_SPALL.alias()
        # ph = SDSSV_Plateholes.alias()
        # phm = SDSSV_Plateholes_Meta.alias()

        # g2 = Gaia_DR2.alias()
        t = Gaia_unWISE_AGN.alias()

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
                ssV.specobjid.is_null()
            )
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
                ssph.pkey.is_null()
            )
        )

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get('priority', 10000)))
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        inertial = peewee.Value(True)
        cadence = peewee.Value(self.parameters['cadence'])
        instrument = peewee.Value(self.instrument)

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

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

        bp_rp = t.bp - t.rp
        g = (t.g + coeffs['g0_p'] + coeffs['g1_p'] * bp_rp + coeffs['g2_p'] * bp_rp * bp_rp +
             coeffs['g3_p'] * bp_rp * bp_rp * bp_rp)
        r = (t.g + coeffs['r0_p'] + coeffs['r1_p'] * bp_rp + coeffs['r2_p'] * bp_rp * bp_rp +
             coeffs['r3_p'] * bp_rp * bp_rp * bp_rp)
        i = (t.g + coeffs['i0_p'] + coeffs['i1_p'] * bp_rp + coeffs['i2_p'] * bp_rp * bp_rp +
             coeffs['i3_p'] * bp_rp * bp_rp * bp_rp)
        z = (t.g + coeffs['z0_p'] + coeffs['z1_p'] * bp_rp + coeffs['z2_p'] * bp_rp * bp_rp +
             coeffs['z3_p'] * bp_rp * bp_rp * bp_rp)

        # validity checks - set limits semi-manually
        bp_rp_min = 0.0
        bp_rp_max = 1.8
        valid = (t.g.between(0.1, 29.9) &
                 t.bp.between(0.1, 29.9) &
                 t.rp.between(0.1, 29.9) &
                 bp_rp.between(bp_rp_min, bp_rp_max))

        opt_prov = peewee.Case(None, ((valid, 'sdss_psfmag_from_gaiadr2'),), 'undefined')
        magnitude_g = peewee.Case(None, ((valid, g),), 'NaN')
        magnitude_r = peewee.Case(None, ((valid, r),), 'NaN')
        magnitude_i = peewee.Case(None, ((valid, i),), 'NaN')
        magnitude_z = peewee.Case(None, ((valid, z),), 'NaN')

        # Create temporary tables for the base query and the Q3C cross-match
        # tables.

        bquery = (
            c.select(
                c.catalogid,
                c.ra,   # extra
                c.dec,   # extra
                t.gaia_sourceid,   # extra
                t.unwise_objid,   # extra
                priority.alias('priority'),
                value.alias('value'),
                inertial.alias('inertial'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                opt_prov.alias('optical_prov'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                t.g.alias('gaia_g'),
                t.bp.alias('bp'),
                t.rp.alias('rp'),
                t.w1.alias('gua_w1'),   # extra
                t.w2.alias('gua_w2'),   # extra
                t.prob_rf.alias('gua_prob_rf'),   # extra
                t.phot_z.alias('gua_phot_z'),   # extra
                # rely on the centralised magnitude routines for 'real' griz, bp,rp,gaia_g
            )
            .join(c2tic)
            .join(tic)
            # .join(g2)    # can skip this join using the gaia_int from the TIC
            # .join(t, on=(g2.source_id == t.gaia_sourceid))
            .join(t, on=(tic.gaia_int == t.gaia_sourceid))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid)
                    # (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk)
                    # (c2s2020.version_id == version_id)
                )
            )
            # finished joining the spectroscopy
            .where(
                c.version_id == version_id,
                # c2tic.version_id == version_id,
                c2tic.best >> True,
            )
            .where(
                (t.prob_rf >= self.parameters['prob_rf_min']),
                (t.g >= self.parameters['mag_g_min']),
                (t.rp >= self.parameters['mag_rp_min']),
                (
                    (t.g < self.parameters['mag_g_max']) |
                    (t.rp < self.parameters['mag_rp_max'])
                ),
            )
            # then reject any GUA targets with existing good DR16+SDSS-V spectroscopy
            .where(
                s16.c.specobjid.is_null(True),
                s2020.c.pk.is_null(True)
            )
            # avoid duplicates - trust the gaia ids in the GUA parent sample
            .distinct([t.gaia_sourceid])
        )

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            bquery = (bquery
                      .where(peewee.fn.q3c_radial_query(c.ra,
                                                        c.dec,
                                                        query_region[0],
                                                        query_region[1],
                                                        query_region[2])))

        self.log.debug('Creating temporary table for base query ...')
        bquery.create_table(self.name + '_bquery', temporary=True)
        self.database.execute_sql(f'CREATE INDEX ON {self.name}_bquery (ra, dec)')
        self.database.execute_sql(f'ANALYZE {self.name}_bquery')

        sph.create_table(self.name + '_sph', temporary=True)
        self.database.execute_sql(f'CREATE INDEX ON {self.name}_sph (target_ra, target_dec)')
        self.database.execute_sql(f'ANALYZE {self.name}_sph')

        sV.create_table(self.name + '_sv', temporary=True)
        self.database.execute_sql(f'CREATE INDEX ON {self.name}_sv (plug_ra, plug_dec)')
        self.database.execute_sql(f'ANALYZE {self.name}_sv')

        bquery_table = peewee.Table(f'{self.name}_bquery', alias='bquery')
        sph_table = peewee.Table(f'{self.name}_sph')
        sV_table = peewee.Table(f'{self.name}_sv')

        query = (
            bquery_table
            .select(peewee.SQL('bquery.*'))
            .join(
                sV_table, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(bquery_table.c.ra, bquery_table.c.dec,
                                sV_table.c.plug_ra, sV_table.c.plug_dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph_table, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(bquery_table.c.ra, bquery_table.c.dec,
                                sph_table.c.target_ra, sph_table.c.target_dec,
                                match_radius_spectro)
                )
            )
            # then reject any GUA targets with existing good SDSS-V spectroscopy or a platehole
            .where(
                sV_table.c.specobjid.is_null(True),
                sph_table.c.pkey.is_null(True),
            )
        )

        return query


class BhmGuaDarkCarton(BhmGuaBaseCarton):
    '''
    -------  bhm_gua_dark   ------
    SQL as above plus:

        AND ( gua.g > 16.x AND gua.rp > 16.x)
    '''
    name = 'bhm_gua_dark'


class BhmGuaBrightCarton(BhmGuaBaseCarton):
    '''
    ------  bhm_gua_bright   ------
    SQL as above plus:

        AND ( gua.g < 18.x OR gua.rp < 18.x)

    '''
    name = 'bhm_gua_bright'
