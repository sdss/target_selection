#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-31
# @Filename: bhm_rm.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn

from target_selection.cartons.base import BaseCarton
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    BHM_RM_v0_2,
    CatalogToBHM_RM_v0,
    BHM_RM_Tweaks,
    SDSSV_BOSS_SPALL,
    SDSSV_Plateholes,
    SDSSV_Plateholes_Meta,
)


#  This module provides the following BHM cartons in v05:
#  bhm_rm_core
#  bhm_rm_known_spec
#  bhm_rm_var
#  bhm_rm_ancillary

pmsig_min = -3.0
plxsig_min = -3.0


class BhmRmBaseCarton(BaseCarton):
    '''
    This class provides common setting and the masking routines used by all RM cartons
    '''

    name = 'bhm_rm_base'
    base_name = 'bhm_rm_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_rm'
    instrument = 'BOSS'
    tile = False
    priority = None
    inertial = True
    alias_c = None
    alias_t = None
    alias_tw = None

    def get_fieldlist(self):
        '''Read the RM field centres from the yaml'''
        fieldlist = []
        base_parameters = self.config['parameters'].get(self.base_name, None)
        if base_parameters:
            fieldlist = base_parameters['fieldlist']
        return fieldlist

    def append_spatial_query(self, query, fieldlist):
        '''extend the peewee query using a list of field centres'''
        if fieldlist is None:
            return query
        elif len(fieldlist) == 0:
            return query

        q = False
        for f in fieldlist:
            q = (q | peewee.fn.q3c_radial_query(self.alias_c.ra,
                                                self.alias_c.dec,
                                                f['racen'],
                                                f['deccen'],
                                                f['radius']))
        return query.where(q)

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2t = CatalogToBHM_RM_v0.alias()
        t = BHM_RM_v0_2.alias()
        stw = BHM_RM_Tweaks.alias()
        self.alias_c = c
        self.alias_t = t

        fieldlist = self.get_fieldlist()

        tw = (
            stw.
            select(
                stw.pkey.alias('pkey'),
                stw.ra.alias('ra'),
                stw.dec.alias('dec'),
                stw.rm_suitability.alias('rm_suitability'),
            ).
            where(
                (stw.date_set == '30-Nov-2020') |
                (stw.date_set == '25-May-2021')
            )
        )
        self.alias_tw = tw
        # #########################################################################
        # prepare the spectroscopy catalogues

        # SDSS-V spAll - select only objects we want to exclude on
        # the basis of their pipeline classifications
        # Currently this is only for secure STARs in the COSMOS field
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
                fn.rank().over(partition_by=[ssV.catalogid],
                               order_by=[ssV.sn_median_all.desc()]).alias('sn_rank'),
            )
            .where(
                ssV.programname.contains('RM'),
                ssV.firstcarton.contains('bhm_rm_'),
                ssV.class_ == 'STAR',
                ssV.zwarning == 0,
                ssV.sn_median_all > 2.0,
                # select only COSMOS plates
                ssV.plate << [15038, 15070, 15071, 15252, 15253, 15289]
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped and that have firstcarton ~ 'bhm_rm_'
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
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
            .where(
                ssph.holetype == 'BOSS_SHARED',
                ssph.sourcetype == 'SCI',
                ssph.firstcarton.contains('bhm_rm_'),
                ssphm.isvalid > 0,
            )
            .distinct([ssph.catalogid])
            .alias('sph')
        )

        # fold in tiers of magnitude-based priority
        priority_mag_step = 0.5
        priority_mag_bright = 17.0
        priority_mag_faint = 22.0
        priority_mag_bright_known_spec = 20.5
        priority_floor = self.parameters.get('priority', 10000)
        priority1 = peewee.Case(
            None,
            (
                (
                    (t.mi <= priority_mag_bright),
                    priority_floor + 0
                ),
                (
                    (
                        (self.name == 'bhm_rm_known_spec') &
                        ~(t.field_name.contains('SDSS-RM')) &
                        (t.mi <= priority_mag_bright_known_spec)
                    ),
                    priority_floor + 0
                ),
                (
                    (t.mi <= priority_mag_faint),
                    priority_floor +
                    5 * (1 + peewee.fn.floor((t.mi - priority_mag_bright) /
                                             priority_mag_step).cast('int'))
                ),
                (
                    (t.mi > priority_mag_faint),
                    priority_floor + 95
                ),
            ),
            None
        )
        # # this secondary priority rule is based on whether this target was
        # # assigned a platehole during the SDSSV plate programme
        # # boost the priorities of those targets that were put onto plates
        # priority2 = peewee.Case(
        #     None,
        #     (
        #         (sph.c.pkey.is_null(False), -100),
        #         (sph.c.pkey.is_null(True), 0),
        #     ),
        #     None
        # )

        # this secondary priority rule boosts the priority of targets that
        # have rm_suitability >= 1 in the bhm_rm_tweaks table
        priority2 = peewee.Case(None, ((tw.c.rm_suitability >= 1, -100), ), 0)

        # combine the two priorities
        priority = priority1 + priority2

        # this just checks if this target was
        # assigned a platehole during the SDSSV plate programme
        # for information only - no action taken
        in_SDSSV_plates = peewee.Case(
            None,
            (
                (sph.c.pkey.is_null(False), True),
            ),
            False
        ).cast('bool')

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast('bool')
        match_radius_spectro = 1.0 / 3600.0

        # This is the scheme used in v0
        cadence_v0 = peewee.Case(None,
                                 (
                                     (t.field_name.contains('S-CVZ'), 'bhm_rm_lite5_100x8'),
                                 ),
                                 'bhm_rm_174x8')

        # this gives the new names for the same cadences assumed in v0
        cadence_v0p5 = peewee.Case(None,
                                   (
                                       (t.field_name.contains('S-CVZ'), 'dark_100x8'),
                                   ),
                                   'dark_174x8')

        # the following will replace old generic cadences when relevant table has been populated
        # TODO - replace when correct cadences are loaded
        cadence_v1p0 = peewee.Case(None,
                                   (
                                       (t.field_name.contains('SDSS-RM'), 'bhm_rm_sdss-rm'),
                                       (t.field_name.contains('COSMOS'), 'bhm_rm_cosmos'),
                                       (t.field_name.contains('XMM-LSS'), 'bhm_rm_xmm-lss'),
                                       (t.field_name.contains('S-CVZ'), 'bhm_rm_cvz-s'),
                                       (t.field_name.contains('CDFS'), 'bhm_rm_cdfs'),
                                       (t.field_name.contains('ELIAS-S1'), 'bhm_rm_elias-s1'),
                                   ),
                                   'dark_174x8')

        # Photometric precedence: DES>PS1>SDSS(>Gaia)>NSC.
        opt_prov = peewee.Case(None,
                               (
                                   (t.sdss == 1, 'sdss_psfmag'),
                                   (t.des == 1, 'psfmag'),
                                   (t.ps1 == 1, 'ps_psfmag'),
                                   (t.optical_survey == 'Gaia', 'other'),
                                   (t.nsc == 1, 'psfmag'),
                               ),
                               'other')

        magnitude_g = peewee.Case(None,
                                  (
                                      ((t.sdss == 1) & (t.psfmag_sdss[1] > 0.0), t.psfmag_sdss[1]),
                                      ((t.des == 1) & (t.psfmag_des[0] > 0.0), t.psfmag_des[0]),
                                      ((t.ps1 == 1) & (t.psfmag_ps1[0] > 0.0), t.psfmag_ps1[0]),
                                      ((t.optical_survey == 'Gaia') & (t.mag_gaia[0] > 0.0),
                                       t.mag_gaia[0]),  # just using gaia G for now
                                      ((t.nsc == 1) & (t.mag_nsc[0] > 0.0), t.mag_nsc[0]),
                                  ),
                                  99.9)  # should never get here
        magnitude_r = peewee.Case(None,
                                  (
                                      ((t.sdss == 1) & (t.psfmag_sdss[2] > 0.0), t.psfmag_sdss[2]),
                                      ((t.des == 1) & (t.psfmag_des[1] > 0.0), t.psfmag_des[1]),
                                      ((t.ps1 == 1) & (t.psfmag_ps1[1] > 0.0), t.psfmag_ps1[1]),
                                      ((t.nsc == 1) & (t.mag_nsc[1] > 0.0), t.mag_nsc[1]),
                                  ),
                                  99.9)  # should never get here
        magnitude_i = peewee.Case(None,
                                  (
                                      ((t.sdss == 1) & (t.psfmag_sdss[3] > 0.0), t.psfmag_sdss[3]),
                                      ((t.des == 1) & (t.psfmag_des[2] > 0.0), t.psfmag_des[2]),
                                      ((t.ps1 == 1) & (t.psfmag_ps1[2] > 0.0), t.psfmag_ps1[2]),
                                      ((t.nsc == 1) & (t.mag_nsc[2] > 0.0), t.mag_nsc[2]),
                                      (t.mi > 0.0, t.mi),
                                      ((t.optical_survey == 'Gaia') & (t.mag_gaia[2] > 0.0),
                                       t.mag_gaia[2]),  # just using gaia RP for now
                                  ),
                                  99.9)  # should never get here
        magnitude_z = peewee.Case(None,
                                  (
                                      ((t.sdss == 1) & (t.psfmag_sdss[4] > 0.0), t.psfmag_sdss[4]),
                                      ((t.des == 1) & (t.psfmag_des[3] > 0.0), t.psfmag_des[3]),
                                      ((t.ps1 == 1) & (t.psfmag_ps1[3] > 0.0), t.psfmag_ps1[3]),
                                      ((t.nsc == 1) & (t.mag_nsc[3] > 0.0), t.mag_nsc[3]),
                                  ),
                                  99.9)  # should never get here

        query = (
            c.select(
                c.catalogid,
                c.ra,  # extra
                c.dec,  # extra
                t.field_name.alias('rm_field_name'),  # extra
                t.pk.alias('rm_pk'),  # extra
                instrument.alias('instrument'),
                priority.alias('priority'),
                priority1.alias('priority1'),
                priority2.alias('priority2'),
                value.alias('value'),
                cadence_v0p5.alias('cadence'),
                cadence_v0.alias('cadence_v0'),  # extra
                cadence_v0p5.alias('cadence_v0p5'),  # extra
                cadence_v1p0.alias('cadence_v1p0'),  # extra
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                opt_prov.alias('optical_prov'),
                inertial.alias('inertial'),
                t.optical_survey.alias('optical_survey'),  # extra
                c2t.best.alias("c2t_best"),  # extra
                in_SDSSV_plates.alias('in_SDSSV_plates'),  # extra
                tw.c.rm_suitability.cast('int').alias('rm_suitability'),  # extra
            )
            .join(c2t)
            # An explicit join is needed because we are using c2t for Catalog_to_BHM_RM_v0
            # rather than a native c2t for Catalog_to_BHM_RM_v0_2
            .join(t, on=(c2t.target_id == t.pk))
            .where(
                c.version_id == version_id,
                c2t.version_id == version_id,
                # c2t.best >> True   # TODO check if this is dropping RM targets
                #                    # like it does for AQMES
            )
            .where
            (
                (
                    (t.mi >= self.parameters['mag_i_min']) &
                    (t.mi < self.parameters['mag_i_max'])
                ) |
                (
                    # S-CVZ targets often have only Gaia photom
                    (t.field_name.contains('S-CVZ')) &
                    (t.mg >= self.parameters['mag_g_min_cvz_s']) &
                    (t.mg < self.parameters['mag_g_max_cvz_s'])
                )
            )
            .switch(c)
            .join(
                tw,
                JOIN.LEFT_OUTER,
                on=(fn.q3c_join(tw.c.ra, tw.c.dec,
                                c.ra, c.dec,
                                match_radius_spectro))
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (sV.c.sn_rank == 1)   # only consider the best spectrum per object
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
            .where(
                # Reject any objects where the highest SNR spectrum for
                # this target in sdssv_boss_spall is classified as STAR
                sV.c.specobjid.is_null(True),
                #
                # Reject any targets that are flagged as being unsuitable for RM in bhm_rm_tweaks
                # bhm_rm_tweaks.rm_suitability==0 means:
                # 'target is probably unsuitable for RM, do not observe in the future'
                (
                    tw.c.pkey.is_null(True) |
                    (tw.c.rm_suitability != 0)
                )
            )
            .distinct([t.pk])   # avoid duplicates - trust the RM parent sample
                                # - only needed if NOT using c2t.best = True condition
        )
        query = self.append_spatial_query(query, fieldlist)

        return query


class BhmRmCoreCarton(BhmRmBaseCarton):
    '''
    bhm_rm_core: select all photometric QSO targets with the
                 likelihood method (Skewt), flux-limited to 21.5 in i-band PSF mag

    SELECT * FROM bhm_rm
            WHERE skewt_qso = 1
            AND WHERE mi BETWEEN 15.0 AND 21.5
            AND WHERE pmsig < 5.0
            AND WHERE plxsig < 5.0


    also require t.skewt_qso_prior == 1 in CVZ-S
    '''

    name = 'bhm_rm_core'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.skewt_qso == 1),
            (
                ~(t.field_name.contains('S-CVZ')) |
                ((t.field_name.contains('S-CVZ')) & (t.skewt_qso_prior == 1))
            ),
            (t.pmsig < self.parameters['pmsig_max']),
            (t.plxsig < self.parameters['plxsig_max']),
            (
                # catch bad photometry - require at least gri detections in at least one system
                # in all fields - except in S-CVZ where targets often have only Gaia photom
                (
                    (t.sdss == 1) &
                    (t.psfmag_sdss[1] > 0.0) &
                    (t.psfmag_sdss[2] > 0.0) &
                    (t.psfmag_sdss[3] > 0.0)
                ) |
                (
                    (t.ps1 == 1) &
                    (t.psfmag_ps1[0] > 0.0) &
                    (t.psfmag_ps1[1] > 0.0) &
                    (t.psfmag_ps1[2] > 0.0)
                ) |
                (
                    (t.des == 1) &
                    (t.psfmag_des[0] > 0.0) &
                    (t.psfmag_des[1] > 0.0) &
                    (t.psfmag_des[2] > 0.0)
                ) |
                (
                    (t.nsc == 1) &
                    (t.mag_nsc[0] > 0.0) &
                    (t.mag_nsc[1] > 0.0) &
                    (t.mag_nsc[2] > 0.0)
                ) |
                (
                    (t.field_name.contains('S-CVZ')) &
                    (t.mg > 0.0)
                )
            ),
            #  (t.pmsig > pmsig_min) &  # this catches cases with NULL=-9
            #  (t.plxsig > plxsig_min) &
            ~(t.field_name.contains('SDSS-RM')),  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmKnownSpecCarton(BhmRmBaseCarton):
    '''
    bhm_rm_known_spec:  select all spectroscopically confirmed QSOs where redshift is extragalactic

    SELECT * FROM bhm_rm
        WHERE specz > 0.005
        AND WHERE mi BETWEEN 15.0 AND 21.5   # <- exact limits TBD

    For SDSS-RM field select only QSOs that were part of the original SDSS-III/IV RM program
    Do that via a bit in the spec_bitmask field. Also restrict to mi < 21.0
    '''

    name = 'bhm_rm_known_spec'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        tw = self.alias_tw
        spec_bitmask_sdss_rm_qso = 2**3
        query = query.where(
            (
                (t.spec_q == 1) |
                (
                    (tw.c.pkey.is_null(False)) &
                    (tw.c.rm_suitability == 1)
                )

            ),
            (
                ~(t.field_name.contains('SDSS-RM')) |
                (
                    # include extra constraints on SDSS-RM targets
                    (t.spec_bitmask.bin_and(spec_bitmask_sdss_rm_qso) != 0) &
                    (t.mi < self.parameters['mag_i_max_sdss_rm'])
                )
            ),
            (
                ~(t.field_name.contains('COSMOS')) |
                (
                    # include extra constraints on COSMOS targets
                    (t.mi < self.parameters['mag_i_max_cosmos'])
                )
            ),
            (
                ~(t.field_name.contains('XMM-LSS')) |
                (
                    # include extra constraints on XMM-LSS targets
                    (t.mi < self.parameters['mag_i_max_xmm_lss'])
                )
            ),
            (
                (t.specz.between(self.parameters['specz_min'],
                                 self.parameters['specz_max'])) |
                (tw.c.rm_suitability == 1)
                # allow this here because recently observed QSOs will
                # not yet have specz in BHM_RM_v0_2
            ),
        )

        return query


class BhmRmVarCarton(BhmRmBaseCarton):
    '''bhm_rm_var: selected based on g-band variability > 0.05 mag
                   and bright enough to be detected by Gaia (G<~21)

    SELECT * FROM bhm_rm
        WHERE ( (des_var_sn[0] > 5.0 AND des_var_rms[0] > 0.05)  OR
                (ps1_var_sn[0]>5.0 AND ps1_var_rms[0]>0.05))
        AND WHERE mi BETWEEN 15.0 AND 21.5    # <- exact limits TBD
        AND WHERE pmsig < x.x
        AND WHERE plxsig < x.x
        AND WHERE gaia = 1

    #debug select t.pk, t.ra, t.dec, t.mi, t.psfmag_sdss[4] as
    psfmag_i,t.pmsig,t.ps1_var_sn[1],t.ps1_var_rms[1],t.des_var_sn[1],t.des_var_rms[1]
    from bhm_rm_v0 as t where (t.gaia = 1 AND t.mi < 21.5 AND t.pmsig
    < 5.0 AND t.plxsig < 5.0 AND t.ps1_var_sn[1] > 5.0 AND
    t.ps1_var_rms[1] > 0.05 ) limit 10;

    '''

    name = 'bhm_rm_var'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (
                (
                    (t.des_var_sn[0] > self.parameters['des_var_sn_min']) &
                    (t.des_var_rms[0] > self.parameters['des_var_rms_min'])
                ) |
                (
                    (t.ps1_var_sn[0] > self.parameters['ps1_var_sn_min']) &
                    (t.ps1_var_rms[0] > self.parameters['ps1_var_rms_min'])
                )
            ) &
            (t.source_id_gaia > 0) &
            (t.pmsig < self.parameters['pmsig_max']) &
            (t.plxsig < self.parameters['plxsig_max']) &
            (t.pmsig > pmsig_min) &  # this catches cases with NULL=-9
            (t.plxsig > plxsig_min) &
            ~(t.field_name.contains('SDSS-RM'))  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmAncillaryCarton(BhmRmBaseCarton):
    '''
    bhm_rm_ancillary: from the Gaia_unWISE AGN catalog or the XDQSO catalog,
                      but requiring no proper motion/parallax detection from Gaia DR2

    SELECT * FROM bhm_rm
        WHERE photo_bitmask & (2^0+2^1) != 0
        AND WHERE mi BETWEEN 15.0 AND 21.5    ← TBD
        AND WHERE pmsig < 5.0
        AND WHERE plxsig < 5.0
    '''

    name = 'bhm_rm_ancillary'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t

        query = query.where(
            (t.photo_bitmask.bin_and(self.parameters['photo_bitmask']) != 0) &
            (t.pmsig < self.parameters['pmsig_max']) &
            (t.plxsig < self.parameters['plxsig_max']) &
            (t.pmsig > pmsig_min) &  # this catches cases with NULL=-9
            (t.plxsig > plxsig_min) &
            ~(t.field_name.contains('SDSS-RM'))  # ignore this carton in the SDSS-RM field
        )

        return query
