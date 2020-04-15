#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-31
# @Filename: bhm_rm.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

import peewee
import sdssdb

from sdssdb.peewee.sdss5db import catalogdb

from target_selection.cartons.base import BaseCarton
from target_selection.cartons.skymask import SkyMask
import pkg_resources


# Carton descriptions here:
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ReverberationMapping(RM)

## This module provides the following BHM cartons:
#  BHM_RM_CORE
#  BHM_RM_KNOWN_SPEC
#  BHM_RM_VAR
#  BHM_RM_ANCILLARY
gaia_epoch = 2015.5


class BhmRmCoreCarton(BaseCarton):

    name = 'bhm_rm_core'
    category = 'science'
    survey = 'BHM'
    cadence = 'bhm_aqmes_rm_174x8'
    tile = False
    priority = 1

    def build_query(self):
        '''
        bhm_rm_core: select all photometric QSO targets with the likelihood method (Skewt), flux-limited to 21.5 in i-band PSF mag

        SELECT * FROM bhm_rm
                WHERE skewt_qso = 1
                AND WHERE mi BETWEEN 15.0 AND 21.5
                AND WHERE pmsig < 5.0
                AND WHERE plxsig < 5.0
        '''
        t = catalogdb.BHM_RM_V0.alias()

        query = (
            t.select(t.pk.alias('catalog_id'),
                     t.ra,
                     t.dec,
                     Value(0.0).alias('pmra'),
                     Value(0.0).alias('pmdec'),
                     Value(gaia_epoch).alias('epoch'),
                     t.mi.alias('mag_i'),
                     Value(self.priority).alias(priority))
            .where((t.skewt_qso == 1) &
                   (t.mi <  self.config['i_mag_max']) &
                   (t.mi >  self.config['i_mag_min']) &
                   (t.pmsig <  self.config['pmsig_max']) &
                   (t.plxsig <  self.config['plxsig_min'])))

        print(f"This query will return nrows={query.count()}")

        return query


class BhmRmKnownSpecCarton(BaseCarton):

    name = 'bhm_rm_known_spec'
    category = 'science'
    survey = 'BHM'
    cadence = 'bhm_aqmes_rm_174x8'
    tile = False
    priority = 1

    def build_query(self):
        '''
        bhm_rm_known_spec:  select all spectroscopically confirmed QSOs where redshift is extragalactic

        SELECT * FROM bhm_rm
            WHERE specz > 0.005
            AND WHERE mi BETWEEN 15.0 AND 21.5   # <- exact limits TBD

        '''
        t = catalogdb.BHM_RM_V0.alias()

        query = (
            t.select(t.pk.alias('catalog_id'),
                     t.ra,
                     t.dec,
                     Value(0.0).alias('pmra'),
                     Value(0.0).alias('pmdec'),
                     Value(gaia_epoch).alias('epoch'),
                     t.mi.alias('mag_i'),
                     Value(self.priority).alias(priority))
            .where((t.specz > self.config['specz_min']) &
                   (t.mi <  self.config['i_mag_max']) &
                   (t.mi >  self.config['i_mag_min'])))

        print(f"This query will return nrows={query.count()}")

        return query




class BhmRmVar(BaseCarton):

    name = 'bhm_rm_var'
    category = 'science'
    survey = 'BHM'
    cadence = 'bhm_aqmes_rm_174x8'
    tile = False
    priority = 1

    def build_query(self):
        '''
        bhm_rm_var: selected based on g-band variability > 0.05 mag and bright enough to be detected by Gaia (G<~21)

        SELECT * FROM bhm_rm
            WHERE ( (des_var_sn[0] > 5.0 AND des_var_rms[0] > 0.05) OR (ps1_var_sn[0]>5.0 AND ps1_var_rms[0]>0.05))
            AND WHERE mi BETWEEN 15.0 AND 21.5    # <- exact limits TBD
            AND WHERE pmsig < 5.0
            AND WHERE plxsig < 5.0
            AND WHERE gaia = 1

        '''
        t = catalogdb.BHM_RM_V0.alias()

        query = (
            t.select(t.pk.alias('catalog_id'),
                     t.ra,
                     t.dec,
                     Value(0.0).alias('pmra'),
                     Value(0.0).alias('pmdec'),
                     Value(gaia_epoch).alias('epoch'),
                     t.mi.alias('mag_i'),
                     Value(self.priority).alias(priority))
            .where(
                (
                    (
                        (t.des_var_sn[0] > self.config['des_var_sn_min']) &
                        (t.des_var_rms[0] > self.config['des_var_rms_min'])
                    ) |
                (
                    (t.ps1_var_sn[0] > self.config['ps1_var_sn_min']) &
                    (t.ps1_var_rms[0] > self.config['ps1_var_rms_min'])
                )
                ) &
                (t.mi <  self.config['i_mag_max']) &
                (t.mi >  self.config['i_mag_min']) &
                # TBD                (t.mg <  self.config['g_mag_max']) &
                (t.gaia == 1) &
                (t.pmsig <  self.config['pmsig_max']) &
                (t.plxsig <  self.config['plxsig_min'])
            )
        )

        print(f"This query will return nrows={query.count()}")

        return query



class BhmRmAncillaryCarton(BaseCarton):

    name = 'bhm_rm_ancillary'
    category = 'science'
    survey = 'BHM'
    cadence = 'bhm_aqmes_rm_174x8'
    tile = False
    priority = 1

    def build_query(self):
        '''
        bhm_rm_ancillary: from the Gaia_unWISE AGN catalog or the XDQSO catalog,
                          but requiring no proper motion/parallax detection from Gaia DR2

        SELECT * FROM bhm_rm
            WHERE photo_bitmask & (2^0+2^1) != 0
            AND WHERE mi BETWEEN 15.0 AND 21.5    ← TBD
            AND WHERE pmsig < 5.0
            AND WHERE plxsig < 5.0
        '''
        t = catalogdb.BHM_RM_V0.alias()

        query = (
            t.select(t.pk.alias('catalog_id'),
                     t.ra,
                     t.dec,
                     Value(0.0).alias('pmra'),
                     Value(0.0).alias('pmdec'),
                     Value(gaia_epoch).alias('epoch'),
                     t.mi.alias('mag_i'),
                     Value(self.priority).alias(priority))
            .where((t.photo_bitmask.bin_and(self.config['photo_bitmask']) != 0 ) &
                   (t.mi <  self.config['i_mag_max']) &
                   (t.mi >  self.config['i_mag_min']) &
                   (t.pmsig <  self.config['pmsig_max']) &
                   (t.plxsig <  self.config['plxsig_min'])
            )
        )

        print(f"This query will return nrows={query.count()}")

        return query
