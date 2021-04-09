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
# import sdssdb

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             BHM_CSC,
                                             CatalogToBHM_CSC)

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

from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ChandraSourceCatalogue(CSC)-inprogress  # noqa: E501

# This module provides the following BHM cartons:
#   bhm_csc_boss_dark
#   bhm_csc_boss_bright
#   bhm_csc_apogee


class BhmCscBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for all CSC cartons
    To be sub-classed, not typically to be called directly.

    SELECT DISTINCT ON (t.pk)
      c.catalogid, c.lead, c.ra, c.dec,
      q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
      t.mag_g, t.mag_r, t.mag_i
    FROM catalog AS c
    JOIN catalog_to_bhm_csc AS c2t
      ON c.catalogid = c2t.catalogid
    JOIN bhm_csc AS t
      ON c2t.target_id = t.pk
    WHERE
          c.version_id = 21
      AND c2t.version_id = 21
      AND c2t.best IS TRUE
      AND t.spectrograph = 'BOSS'
      AND t.mag_i BETWEEN 16.0 AND 22.0
    ORDER BY t.pk
    LIMIT 20;


    # check how many csc targets have positions defined by the  bhm_csc.oir_ra,oir_dec fields
    # these may actually be Chandra X-ray coords and not optical positions!

    SELECT
        count(*) as count,
        AVG(sep) as mean_sep,
        MIN(sep) as min_sep,
        MAX(sep) as max_sep,
        AVG(mag_r) as mean_mag_r,
        MIN(mag_r) as min_mag_r,
        MAX(mag_r) as max_mag_r
    FROM (
        SELECT DISTINCT ON (t.pk)
           c.catalogid, c.lead, c.ra, c.dec,
          q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
          t.mag_g, t.mag_r, t.mag_i
        FROM catalog AS c
        JOIN catalog_to_bhm_csc AS c2t
          ON c.catalogid = c2t.catalogid
        JOIN bhm_csc AS t
          ON c2t.target_id = t.pk
        WHERE
              c.version_id = 21
          AND c2t.version_id = 21
          AND c2t.best IS TRUE
          AND t.spectrograph = 'BOSS'
          AND t.mag_i BETWEEN 16.0 AND 22.0
        ORDER BY t.pk
        ) AS subq
      WHERE subq.lead = 'bhm_csc';

    #
    SELECT DISTINCT ON (t.pk)
          c.catalogid, c.lead, c.ra, c.dec,
          q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
          t.mag_g, t.mag_r, t.mag_i
        INTO sandbox.temp_td_csc_subset
        FROM catalog AS c
        JOIN catalog_to_bhm_csc AS c2t
          ON c.catalogid = c2t.catalogid
        JOIN bhm_csc AS t
          ON c2t.target_id = t.pk
        WHERE
              c.version_id = 21
          AND c2t.version_id = 21
          AND c2t.best IS TRUE
          AND t.spectrograph = 'BOSS'
          AND t.mag_i BETWEEN 16.0 AND 22.0
          AND c.lead = 'bhm_csc'
        ORDER BY t.pk;



    '''

    name = 'bhm_csc_base'
    base_name = 'bhm_csc_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_csc'
    tile = False
    priority = None
    alias_t = None
    instrument = None
    this_cadence = None

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2t = CatalogToBHM_CSC.alias()
        t = BHM_CSC.alias()
        self.alias_t = t
        # c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        # s16 = SDSS_DR16_SpecObj.alias()
        # s2020 = BHM_eFEDS_Veto.alias()
        # sV = SDSSV_BOSS_SPALL.alias()
        # ph = SDSSV_Plateholes.alias()
        # phm = SDSSV_Plateholes_Meta.alias()

        # set the Carton priority+values here - read from yaml
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        instrument = peewee.Value(self.instrument)
        cadence = peewee.Value(self.this_cadence)
        opt_prov = peewee.Value('ps1_psfmag')

        if (self.instrument == 'BOSS'):

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
            priority = (
                peewee.Value(self.parameters['priority_floor']) +
                priority_1 * dpriority_has_spec
            )
        else:
            priority = peewee.Value(self.parameters['priority_floor'])

        # Process the bhm_csc.[g,r,i,z,h] magnitudes to deal with zeros
        magnitude_g = peewee.Case(None, ((t.mag_g <= 0.0, None),), t.mag_g).cast('float')
        magnitude_r = peewee.Case(None, ((t.mag_r <= 0.0, None),), t.mag_r).cast('float')
        magnitude_i = peewee.Case(None, ((t.mag_i <= 0.0, None),), t.mag_i).cast('float')
        magnitude_z = peewee.Case(None, ((t.mag_z <= 0.0, None),), t.mag_z).cast('float')
        magnitude_h = peewee.Case(None, ((t.mag_h <= 0.0, None),), t.mag_h).cast('float')

        query = (
            c.select(
                c.catalogid,
                t.cxo_name,   # extra
                t.pk.alias('csc_pk'),   # extra
                c.ra,  # extra
                c.dec,  # extra
                priority.alias('priority'),
                value.alias('value'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                magnitude_h.alias('h'),
                opt_prov.alias('opt_prov'),
                t.oir_ra.alias('csc_oir_ra'),
                t.oir_dec.alias('csc_oir_dec'),
            )
            .join(c2t)
            .join(t)
            .where(
                c.version_id == version_id,
                c2t.version_id == version_id,
                c2t.best >> True
            )
            # .distinct([c2t.target_id])  # avoid duplicates - trust the CSC parent sample,
            .distinct([c.catalogid])  # avoid duplicates - trust the catalogid,
            .where
            (
                t.spectrograph == self.instrument
            )
        )

        if (self.instrument == 'BOSS'):
            # Append the spectro query
            query = (
                query
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
            )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query


class BhmCscBossDarkCarton(BhmCscBaseCarton):
    '''
    SELECT * from bhm_csc AS c
    WHERE c.spectrograph = "BOSS"
       ND c.mag_i BETWEEN 17.x AND 21.x
    '''
    name = 'bhm_csc_boss_dark'
    this_cadence = 'dark_1x4'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i < self.parameters['mag_i_max'])
        )
        return query


class BhmCscBossBrightCarton(BhmCscBaseCarton):
    '''
    SELECT * from bhm_csc AS c
    WHERE c.spectrograph = "BOSS"
      AND c.mag_i BETWEEN 1x.0 AND 18.x
    '''
    name = 'bhm_csc_boss_bright'
    this_cadence = 'bright_1x1'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i < self.parameters['mag_i_max'])
        )
        return query


class BhmCscApogeeCarton(BhmCscBaseCarton):
    '''
    SELECT * from bhm_csc AS c
    WHERE c.spectrograph = "APOGEE"
      AND c.mag_H BETWEEN 1x.x AND 1x.x
    '''
    name = 'bhm_csc_apogee'
    this_cadence = 'bright_3x1'
    instrument = 'APOGEE'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_h >= self.parameters['mag_h_min']) &
            (t.mag_h < self.parameters['mag_h_max'])
        )
        return query

#
# '''
# Exporting from the temp table
#
# \copy (SELECT * FROM sandbox.temp_bhm_csc_boss_dark)  TO 'bhm_csc_boss_dark.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_csc_boss_bright)
#      TO 'bhm_csc_boss_bright.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_csc_apogee)  TO 'bhm_csc_apogee.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_td_csc_subset)  TO 'temp_td_csc_subset.csv' with csv header
# for F in bhm_csc_*.csv; do
#    stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic;
# done
#
# '''
