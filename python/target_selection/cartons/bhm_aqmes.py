#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_aqmes.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn
from astropy.io import fits
import pkg_resources


from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    SDSS_DR16_SpecObj,
    SDSS_DR16_QSO,
    CatalogToSDSS_DR16_SpecObj,
    SDSSV_Plateholes,
    SDSSV_Plateholes_Meta,
)


from target_selection.cartons.base import BaseCarton

# this should probably live in a better place
radius_apo = 1.49  # degrees

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-AQMES-medium-inprogress   # noqa: E501

# This provides the following BHM cartons:
#  bhm_aqmes_med
#  bhm_aqmes_med-faint
#  bhm_aqmes_wide2
#  bhm_aqmes_wide2-faint
#  # bhm_aqmes_wide3        # dumped in v0.5
#  # bhm_aqmes_wide3-faint  # dumped in v0.5
#  bhm_aqmes_bonus-dark
#  bhm_aqmes_bonus-bright

# how do we relate the cadence names in v0.5 to cadence names in v0?

cadence_map_v0p5_to_v0 = {
    'dark_10x4': 'bhm_aqmes_medium_10x4',
    'dark_10x4_4yr': 'bhm_aqmes_medium_10x4',
    # 'dark_3x4': 'bhm_aqmes_wide_3x4',
    'dark_2x4': 'bhm_aqmes_wide_2x4',
    'dark_1x4': 'bhm_spiders_1x4',
    'bright_3x1': 'bhm_boss_bright_3x1',
}


class BhmAqmesBaseCarton(BaseCarton):

    '''
    Parent class that provides the underlying selections for all AQMES cartons
    '''

    name = 'bhm_aqmes_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_aqmes'
    instrument = 'BOSS'
    inertial = True
    tile = False
    priority = None
    alias_c = None
    alias_t = None
    alias_c2s = None
    cadence = None
    cadence_v0p5 = None

    # read the AQMES field centres from a fits file and convert to a list of dicts
    def get_fieldlist(self):
        stub = self.parameters.get('fieldlist', None)
        if stub is None or stub == '' or stub == 'None':
            return None

        # filename = pkg_resources.resource_filename( __name__, stub)
        filename = pkg_resources.resource_filename('target_selection', stub)
        assert len(filename) > 0

        try:
            hdul = fits.open(filename)
        except BaseException:
            raise Exception(f"Failed to find/open fieldlist file: {filename}")

        assert len(hdul[1].data) > 0

        # choose the correct subset of fields based on the cadence name and form a list of dicts
        # we have to use the v0 cadence names though
        assert self.cadence_v0p5 in cadence_map_v0p5_to_v0
        cadence_v0 = cadence_map_v0p5_to_v0[self.cadence_v0p5]

        try:
            fieldlist = [
                {'racen': r['RACEN'],
                 'deccen': r['DECCEN'],
                 'radius': radius_apo, }
                for r in hdul[1].data
                if r['CADENCE'] == cadence_v0
            ]
        except BaseException:
            raise Exception(f"Error interpreting contents of fieldlist file: {filename}")

        assert len(fieldlist) > 0

        return fieldlist

    def append_spatial_query(self, query, cte, fieldlist):
        '''Extend the peewee query using a list of field centres'''
        if fieldlist is None:
            return query
        elif len(fieldlist) == 0:
            return query

        q = False
        for f in fieldlist:
            q = (q | peewee.fn.q3c_radial_query(cte.c.ra,
                                                cte.c.dec,
                                                f['racen'],
                                                f['deccen'],
                                                f['radius']))
        return query.where(q)

    # main query
    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2s = CatalogToSDSS_DR16_SpecObj.alias()
        s = SDSS_DR16_SpecObj.alias()
        t = SDSS_DR16_QSO.alias()
        self.alias_c = c
        self.alias_t = t
        self.alias_c2s = c2s

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped and that have firstcarton ~ 'bhm_aqmes_'
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
                ssph.firstcarton.contains('bhm_aqmes_'),
                ssphm.isvalid > 0,
            )
            .distinct([ssph.catalogid])
            .alias('sph')
        )

        # set the Carton priority+values here - read from yaml
        priority_floor = peewee.Value(int(self.parameters.get('priority', 999999)))
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast('bool')
        opt_prov = peewee.Value('sdss_psfmag')
        cadence_v0 = peewee.Value(cadence_map_v0p5_to_v0[self.cadence_v0p5]).cast('text')
        # cadence = peewee.Value(cadence_v0)
        cadence = peewee.Value(self.cadence_v0p5).cast('text')

        # # this is DEBUG until the new v0.5 cadences exist in the DB
        # # - doesn't work because self.cadence is checked before this point
        # # - so give up until targetdb.cadence is populated
        # assert self.cadence in cadence_map_v0p5_to_v0
        # v0_cadence = cadence_map_v0p5_to_v0[self.cadence]
        # cadence = peewee.Value(v0_cadence).alias('cadence')

        match_radius_spectro = 1.0 / 3600.0

        priority_boost = peewee.Case(
            None,
            (
                (sph.c.pkey.is_null(False), 0),  # has a platehole entry
                (sph.c.pkey.is_null(True), 1),   # not in plate programme
            ),
            None
        )
        priority = priority_floor + priority_boost

        magnitude_sdss_g = peewee.Case(
            None, ((t.psfmag[1].between(0.1, 29.9), t.psfmag[1]),), 'NaN').cast('float')
        magnitude_sdss_r = peewee.Case(
            None, ((t.psfmag[2].between(0.1, 29.9), t.psfmag[2]),), 'NaN').cast('float')
        magnitude_sdss_i = peewee.Case(
            None, ((t.psfmag[3].between(0.1, 29.9), t.psfmag[3]),), 'NaN').cast('float')
        magnitude_sdss_z = peewee.Case(
            None, ((t.psfmag[4].between(0.1, 29.9), t.psfmag[4]),), 'NaN').cast('float')
        magnitude_gaia_g = peewee.Case(
            None, ((t.gaia_g_mag.between(0.1, 29.9), t.gaia_g_mag),), 'NaN').cast('float')
        magnitude_gaia_bp = peewee.Case(
            None, ((t.gaia_bp_mag.between(0.1, 29.9), t.gaia_bp_mag),), 'NaN').cast('float')
        magnitude_gaia_rp = peewee.Case(
            None, ((t.gaia_rp_mag.between(0.1, 29.9), t.gaia_rp_mag),), 'NaN').cast('float')

        bquery = (
            c.select(
                c.catalogid,
                t.pk.alias('dr16q_pk'),  # extra
                s.specobjid.cast('text').alias('dr16_specobjid'),  # extra
                c.ra,   # extra
                c.dec,   # extra
                priority.alias('priority'),
                value.alias('value'),
                inertial.alias('inertial'),
                instrument.alias('instrument'),
                cadence.alias('cadence'),
                cadence_v0.alias('cadence_v0'),
                opt_prov.alias('optical_prov'),
                magnitude_sdss_g.alias('g'),
                magnitude_sdss_r.alias('r'),
                magnitude_sdss_i.alias('i'),
                magnitude_sdss_z.alias('z'),
                magnitude_gaia_g.alias('gaia_g'),
                magnitude_gaia_bp.alias('bp'),
                magnitude_gaia_rp.alias('rp'),
                t.plate.alias('dr16q_plate'),   # extra
                t.mjd.alias('dr16q_mjd'),   # extra
                t.fiberid.alias('dr16q_fiberid'),   # extra
                t.ra.alias("dr16q_ra"),   # extra
                t.dec.alias("dr16q_dec"),   # extra
                t.gaia_ra.alias("dr16q_gaia_ra"),   # extra
                t.gaia_dec.alias("dr16q_gaia_dec"),   # extra
                t.sdss2gaia_sep.alias("dr16q_sdss2gaia_sep"),   # extra
                t.z.alias("dr16q_redshift"),   # extra
                c2s.best.alias("c2s_best"),  # extra
            )
            .join(c2s)
            .join(s)
            .join(
                t,
                on=((s.plate == t.plate) &
                    (s.mjd == t.mjd) &
                    (s.fiberid == t.fiberid))
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
                c.version_id == version_id,
                c2s.version_id == version_id,
                # c2s.best >> True,   # TODO check this is working in v0.5
                #                     # - this condition killed many AQMES
                #                     #   targets in v0 cross-match
            )
            .where
            (
                t.psfmag[3] >= self.parameters['mag_i_min'],
                t.psfmag[3] < self.parameters['mag_i_max'],
                # (t.z >= self.parameters['redshift_min']), # not needed
                # (t.z <= self.parameters['redshift_max']),
            )
            # .distinct([t.pk])   # avoid duplicates - trust the QSO parent sample
            .distinct([c.catalogid])   # avoid duplicates - trust the catalog
            .cte('bquery', materialized=True)
        )

        query = bquery.select(peewee.SQL('bquery.*'))
        query = self.append_spatial_query(query, bquery, self.get_fieldlist())
        query = query.with_cte(bquery)

        return query


# -------AQMES medium section ------ #


class BhmAqmesMedCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso
    WHERE  psfmag_i BETWEEN 16.x AND 19.1
    AND   {target lies in spatial selection}
    '''
    name = 'bhm_aqmes_med'
    cadence_v0p5 = 'dark_10x4_4yr'

    # TD's note to self:
    # add something like the following if want to add carton-specific selections
    #    def build_query(self, version_id, query_region=None):
    #        query = super().build_query(version_id, query_region)
    #        query = query.where( # .... add extra terms here
    #        )
    #        return query


class BhmAqmesMedFaintCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso
    WHERE  psfmag_i BETWEEN 19.1 AND 21.0
    AND   {target lies in spatial selection}
    '''
    name = 'bhm_aqmes_med_faint'
    cadence_v0p5 = 'dark_10x4_4yr'
    program = 'bhm_filler'

# -------AQMES medium section ----- #
#
#

#
# -------AQMES wide section ------ #


# class BhmAqmesWide3Carton(BhmAqmesBaseCarton):
#     '''
#     SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 16.0 AND 19.1
#     '''
#     name = 'bhm_aqmes_wide3'
#     cadence_v0p5 = 'dark_3x4'
#
#
# class BhmAqmesWide3FaintCarton(BhmAqmesBaseCarton):
#     '''
#     SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 19.1 AND 21.0
#     '''
#     name = 'bhm_aqmes_wide3_faint'
#     cadence_v0p5 = 'dark_3x4'
#     program = 'bhm_filler'


class BhmAqmesWide2Carton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 16.0 AND 19.1
    '''
    name = 'bhm_aqmes_wide2'
    cadence_v0p5 = 'dark_2x4'


class BhmAqmesWide2FaintCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 19.1 AND 21.0
    '''
    name = 'bhm_aqmes_wide2_faint'
    cadence_v0p5 = 'dark_2x4'
    program = 'bhm_filler'

# -------AQMES wide section ------ #


# -------AQMES bonus section ------ #

class BhmAqmesBonusCoreCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 16.0 AND 19.1
    {NO spatial constraint}
    '''
    name = 'bhm_aqmes_bonus_core'
    cadence_v0p5 = 'dark_1x4'
    program = 'bhm_filler'


class BhmAqmesBonusFaintCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 19.1 AND 21.0
    '''
    name = 'bhm_aqmes_bonus_faint'
    cadence_v0p5 = 'dark_1x4'
    program = 'bhm_filler'


class BhmAqmesBonusBrightCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr16_qso WHERE psfmag_i BETWEEN 14.0 AND 18.0
    '''
    name = 'bhm_aqmes_bonus_bright'
    cadence_v0p5 = 'bright_3x1'
    program = 'bhm_filler'

# ------- AQMES bonus section ------ #
