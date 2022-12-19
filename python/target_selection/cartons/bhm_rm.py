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
    BHM_RM_v1,
    CatalogToLegacy_Survey_DR8,
    CatalogToLegacy_Survey_DR10,
    CatalogToGaia_DR3,
)


#  This module provides the following BHM cartons in v05:
#  bhm_rm_core
#  bhm_rm_known_spec
#  bhm_rm_var
#  bhm_rm_ancillary

#pmsig_min = -3.0
#plxsig_min = -3.0


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
        c2ls8 = CatalogToLegacy_Survey_DR8.alias()
        c2ls10 = CatalogToLegacy_Survey_DR10.alias()
        c2g3 = CatalogToGaia_DR3.alias()
        t = BHM_RM_v1.alias()
        self.alias_c = c
        self.alias_t = t

        fieldlist = self.get_fieldlist()

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
                    (t.mag_i <= priority_mag_bright),
                    priority_floor + 0
                ),
                (
                    (
                        (self.name == 'bhm_rm_known_spec') &
                        ~(t.rm_field_name.contains('SDSS-RM')) &
                        (t.mag_i <= priority_mag_bright_known_spec)
                    ),
                    priority_floor + 0
                ),
                (
                    (t.mag_i <= priority_mag_faint),
                    priority_floor +
                    5 * (1 + peewee.fn.floor((t.mag_i - priority_mag_bright) /
                                             priority_mag_step).cast('int'))
                ),
                (
                    (t.mag_i > priority_mag_faint),
                    priority_floor + 95
                ),
            ),
            None
        )

        # combine the priorities
        priority = priority1

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast('bool')

        # This is the scheme used in v0
        cadence_v0 = peewee.Case(None,
                                 (
                                     (t.rm_field_name.contains('S-CVZ'), 'bhm_rm_lite5_100x8'),
                                 ),
                                 'bhm_rm_174x8')

        # this gives the new names for the same cadences assumed in v0
        cadence_v0p5 = peewee.Case(None,
                                   (
                                       (t.rm_field_name.contains('S-CVZ'), 'dark_100x8'),
                                   ),
                                   'dark_174x8')

        # the following will replace old generic cadences when relevant table has been populated
        # TODO - replace when correct cadences are loaded
        cadence_v1p0 = peewee.Case(None,
                                   (
                                       (t.rm_field_name.contains('SDSS-RM'), 'bhm_rm_sdss-rm'),
                                       (t.rm_field_name.contains('COSMOS'), 'bhm_rm_cosmos'),
                                       (t.rm_field_name.contains('XMM-LSS'), 'bhm_rm_xmm-lss'),
                                       (t.rm_field_name.contains('S-CVZ'), 'bhm_rm_cvz-s'),
                                       (t.rm_field_name.contains('CDFS'), 'bhm_rm_cdfs'),
                                       (t.rm_field_name.contains('ELIAS-S1'), 'bhm_rm_elias-s1'),
                                   ),
                                   'dark_174x8')

        opt_prov = peewee.Value("psfmag")

        query = (
            t.select(
                c.catalogid,
                c2ls10.catalogid.alias('c2ls10_catalogid'),  # extra
                c2g3.catalogid.alias('c2g3_catalogid'),  # extra
                c2ls8.catalogid.alias('c2ls8_catalogid'),  # extra
                c.ra,  # extra
                c.dec,  # extra
                t.rm_field_name.alias('rm_field_name'),  # extra
                t.pk.alias('rm_pk'),  # extra
                instrument.alias('instrument'),
                priority.alias('priority'),
                value.alias('value'),
                cadence_v0p5.alias('cadence'),
                cadence_v0.alias('cadence_v0'),  # extra
                cadence_v0p5.alias('cadence_v0p5'),  # extra
                cadence_v1p0.alias('cadence_v1p0'),  # extra
                t.mag_g.alias('g'),
                t.mag_r.alias('r'),
                t.mag_i.alias('i'),
                t.mag_z.alias('z'),
                t.gaia_g.alias('gaia_mag'),
                t.gaia_bp.alias('gaia_bp'),
                t.gaia_rp.alias('gaia_rp'),
                opt_prov.alias('optical_prov'),
                inertial.alias('inertial'),
                # c2t.best.alias("c2t_best"),  # extra
            )
            .join(c2ls10, JOIN.LEFT_OUTER,
                  on=(c2ls10.target_id == t.ls_id_dr10))
            .join(c2g3, JOIN.LEFT_OUTER,
                  on=(c2g3.target_id == t.gaia_dr3_source_id))
            .join(c2ls8, JOIN.LEFT_OUTER,
                  on=(c2ls8.target_id == t.ls_id_dr8))
            .join(c, on=(fn.coalesce(c2ls10.catalogid, c2g3.catalogid, c2ls8.catalogid)))
            .where(
                c.version_id == version_id,
                fn.coalesce(c2ls10.version_id, version_id) == version_id,
                fn.coalesce(c2g3.version_id, version_id) == version_id,
                fn.coalesce(c2ls8.version_id, version_id) == version_id,
                # fn.coalesce(c2ls10.best,True) >> True   # TODO check if this is dropping RM
                #                                         # targets like it does for AQMES
            )
            .where
            (
                (
                    (t.mag_i >= self.parameters['mag_i_min']) &
                    (t.mag_i < self.parameters['mag_i_max'])
                ) |
                (
                    # S-CVZ targets often have only Gaia photom
                    (t.rm_field_name.contains('S-CVZ')) &
                    (t.gaia_g >= self.parameters['mag_g_min_cvz_s']) &
                    (t.gaia_g < self.parameters['mag_g_max_cvz_s'])
                )
            )
            .where(
                # Reject any objects where the t.rm_unsuitable flag is set
                t.rm_unsuitable >> True,
            )
            # .distinct([t.pk])   # avoid duplicates - trust the RM parent sample
            # - only needed if NOT using c2t.best = True condition
        )
        query = self.append_spatial_query(query, fieldlist)

        return query


class BhmRmCoreCarton(BhmRmBaseCarton):

    name = 'bhm_rm_core'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.rm_core >> True),
            ~(t.rm_field_name.contains('SDSS-RM')),  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmKnownSpecCarton(BhmRmBaseCarton):
    '''
    bhm_rm_known_spec:  select all spectroscopically confirmed QSOs where redshift is extragalactic
    '''

    name = 'bhm_rm_known_spec'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (
                (t.rm_known_spec >> True),
            ),
            (
                ~(t.rm_field_name.contains('SDSS-RM')) |
                (
                    # include extra constraints on SDSS-RM targets
                    (t.mag_i < self.parameters['mag_i_max_sdss_rm'])
                )
            ),
            (
                ~(t.rm_field_name.contains('COSMOS')) |
                (
                    # include extra constraints on COSMOS targets
                    (t.mag_i < self.parameters['mag_i_max_cosmos'])
                )
            ),
            (
                ~(t.rm_field_name.contains('XMM-LSS')) |
                (
                    # include extra constraints on XMM-LSS targets
                    (t.mag_i < self.parameters['mag_i_max_xmm_lss'])
                )
            ),
        )

        return query


class BhmRmVarCarton(BhmRmBaseCarton):
    '''bhm_rm_var: selected based on g-band variability > 0.05 mag
                   and bright enough to be detected by Gaia (G<~21)
    '''

    name = 'bhm_rm_var'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.rm_var >> True),
            ~(t.rm_field_name.contains('SDSS-RM'))  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmAncillaryCarton(BhmRmBaseCarton):
    '''
    bhm_rm_ancillary: from the Gaia_unWISE AGN catalog or the XDQSO catalog,
                      but requiring no proper motion/parallax detection from Gaia DR2
    '''

    name = 'bhm_rm_ancillary'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t

        query = query.where(
            (t.rm_ancillary >> True),
            ~(t.rm_field_name.contains('SDSS-RM'))  # ignore this carton in the SDSS-RM field
        )

        return query
