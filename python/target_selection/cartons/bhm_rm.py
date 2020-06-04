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

sdssdb.peewee.sdss5db.database.set_profile('tunnel_operations')

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, BHM_RM_v0, CatalogToBHM_RM_v0)

from target_selection.cartons.base import BaseCarton
#from target_selection.cartons.skymask import SkyMask
#import pkg_resources


# Carton descriptions here:
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ReverberationMapping(RM)

## This module provides the following BHM cartons:
#  BHM_RM_CORE
#  BHM_RM_KNOWN_SPEC
#  BHM_RM_VAR
#  BHM_RM_ANCILLARY
#gaia_epoch = 2015.5

'''
t = BHM_RM_v0.alias()
for f in t._meta.fields:
    print (f)
'''

# sql, param = qq.sql()
# print (sql.replace("%s","{}").format(*param))


class BhmRmBaseCarton(BaseCarton):
    '''
    This class provides common setting and the masking routines used by all RM cartons
    '''

    name = 'bhm_rm_base'
    base_name = 'bhm_rm_base'
    category = 'science'
    mapper = 'BHM'
    cadence = 'bhm_rm_174x8'
    tile = False
    priority = 1000
    alias_t = None

    # form a peewee sub-query using the RM field centres
    def append_query_region(self, cat_alias, query):
        base_parameters = self.config['parameters'].get(self.base_name, None)
        q = False
        for field in base_parameters['fields']:
            q = ( q | peewee.fn.q3c_radial_query(cat_alias.ra,
                                                 cat_alias.dec,
                                                 field['racen'],
                                                 field['deccen'],
                                                 field['radius']))
        return query.where(q)


    def build_query(self):
        t = BHM_RM_v0.alias()
        self.alias_t = t
        c2t = CatalogToBHM_RM_v0.alias()

        query = (
            c
            .select(c.catalogid,
                    c.ra,
                    c.dec,
                    t.mi.alias('mag_i')
            )
            .join(c2t)
            .join(t)
            )
        query = self.append_query_region(c, query)

        return query







class BhmRmCoreCarton(BhmRmBaseCarton):
    '''
    bhm_rm_core: select all photometric QSO targets with the likelihood method (Skewt), flux-limited to 21.5 in i-band PSF mag

    SELECT * FROM bhm_rm
            WHERE skewt_qso = 1
            AND WHERE mi BETWEEN 15.0 AND 21.5
            AND WHERE pmsig < 5.0
            AND WHERE plxsig < 5.0
    '''

    name = 'bhm_rm_core'
    priority = 1002

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
#        c = Catalog.alias()
#        t = BHM_RM_v0.alias()
#        c2t = CatalogToBHM_RM_v0.alias()

#        query = (
#            c
#            .select(c.catalogid,
#                    c.ra,
#                    c.dec,
#                    t.mi.alias('mag_i')
#            )
#            .join(c2t)
#            .join(t)
        query = query.where
        (
            (t.skewt_qso == 1),
            (t.mi <  self.parameters['i_mag_max']),
            (t.mi >  self.parameters['i_mag_min']),
            (t.pmsig <  self.parameters['pmsig_max']),
            (t.plxsig <  self.parameters['plxsig_min']),
        )

        print(f"This query will return nrows={query.count()}")

        return query


class BhmRmKnownSpecCarton(BhmRmBaseCarton):
    '''
    bhm_rm_known_spec:  select all spectroscopically confirmed QSOs where redshift is extragalactic

    # recommended by q3c notes:
    set enable_mergejoin to off;
    set enable_seqscan to off;

    SELECT * FROM bhm_rm
        WHERE specz > 0.005
        AND WHERE mi BETWEEN 15.0 AND 21.5   # <- exact limits TBD

    # functional SQL for one field:
    ##SELECT c.*,t.specz,t.mi FROM catalog AS C INNER JOIN catalog_to_bhm_rm_v0 AS c2t ON c.catalogid = c2t.catalogid INNER JOIN bhm_rm_v0 AS t ON c2t.target_id = t.pk WHERE (q3c_radial_query(c.ra, c.dec, 36.45, -4.6, 1.49) AND c.version_id = 13 AND t.specz > 0.005 AND mi < 21.1 AND mi > 16.0);
    SELECT DISTINCT ON (t.pk) c.*,t.mg,t.mi,t.specz,t.spec_q,t.spec_strmask,specz_ref,q3c_dist(c.ra,c.dec,t.ra,t.dec)*3600.::FLOAT as sep FROM catalog AS C INNER JOIN catalog_to_bhm_rm_v0 AS c2t ON c.catalogid = c2t.catalogid INNER JOIN bhm_rm_v0 AS t ON c2t.target_id = t.pk WHERE (q3c_radial_query(c.ra, c.dec, 36.45, -4.6, 1.49) AND c.version_id = 13 AND t.specz > 0.005 AND mi < 21.1 AND mi > 16.0);
    \copy (SELECT DISTINCT ON (t.pk) c.*,t.mg,t.mi,t.specz,t.spec_q,t.spec_strmask,specz_ref,q3c_dist(c.ra,c.dec,t.ra,t.dec)*3600.::FLOAT as sep FROM catalog AS C INNER JOIN catalog_to_bhm_rm_v0 AS c2t ON c.catalogid = c2t.catalogid INNER JOIN bhm_rm_v0 AS t ON c2t.target_id = t.pk WHERE (q3c_radial_query(c.ra, c.dec, 35.5, -4.6, 1.49) AND c.version_id = 13 AND t.specz > 0.005 AND mi < 22.0 AND mi > 16.0)) to ./test_faint.csv csv header;
    '''

    name = 'bhm_rm_known_spec'

    priority = 1001

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
#        c = Catalog.alias()
#        t = BHM_RM_v0.alias()
#        c2t = CatalogToBHM_RM_v0.alias()

#        query = (
#            c
#            .select(c.catalogid,
#                    c.ra,
#                    c.dec,
#                    t.mi.alias('mag_i')
#             )
#            .join(c2t)
#            .join(t)
        query = query.where
        (
            (t.specz > self.parameters['specz_min']),
            (t.mi <  self.parameters['i_mag_max']),
            (t.mi >  self.parameters['i_mag_min'])
        )

        print(f"This query will return nrows={query.count()}")

        return query




class BhmRmVar(BhmRmBaseCarton):
    '''
    bhm_rm_var: selected based on g-band variability > 0.05 mag and bright enough to be detected by Gaia (G<~21)

    SELECT * FROM bhm_rm
        WHERE ( (des_var_sn[0] > 5.0 AND des_var_rms[0] > 0.05) OR (ps1_var_sn[0]>5.0 AND ps1_var_rms[0]>0.05))
        AND WHERE mi BETWEEN 15.0 AND 21.5    # <- exact limits TBD
        AND WHERE pmsig < 5.0
        AND WHERE plxsig < 5.0
        AND WHERE gaia = 1

    '''

    name = 'bhm_rm_var'

    priority = 1003

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
#        c = Catalog.alias()
#        t = BHM_RM_v0.alias()
#        c2t = CatalogToBHM_RM_v0.alias()
#
#        query = (
#            c
#            .select(c.catalogid,
#                    c.ra,
#                    c.dec,
#                    t.mi.alias('mag_i')
#             )
#            .join(c2t)
#            .join(t)
        query = query.where
        (
            (
                (
                    (t.des_var_sn[0] > self.parameters['des_var_sn_min']) &
                    (t.des_var_rms[0] > self.parameters['des_var_rms_min'])
                ) |
                (
                    (t.ps1_var_sn[0] > self.parameters['ps1_var_sn_min']) &
                    (t.ps1_var_rms[0] > self.parameters['ps1_var_rms_min'])
                )
            ),
            (t.mi <  self.parameters['i_mag_max']),
            (t.mi >  self.parameters['i_mag_min']),
            (t.gaia == 1),
            (t.pmsig <  self.parameters['pmsig_max']),
            (t.plxsig <  self.parameters['plxsig_min']),
            # (t.mg <  self.parameters['g_mag_max']) &  # TBD
        )


        print(f"This query will return nrows={query.count()}")

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

    priority = 1004

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
#        c = Catalog.alias()
#        t = BHM_RM_v0.alias()
#        c2t = CatalogToBHM_RM_v0.alias()
#
#        query = (
#            c
#            .select(c.catalogid,
#                    c.ra,
#                    c.dec,
#                    t.mi.alias('mag_i')
#             )
#            .join(c2t)
#            .join(t)

        query = query.where
        (
            (t.photo_bitmask.bin_and(self.parameters['photo_bitmask']) != 0 ),
            (t.mi <  self.parameters['i_mag_max']),
            (t.mi >  self.parameters['i_mag_min']),
            (t.pmsig <  self.parameters['pmsig_max']),
            (t.plxsig <  self.parameters['plxsig_min']),
        )

        print(f"This query will return nrows={query.count()}")

        return query
