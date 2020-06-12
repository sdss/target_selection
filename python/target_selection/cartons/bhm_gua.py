#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_gua.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import peewee
import sdssdb
#from astropy.io import fits
#import pkg_resources

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToGaia_unWISE_AGN,
#                                             CatalogToSDSS_DR16_SpecObj,
                                             Gaia_unWISE_AGN,
                                             SDSS_DR16_SpecObj)


from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms

## This module provides the following BHM cartons:
# bhm_gaia_unwise_agn_dark
# bhm_gaia_unwise_agn_bright


'''
  pseudo-SQL

    both Cartons:
        SELECT * FROM gaia_unwise_agn AS gua
        LEFT JOIN sdss_specobj_dr16 AS so ON  q3c_join(gua.ra, gua.dec, so.plug_ra, so.plug_dec, 1.0)
        AND WHERE gua.prob_rf > 0.8
        AND WHERE (so.plug_ra = NULL OR so.zwarning != 0 OR so.sn_median_all < x.x OR so.z_err > 0.0xx )

    bhm_gaia_unwise_agn_dark
        AND WHERE ( gua.g > 16.5 AND gua.rp > 16.5)

    bhm_gaia_unwise_agn_dark
        AND WHERE (gua.g < 18.5 OR gua.rp < 18.5)
'''



class BhmGuaBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for both Gaia UnWISE AGN cartons
    To be sub-classed, not to be called directly.
    '''

    name = 'bhm_gua_base'
    category = 'science'
    mapper = 'BHM'
    program = 'BHM-GUA'
    tile = False
    priority = None
    alias_c = None
    alias_t = None

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2t = CatalogToGaia_unWISE_AGN.alias()
        t = CatalogToGaia_unWISE_AGN.alias()
        self.alias_c = c
        self.alias_t = t

        # set the Carton priority+values here - read from yaml
        target_priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')
        target_value = peewee.Value(float(self.parameters.get('value', 1.0))).alias('value')
        pmra = peewee.Value(0.0).alias('pmra')
        pmdec = peewee.Value(0.0).alias('pmdec')

        query = (
            c
            .select(c.catalogid,
                    target_priority,
                    pmra,
                    pmdec,
                    target_value,
                    t.g.alias('magnitude_g'),
            )
            .join(c2t)
            .join(t)
            .where(c.version_id == version_id,
                   c2t.version_id == version_id)
            .distinct([t.gaia_sourceid])   # avoid duplicates - trust the GUA parent sample
            .where
            (
                (t.g >= self.parameters['mag_g_min']),
                (t.g >= self.parameters['mag_g_min']),
            )
        )

        return query


#-------  bhm_gaia_unwise_agn_dark   ------ #

class BhmGuaDarkCarton(BhmGuaBaseCarton):
    '''
        AND WHERE ( gua.g > 16.5 AND gua.rp > 16.5)

    '''
    name = 'bhm_gaia_unwise_agn_dark'
    cadence = 'bhm_spiders_1x4'


#    def build_query(self, version_id, query_region=None):
#        query = super().build_query(version_id, query_region)
#        query = query.where(
#            (t.mag_i >= self.parameters['mag_i_min']) &
#            (t.mag_i <  self.parameters['mag_i_max']) &
#            (t.spectrograph = self.instrument_code)
#        )
#        return query

#-------bhm_csc_boss_bright ------ #

#-# class BhmCscBossBrightCarton(BhmCscBaseCarton):
#-#     '''
#-#         SELECT * from bhm_csc AS c WHERE c.spectrograph = "BOSS" AND WHERE c.mag_i BETWEEN 1x.0 AND 18.x
#-#     '''
#-#     name = 'bhm_csc_boss_bright'
#-#     cadence = 'bhm_csc_boss_1x1'
#-#     instrument_code = 'BOSS'
#-#
#-#     def build_query(self, version_id, query_region=None):
#-#         query = super().build_query(version_id, query_region)
#-#         query = query.where(
#-#             (t.mag_i >= self.parameters['mag_i_min']) &
#-#             (t.mag_i <  self.parameters['mag_i_max']) &
#-#             (t.spectrograph = self.instrument_code)
#-#         )
#-#         return query
#-#
#-# #-------bhm_csc_apogee ------ #
#-#
#-# class BhmCscApogeeCarton(BhmCscBaseCarton):
#-#     '''
#-#         SELECT * from bhm_csc AS c WHERE c.spectrograph = "APOGEE" AND WHERE c.mag_H BETWEEN 1x.x AND 1x.x
#-#     '''
#-#     name = 'bhm_csc_apogee'
#-#     cadence = 'bhm_csc_apogee_3x1'
#-#     instrument_code = 'APOGEE'
#-#
#-#     def build_query(self, version_id, query_region=None):
#-#         query = super().build_query(version_id, query_region)
#-#         query = query.where(
#-#             (t.mag_h >= self.parameters['mag_h_min']) &
#-#             (t.mag_h <  self.parameters['mag_h_max']) &
#-#             (t.spectrograph = self.instrument_code)
#-#         )
#-#         return query
#-#
