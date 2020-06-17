#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_gua.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import peewee
import sdssdb
from peewee import JOIN

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToTIC_v8,
                                             TIC_v8,
                                             Gaia_DR2,
####                                             CatalogToGaia_unWISE_AGN,
#TODO                                             CatalogToSDSS_DR16_SpecObj,
                                             SDSS_DR16_SpecObj,
                                             Gaia_unWISE_AGN)



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
        AND WHERE (so.specobjid = NULL OR so.zwarning != 0 OR so.sn_median_all < x.x OR so.z_err > 0.0xx )

    bhm_gaia_unwise_agn_dark
        AND WHERE ( gua.g > 16.5 AND gua.rp > 16.5 AND (gua.g < 21.2 OR gua.rp < 21.5 ) )

    bhm_gaia_unwise_agn_bright
        AND WHERE ( gua.g > 13.5 AND gua.rp > 13.5 AND (gua.g < 18.0 OR gua.rp < 18.0) )
'''



class BhmGuaBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for both Gaia UnWISE AGN cartons
    To be sub-classed, not to be called directly.

    To get from Catalog to GUA we join tables via :
    Catalog -> CatalogToTIC_v8 -> Gaia_DR2 -> Gaia_unWISE_AGN
    '''

    name = 'bhm_gaia_unwise_agn_base'
    category = 'science'
    mapper = 'BHM'
    program = 'BHM-GUA'
    tile = False
    priority = None

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2tic = CatalogToTIC_v8.alias()
        tic = TIC_v8.alias()
#TODO        c2s = CatalogToSDSS_DR16_SpecObj.alias()
        g = Gaia_DR2.alias()
        t = Gaia_unWISE_AGN.alias()
        s = SDSS_DR16_SpecObj.alias()

        # set the Carton priority+values here - read from yaml
        target_priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')
        target_value = peewee.Value(float(self.parameters.get('value', 1.0))).alias('value')
        pmra = peewee.Value(0.0).alias('pmra')
        pmdec = peewee.Value(0.0).alias('pmdec')
        match_radius_spectro = self.parameters['spec_join_radius']/3600.0

        query = (
            c.select(c.catalogid,
                     target_priority,
                     pmra,
                     pmdec,
                     target_value,
            )
            .join(c2tic)
            .join(tic)
            .join(g)
            .join(t, on=(g.source_id == t.gaia_sourceid))
            .join(s, JOIN.LEFT_OUTER,
                  on=(peewee.fn.q3c_join(c.ra,c.dec,
                                         s.ra,s.dec,
                                         match_radius_spectro) &
                      (s.snmedian >= self.parameters['spec_sn_thresh']) &
                      (s.zwarning == 0) &
                      (s.zerr <= self.parameters['spec_z_err_thresh']) &
                      (s.zerr > 0.0) &
                      (s.scienceprimary > 0)
                      )
            )
            .where(c.version_id == version_id,
                   c2tic.version_id == version_id,
                   c2tic.best == True)
            .where(
                (t.prob_rf >= self.parameters['prob_rf_min']),
                (t.g >= self.parameters['mag_g_min']),
                (t.rp >= self.parameters['mag_rp_min']),
                (
                    (t.g < self.parameters['mag_g_max']) |
                    (t.rp < self.parameters['mag_rp_max'])
                ),
                (s.specobjid.is_null()),
            )
#TODO            .switch(c)
#TODO            .join(c2s)
#TODO            .join(s)
#TODO            .where(
#TODO                (s.specobjid.is_null()) |
#TODO                (s.zwarning != 0 ) |
#TODO                (s.sn_median_all < self.parameters['spec_sn_thresh']) |
#TODO                (s.z_err >  self.parameters['spec_z_err_thresh'])
#TODO            )
            .distinct([t.gaia_sourceid])   # avoid duplicates - trust the GUA parent sample
        )

        return query


#-------  bhm_gaia_unwise_agn_dark   ------ #

class BhmGuaDarkCarton(BhmGuaBaseCarton):
    '''
        AND WHERE ( gua.g > 16.x AND gua.rp > 16.x)
    '''
    name = 'bhm_gaia_unwise_agn_dark'
    cadence = 'bhm_spiders_1x4'

#-------  bhm_gaia_unwise_agn_bright   ------ #

class BhmGuaBrightCarton(BhmGuaBaseCarton):
    '''
        AND WHERE ( gua.g < 18.x OR gua.rp < 18.x)
    '''
    name = 'bhm_gaia_unwise_agn_bright'
    cadence = 'bhm_boss_bright_3x1'
