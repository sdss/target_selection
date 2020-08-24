#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_gua.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# flake8: noqa
# isort: skip_file

import peewee
import sdssdb
from peewee import JOIN

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToTIC_v8,
                                             TIC_v8,
                                             Gaia_DR2,
#this is old and does not work               CatalogToGaia_unWISE_AGN,
                                             CatalogToSDSS_DR16_SpecObj,
                                             SDSS_DR16_SpecObj,
                                             Gaia_unWISE_AGN)



from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms

## This module provides the following BHM cartons:
# bhm_gua_dark
# bhm_gua_bright


'''
  pseudo-SQL

    both Cartons:
        SELECT * FROM gaia_unwise_agn AS gua
        LEFT JOIN sdss_dr16_specobj AS so ON  q3c_join(gua.ra, gua.dec, so.plug_ra, so.plug_dec, 1.0)
        AND WHERE gua.prob_rf > 0.8
        AND WHERE (so.specobjid = NULL OR so.zwarning != 0 OR so.sn_median_all < x.x OR so.z_err > 0.0xx )

    bhm_gua_dark
        AND WHERE ( gua.g > 16.5 AND gua.rp > 16.5 AND (gua.g < 21.2 OR gua.rp < 21.0 ) )

    bhm_gua_bright
        AND WHERE ( gua.g > 13.0 AND gua.rp > 13.5 AND (gua.g < 18.5 OR gua.rp < 18.5) )
'''



class BhmGuaBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for both Gaia UnWISE AGN cartons
    To be sub-classed, not to be called directly.

    To get from Catalog to GUA we join tables via :
    Catalog -> CatalogToTIC_v8 -> Gaia_DR2 -> Gaia_unWISE_AGN
    '''

    name = 'bhm_gua_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_filler'
    tile = False
    priority = None

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        ##############c2t = CatalogToGaia_unWISE_AGN.alias() - deprecated
        c2tic = CatalogToTIC_v8.alias()
        tic = TIC_v8.alias()
        c2s = CatalogToSDSS_DR16_SpecObj.alias()
        g = Gaia_DR2.alias()
        t = Gaia_unWISE_AGN.alias()
        s = SDSS_DR16_SpecObj.alias()

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float').alias('value')
        pmra = peewee.Value(0.0).cast('float').alias('pmra')
        pmdec = peewee.Value(0.0).cast('float').alias('pmdec')
        parallax = peewee.Value(0.0).cast('float').alias('parallax')

        query = (
            c.select(c.catalogid,
                     priority,
                     value,
                     pmra,
                     pmdec,
                     parallax,
#                     t.g.alias('g'),   ##rely on the centralised magnitude routines
#                     t.bp.alias('bp'),
#                     t.rp.alias('rp'),
            )
            .join(c2tic)
            .join(tic)
            .join(g)
            .join(t, on=(g.source_id == t.gaia_sourceid))
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
            )
            .switch(c)
            .join(c2s, JOIN.LEFT_OUTER)
            .join(s, JOIN.LEFT_OUTER)
            .where(
                (s.specobjid.is_null()) |
                (s.zwarning != 0 ) |
                (s.snmedian < self.parameters['spec_sn_thresh']) |
                (s.zerr >  self.parameters['spec_z_err_thresh'])
            )
            .distinct([t.gaia_sourceid])   # avoid duplicates - trust the GUA parent sample
        )

        return query


#-------  bhm_gua_dark   ------ #

class BhmGuaDarkCarton(BhmGuaBaseCarton):
    '''
        AND WHERE ( gua.g > 16.x AND gua.rp > 16.x)
    '''
    name = 'bhm_gua_dark'
    cadence = 'bhm_spiders_1x4'

#-------  bhm_gua_bright   ------ #

class BhmGuaBrightCarton(BhmGuaBaseCarton):
    '''
        AND WHERE ( gua.g < 18.x OR gua.rp < 18.x)
    '''
    name = 'bhm_gua_bright'
    cadence = 'bhm_boss_bright_3x1'




'''
Exporting from the temp table

\copy (SELECT * FROM sandbox.temp_bhm_gaia_unwise_agn_dark)  TO '/home/tdwelly/scratch/targetdb/bhm_gaia_unwise_agn_dark.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_gaia_unwise_agn_bright)  TO '/home/tdwelly/scratch/targetdb/bhm_gaia_unwise_agn_bright.csv' with csv header

for F in bhm_gaia_unwise_agn_*.csv; do   stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic; done
m
'''
