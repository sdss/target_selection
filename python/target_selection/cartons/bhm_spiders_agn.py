#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_spiders_agn.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

import peewee
import sdssdb

from sdssdb.peewee.sdss5db import catalogdb

from target_selection.cartons.base import BaseCarton
from target_selection.cartons.skymask import SkyMask
import pkg_resources


#### some useful snippets:
'''
# example to get the listing of fields from a PeeWee model
print(catalogdb.ErositaAGNMock._meta.fields)

sp = catalogdb.BHM_Spiders_AGN_Superset.alias()
ls = catalogdb.Legacy_Survey_DR8.alias()
q14 = catalogdb.SDSS_DR14_QSO.alias()


for f in sp._meta.fields:
    print (f)
for f in ls._meta.fields:
    print (f)
for f in q14._meta.fields:
    print (f)

'''
####

## This provides the following BHM cartons:

# BHM_SPIDERS_AGN_WIDE
# BHM_SPIDERS_AGN_DEEP
# BHM_SPIDERS_AGN_EFEDS
# EROSITA_POINTLIKE_BRIGHT_BOSS

class BhmSpidersAgnWideCarton(BaseCarton):

    name = 'bhm_spiders_agn_wide'
    category = 'science'
    survey = 'BHM'
    cadence = 'bhm_spiders_1x4'
    tile = False

    # list of masks - move to the config file?
    masks = [
        SkyMask(filename=pkg_resources.resource_filename(
            __name__,
            'masks/eROSITA-DE_exgal_lsdr8_or_psdr2_proc.ply'),
                name="spiders_wide",
                masktype="mangle",
                sense="include",
        ),
        SkyMask(filename=pkg_resources.resource_filename(
            __name__,
            'masks/rsFields-annotated-lco-deep_proc.ply'),
                name="spiders_deep",
                masktype="mangle",
                sense="exclude",
        ),
    ]

    def build_query(self):
        '''
        Pure database level query - generates a super-set of potential targets
        '''
        # get the table name from the config - maybe replace this with a list of options
        #exec(f'tab = catalogdb.{params["catalogdb_table"]}')
        #assert tab is not None, 'Failed to locate catalogdb table'

        sp = catalogdb.BHM_Spiders_AGN_Superset.alias()
        ls = catalogdb.Legacy_Survey_DR8.alias()
        #ps = catalogdb.PanStarrsDr2.alias()

        query = (sp.select(sp.ls_id.alias('catalog_id'),
                            sp.opt_ra.alias('ra'),
                            sp.opt_dec.alias('dec'),
                            sp.opt_pmra.alias('pmra'),
                            sp.opt_pmdec.alias('pmdec'),
                            sp.opt_epoch.alias('epoch'),
                            ls.fiberflux_g.alias('fiberflux_g'),
                            ls.fiberflux_r.alias('fiberflux_r'),
                            ls.fiberflux_z.alias('fiberflux_z'))
                 .join(ls)
                 .where((tab.target_mag_r > self.config['r_mag_min']) &
                        (tab.target_mag_r < self.config['r_mag_max']) &
                        (tab.ero_det_like > self.config['det_like_min'])))

        print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.config['req_ntargets']})")

        return query



    def post_process(self, model, **kwargs):
        # this is where we select on mask location
        # get the id and coords of all of the objects in the temp table
        cat_id = np.array(model.catalog_id[:])
        ra = np.array(model.ra[:])
        dec = np.array(model.dec[:])


        return True















#
# ##############################################################
# ##############################################################
# ##############################################################
# ## MOCK TARGETS ##############################################
# class BhmSpidersAgnWideMockCarton(BaseCarton):
#
#     name = 'bhm_spiders_agn_wide_mock'
#     category = 'science'
#     survey = 'BHM'
#     cadence = 'bhm_spiders_1x4'
#     tile = False
#
#     # list of masks - remove to config file?
#     masks = [
#         {"name": "wide",
#          "type": "mangle",
#          "polarity": "include",
#          "filename": "eROSITA-DE_exgal_lsdr8_or_psdr2_proc.ply",
#         },
#         {"name": "deep",
#          "type": "mangle",
#          "polarity": "exclude",
#          "filename" : "rsFields-annotated-lco-deep_proc.ply",
#         },
#     ]
#
#     def build_query(self):
#         '''
#         Pure database level query - generates a super-set of potential targets
#         '''
#         # get the table name from the config - maybe replace this with a list of options
#         exec(f'tab = catalogdb.{params["catalogdb_table"]}')
#         assert tab is not None, 'Failed to locate catalogdb_table'
#
#         query = (tab.select(tab.gaia_dr2_source_id.alias('catalog_id'),
#                             tab.target_ra.alias('ra'),
#                             tab.target_dec.alias('dec'),
#                             tab.target_pmra.alias('pmra'),
#                             tab.target_pmdec.alias('pmdec'),
#                             tab.target_epoch.alias('epoch'),
#                             tab.target_mag_r.alias('magnitude_g'),
#                             tab.target_mag_r.alias('magnitude_r'),
#                             tab.target_mag_r.alias('magnitude_z'))
#                  .where((tab.target_mag_r > self.config['r_mag_min']) &
#                         (tab.target_mag_r < self.config['r_mag_max']) &
#                         (tab.ero_det_like_0 > self.config['det_like_0_min'])))
#
#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.config['req_ntargets']})")
#
#         return query
#
#
#     def post_process(self, model, **kwargs):
#         # this is where we select on mask location
#         # get the coords of all of the objects oin the temp table
#         ra = np.array(model.ra[:])
#         dec = np.array(model.dec[:])
#
#         return True
#
# ## MOCK TARGETS ##############################################
# ##############################################################
# ##############################################################
# ##############################################################
