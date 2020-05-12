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

#from sdssdb.peewee.sdss5db import catalogdb
from sdssdb.peewee.sdss5db.catalogdb import (Catalog, BHM_Spiders_AGN_Superset, Legacy_Survey_DR8)
#, CatalogToLegacy_Survey_DR8)
#waiting_for_psdr2# ,  PanStarrsDr2)
#waiting_for_psdr2# ,  CatalogToPanStarrsDr2)

from target_selection.cartons.base import BaseCarton
from target_selection.cartons.skymask import SkyMask
from target_selection.mag_flux import *

import pkg_resources


#### some useful snippets:
'''
# example to get the listing of fields from a PeeWee model
print(catalogdb.ErositaAGNMock._meta.fields)

sp = BHM_Spiders_AGN_Superset.alias()
ls = Legacy_Survey_DR8.alias()
q14 = SDSS_DR14_QSO.alias()


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



class BhmSpidersWideBaseCarton(BaseCarton):
    ''' Parent class that provides the mask slections for any SPIDER-wide catalogue'''

    name = 'bhm_spiders_wide'
    category = 'science'
    survey = 'BHM'
    tile = False

    # list of skymasks - move to the config file?
    skymasks = [
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


    def post_process(self, model, **kwargs):
        # this is where we select on mask location
        # get the ids and coords of all of the objects in the temp table
        cat_id = np.array(model.catalog_id[:])
        ra = np.array(model.ra[:])
        dec = np.array(model.dec[:])

        flags = np.ones(len(cat_id), np.bool)

        for sm in self.skymasks:
            sm.apply(lon=ra, lat=dec, flags=flags)

        # not sure what to return at this point - a list of tuples?
        result = [(i,flag) for i,flag in zip(cat_id,flags)]

        return result




class BhmSpidersAgnWideLsCarton(BhmSpidersWideBaseCarton):

    '''
    spiders_agn_wide_ls:

    SELECT * from bhm_spiders_agn_superset AS x
    INNER JOIN legacy_survey_dr8 AS ls ON x.ls_id = ls.ls_id
    WHERE x.ero_version = "version_code_TBD"
    AND WHERE x.ero_det_like > X.X
    AND WHERE x.xmatch_metric > 0.x
    AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
    AND WHERE ls.fiberflux_r < CCC.C
    '''

    name = 'bhm_spiders_agn_wide_ls'
    cadence = 'bhm_spiders_1x4'

    def build_query(self):
        # get the table name from the config - maybe replace this with a list of options
        #exec(f'tab = catalogdb.{params["catalogdb_table"]}')
        #assert tab is not None, 'Failed to locate catalogdb table'

        c = Catalog.alias()
        x = BHM_Spiders_AGN_Superset.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()

        flux_r_max =  AB2nMgy(self.config['r_mag_min'])
        flux_r_min =  AB2nMgy(self.config['r_mag_max'])
        flux_z_min =  AB2nMgy(self.config['z_mag_max'])

        query = (
            c
            .select(c.catalogid,
                    c.ra,
                    c.dec,
                    c.pmra,
                    c.pmdec,
                    x.target_priority.alias('priority'),
                    ls.fiberflux_g.alias('lsfiberflux_g'),
                    ls.fiberflux_r.alias('lsfiberflux_r'),
                    ls.fiberflux_z.alias('lsfiberflux_z'))
            .join(c2ls)
            .join(ls)
            .join(x)
            .where(
                (ls.fibertotflux_r < flux_r_max) &
                ((ls.fiberflux_r   > flux_r_min) | (ls.fiberflux_z > flux_z_min) ) &
                (x.ero_det_like > self.config['det_like_min']) &
                (x.xmatch_metric > self.config['p_any_min'])
            )
        )

        print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.config['req_ntargets']})")

        return query




#waiting_for_psdr2# class BhmSpidersAgnWidePsCarton(BhmSpidersWideBaseCarton):
#waiting_for_psdr2#
#waiting_for_psdr2#     '''
#waiting_for_psdr2#     spiders_agn_wide_ps:
#waiting_for_psdr2#
#waiting_for_psdr2#     SELECT * from bhm_spiders_agn_superset AS x
#waiting_for_psdr2#     INNER JOIN panstarrs_dr2 AS ps ON x.ps1_dr2_objid = ls.ls_id
#waiting_for_psdr2#     WHERE x.ero_version = "version_code_TBD"
#waiting_for_psdr2#     AND WHERE x.ero_det_like > X.X
#waiting_for_psdr2#     AND WHERE x.xmatch_metric > 0.x
#waiting_for_psdr2#     AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
#waiting_for_psdr2#     AND WHERE ls.fiberflux_r < CCC.C
#waiting_for_psdr2#     '''
#waiting_for_psdr2#
#waiting_for_psdr2#     name = 'bhm_spiders_agn_wide_ls'
#waiting_for_psdr2#     cadence = 'bhm_spiders_1x4'
#waiting_for_psdr2#
#waiting_for_psdr2#     def build_query(self):
#waiting_for_psdr2#
#waiting_for_psdr2#         c = Catalog.alias()
#waiting_for_psdr2#         x = BHM_Spiders_AGN_Superset.alias()
#waiting_for_psdr2#         ps = PanStarrsDr2.alias()
#waiting_for_psdr2#         c2ps = CatalogToPanStarrsDr2.alias()
#waiting_for_psdr2#
#waiting_for_psdr2#         flux_r_max =  AB2Jy(self.config['r_mag_min'])
#waiting_for_psdr2#         flux_r_min =  AB2Jy(self.config['r_mag_max'])
#waiting_for_psdr2#         flux_z_min =  AB2Jy(self.config['z_mag_max'])
#waiting_for_psdr2#
#waiting_for_psdr2#         query = (
#waiting_for_psdr2#             c
#waiting_for_psdr2#             .select(c.catalogid,
#waiting_for_psdr2#                     c.ra,
#waiting_for_psdr2#                     c.dec,
#waiting_for_psdr2#                     c.pmra,
#waiting_for_psdr2#                     c.pmdec,
#waiting_for_psdr2#                     x.target_priority.alias('priority'),
#waiting_for_psdr2#                     ps.g_stk_aper_flux.alias('psaperflux_g'),
#waiting_for_psdr2#                     ps.r_stk_aper_flux.alias('psaperflux_r'),
#waiting_for_psdr2#                     ps.i_stk_aper_flux.alias('psaperflux_i'),
#waiting_for_psdr2#                     ps.z_stk_aper_flux.alias('psaperflux_z'))
#waiting_for_psdr2#             .join(ps)
#waiting_for_psdr2#             .join(c2ps)
#waiting_for_psdr2#             .where(
#waiting_for_psdr2#                 (ps.r_stk_aper_flux < flux_r_max) &
#waiting_for_psdr2#                 ((ps.r_stk_aper_flux > flux_r_min) | (ps.z_stk_aper_flux > flux_z_min) ) &
#waiting_for_psdr2#                 (tab.ero_det_like > self.config['det_like_min']) &
#waiting_for_psdr2#                 (tab.xmatch_metric > self.config['p_any_min'])
#waiting_for_psdr2#             )
#waiting_for_psdr2#         )
#waiting_for_psdr2#
#waiting_for_psdr2#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.config['req_ntargets']})")
#waiting_for_psdr2#
#waiting_for_psdr2#         return query















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
