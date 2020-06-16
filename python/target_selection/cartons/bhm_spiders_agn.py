#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_spiders_agn.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

import peewee
from peewee import JOIN
from peewee import fn
import sdssdb


'''
# for testing do this

import peewee
import sdssdb
from sdssdb.peewee.sdss5db.catalogdb import database
database.set_profile('tunnel_operations')
from target_selection.cartons.bhm_spiders_agn import BhmSpidersAgnEfedsCarton
b = BhmSpidersAgnEfedsCarton(targeting_plan='0.1.0-beta.1')
q = b.build_query(version_id=13)
for r in q[:5]: print(r))

'''



from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             BHM_Spiders_AGN_Superset,
                                             Legacy_Survey_DR8,
                                             BHM_eFEDS_Veto,
                                             CatalogToLegacy_Survey_DR8,
#TODO                                             CatalogToSDSS_DR16_SpecObj,
#TODO                                             SDSS_DR16_SpecObj,
                                             )


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

x = BHM_Spiders_AGN_Superset.alias()
ls = Legacy_Survey_DR8.alias()
#q14 = SDSS_DR14_QSO.alias()
c = Catalog.alias()


for f in c._meta.fields:
    print (f)
for f in x._meta.fields:
    print (f)
for f in ls._meta.fields:
    print (f)
#for f in q14._meta.fields:
#    print (f)

'''
####

## This provides the following BHM cartons:

# BHM_SPIDERS_AGN_WIDE
# BHM_SPIDERS_AGN_DEEP
# BHM_SPIDERS_AGN_EFEDS
# EROSITA_POINTLIKE_BRIGHT_BOSS



class BhmSpidersBaseCarton(BaseCarton):
    ''' Parent class that provides the mask selections for any SPIDERs catalogue'''

    name = 'bhm_spiders'
    category = 'science'
    mapper = 'BHM'
    program = 'BHM-SPIDERS'
    tile = False

    # list of skymasks - move to the config file?
    skymasks = [
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



class BhmSpidersWideBaseCarton(BhmSpidersBaseCarton):
    ''' Parent class that provides the mask slections for any SPIDER-wide catalogue'''

    name = 'bhm_spiders_wide'
    category = 'science'
    mapper = 'BHM'
    program = 'BHM-SPIDERS'
    tile = False

    # list of skymasks
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


class BhmSpidersAgnWideLsCarton(BhmSpidersWideBaseCarton):

    '''
    spiders_agn_wide_ls:

    SELECT * from bhm_spiders_agn_superset AS x
    INNER JOIN legacy_survey_dr8 AS ls ON x.ls_id = ls.ls_id
    WHERE x.ero_version = "version_code_TBD"
    AND WHERE x.ero_det_like > X.X
    AND WHERE x.xmatch_metric > 0.x
    AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
    AND WHERE ls.fibertotflux_r < CCC.C
    '''

    name = 'bhm_spiders_agn_wide_ls'
    cadence = 'bhm_spiders_1x4'

    def build_query(self, version_id):
        c = Catalog.alias()
        x = BHM_Spiders_AGN_Superset.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()

        flux_r_max =  AB2nMgy(self.parameters['mag_r_min'])
        flux_r_min =  AB2nMgy(self.parameters['mag_r_max'])
        flux_z_min =  AB2nMgy(self.parameters['mag_z_max'])
        target_value = peewee.Value(self.parameters.get('value', 1.0)).alias('value')

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
                    ls.fiberflux_z.alias('lsfiberflux_z'),
            )
            .join(c2ls)
            .join(ls)
            .join(x)
            .where(c.version_id == version_id,
                   c2ls.version_id == version_id)
            .where(
                (x.ero_version == self.parameters['ero_version'] ),
                (ls.fibertotflux_r < flux_r_max),
                ((ls.fiberflux_r   > flux_r_min) |
                 (ls.fiberflux_z > flux_z_min) ),
                (x.ero_det_like > self.parameters['det_like_min']),
                (x.xmatch_metric > self.parameters['p_any_min']),
            )
        )

        print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")

        return query





#class BhmSpidersAgnEfedsCarton(BhmSpidersBaseCarton):
class BhmSpidersAgnEfedsCarton(BaseCarton):

    '''
    SELECT * from bhm_spiders_agn_superset AS x
    INNER JOIN legacy_survey_dr8 AS ls ON  x.ls_id = ls.ls_id
    LEFT JOIN bhm_efeds_veto AS v ON q3c_join(x.opt_ra,x.opt_dec,v.plug_ra,v.plug_dec,1.0)
    WHERE x.ero_version = "efeds_c940_V2T"
    AND WHERE x.ero_det_like > X.X
    AND WHERE (ls.fiberflux_r > D.D OR
               ls.fiberflux_z > E.E)                           # faint limits
    AND WHERE ls.fibertotflux_r < CCC.C                        # bright limit - use total flux
    AND WHERE (v.plug_ra = NULL OR
               v.sn_median_all < 1.x OR
               v.zwarning > 0 OR
               v.z_err > 2e-3 OR
               v.z_err <= 0.0)
    '''

    category = 'science'
    mapper = 'BHM'
    program = 'BHM-SPIDERS'
    tile = False


    name = 'bhm_spiders_agn_efeds'
    cadence = None # 'bhm_spiders_1x8'

    # config = {
    #     'ero_version':'efeds_c940_V2T',
    #     'mag_r_min': 17.0,
    #     'mag_r_max': 22.5,
    #     'mag_z_max': 21.5,
    #     'det_like_min': 6.0,
    #     'p_any_min': 0.1,
    #     'lr_min': 0.2,
    #     'veto_join_radius': 1.0,
    #     'veto_sn_thresh': 1.0000,
    #     'veto_z_err_thresh': 0.002,
    # }



    def build_query(self, version_id):

        c = Catalog.alias()
        x = BHM_Spiders_AGN_Superset.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()
        v = BHM_eFEDS_Veto.alias()
#TODO        c2s = CatalogToSDSS_dr16_SpecObj.alias()
#TODO        s = SDSS_dr16_SpecObj.alias()

        flux_r_max = AB2nMgy(self.parameters['mag_r_min'])
        flux_r_min = AB2nMgy(self.parameters['mag_r_max'])
        flux_z_min = AB2nMgy(self.parameters['mag_z_max'])

        target_value = peewee.Value(self.parameters.get('value', 1.0)).alias('value')
        match_radius_spectro = self.parameters['veto_join_radius']/3600.0

#        priority_val = Case(None, (
#            (s.specobjid.is_null(), 1510),
#            (s.specNumber.val == 2, 'two'),
#            (Number.val == 3, 'three')),
#                            'a lot')

        query = (
            c
            .select(c.catalogid,
                    (1510 + x.target_priority).alias('priority'),    ## catalog input is always == 1
                    target_value,
                    (22.5-2.5*fn.log10(fn.greatest(1e-10,ls.flux_g))).alias('magnitude_g'),
                    (22.5-2.5*fn.log10(fn.greatest(1e-10,ls.flux_r))).alias('magnitude_r'),
                    (22.5-2.5*fn.log10(fn.greatest(1e-10,ls.flux_z))).alias('magnitude_z'),
            )
            .join(c2ls)
            .join(ls)
            .join(x)
            .join(v, JOIN.LEFT_OUTER,
                  on=fn.q3c_join(c.ra,c.dec,
                                 v.plug_ra,v.plug_dec,
                                 match_radius_spectro))
#TODO            .switch(c)
#TODO            .join(c2s, JOIN.LEFT_OUTER)
#TODO            .join(s, JOIN.LEFT_OUTER)
            .where(c.version_id == version_id,
                   c2ls.version_id == version_id)
            .distinct([ls.ls_id])   # avoid duplicates - trust the ls_id
            .where(
                (x.ero_version == self.parameters['ero_version'] ),
                (
                    (v.plate.is_null()) |
                    (v.sn_median_all < self.parameters['veto_sn_thresh']) |
                    (v.zwarning > 0) |
                    (v.z_err >= self.parameters['veto_z_err_thresh']) |
                    (v.z_err <= 0.0)
                ),
                (ls.fibertotflux_r < flux_r_max),
                (
                    (ls.fiberflux_r >= flux_r_min) |
                    (ls.fiberflux_z >= flux_z_min)
                ),
                (x.ero_det_like > self.parameters['det_like_min']),
                (
                    (
                        (x.xmatch_method == 'XPS-ML/NWAY') &
                        (x.xmatch_metric >= self.parameters['p_any_min'])
                    ) |
                    (
                        (x.xmatch_method == 'XPS-LR') &
                        (x.xmatch_metric >= self.parameters['lr_min'])
                    )
                )
            )
        )

#        print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")

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
#waiting_for_psdr2#     def build_query(self, version_id):
#waiting_for_psdr2#
#waiting_for_psdr2#         c = Catalog.alias()
#waiting_for_psdr2#         x = BHM_Spiders_AGN_Superset.alias()
#waiting_for_psdr2#         ps = PanStarrsDr2.alias()
#waiting_for_psdr2#         c2ps = CatalogToPanStarrsDr2.alias()
#waiting_for_psdr2#
#waiting_for_psdr2#         flux_r_max =  AB2Jy(self.parameters['r_mag_min'])
#waiting_for_psdr2#         flux_r_min =  AB2Jy(self.parameters['r_mag_max'])
#waiting_for_psdr2#         flux_z_min =  AB2Jy(self.parameters['z_mag_max'])
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
#waiting_for_psdr2#                 (tab.ero_det_like > self.parameters['det_like_min']) &
#waiting_for_psdr2#                 (tab.xmatch_metric > self.parameters['p_any_min'])
#waiting_for_psdr2#             )
#waiting_for_psdr2#         )
#waiting_for_psdr2#
#waiting_for_psdr2#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")
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
#                  .where((tab.target_mag_r > self.parameters['r_mag_min']) &
#                         (tab.target_mag_r < self.parameters['r_mag_max']) &
#                         (tab.ero_det_like_0 > self.parameters['det_like_0_min'])))
#
#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")
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




def _test_xmatch_stuff():
    # check numbers of targets in a test patch

    version_id = 11
    search_radius_deg = 0.01
    ra0 = 135.0
    dec0 = 1.0

    c = Catalog.alias()
#    x = BHM_Spiders_AGN_Superset.alias()
    ls = Legacy_Survey_DR8.alias()
    c2ls = CatalogToLegacy_Survey_DR8.alias()


    query = (
        c
        .select(c.catalogid,
                c.ra,
                c.dec,
                c.lead,
                c.version,
                ls.ls_id,
                ls.ra,
                ls.dec,
        )
        .join(c2ls)
        .join(ls)
        .where(c.version_id == version_id,
               c2ls.version_id == version_id)
        .where(
            peewee.fn.q3c_radial_query(c.ra,c.dec,
                                       ra0, dec0,
                                       search_radius_deg)
        )
    )


    query.select().limit(1000).count()
