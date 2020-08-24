#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_spiders_agn.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

# flake8: noqa
# isort: skip_file

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
                                             CatalogToLegacy_Survey_DR8,
                                             Legacy_Survey_DR8,
                                             BHM_Spiders_AGN_Superset,
                                             BHM_eFEDS_Veto,
#TODO                                             CatalogToSDSS_DR16_SpecObj,
                                             SDSS_DR16_SpecObj,
                                             )


#waiting_for_psdr2# ,  PanStarrsDr2)
#waiting_for_psdr2# ,  CatalogToPanStarrsDr2)

from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import psfmag_minus_fiber2mag, AB2nMgy


## This provides the following BHM cartons:
#    bhm_spiders_agn-efeds
# and will eventually provide:
#    bhm_spiders_agn-wide
#    bhm_spiders_agn-deep
# maybe:
#    erosita-pointlike-bright-boss





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
    program = 'bhm_spiders'
    tile = False


    name = 'bhm_spiders_agn-efeds'
    cadence = 'bhm_spiders_1x8'



    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        x = BHM_Spiders_AGN_Superset.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()
        v = BHM_eFEDS_Veto.alias()
#TODO        c2s = CatalogToSDSS_DR16_SpecObj.alias()
        s = SDSS_DR16_SpecObj.alias()

        fiberflux_r_max = AB2nMgy(self.parameters['fibermag_r_min'])
        fiberflux_r_min = AB2nMgy(self.parameters['fibermag_r_max'])
        fiberflux_z_min = AB2nMgy(self.parameters['fibermag_z_max'])

        flux30 = AB2nMgy(30.00)

        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float').alias('value')
        match_radius_spectro = self.parameters['spec_join_radius']/3600.0

        p_f = self.parameters['priority_floor']
        priority = peewee.Case(None,
                               (
                                   ((x.ero_det_like < self.parameters['det_like_for_priority']) & (s.specobjid.is_null(True)), p_f+6),
                                   ((x.ero_det_like < self.parameters['det_like_for_priority']) & (s.specobjid.is_null(False)), p_f+7),
                                   ((x.xmatch_flags == 1 ) & (s.specobjid.is_null(True)), p_f+0),
                                   ((x.xmatch_flags == 0 ) & (s.specobjid.is_null(True)), p_f+1),
                                   ((x.xmatch_flags > 1  ) & (s.specobjid.is_null(True)), p_f+2),
                                   ((x.xmatch_flags == 1 ) & (s.specobjid.is_null(False)), p_f+3),
                                   ((x.xmatch_flags == 0 ) & (s.specobjid.is_null(False)), p_f+4),
                                   ((x.xmatch_flags > 1  ) & (s.specobjid.is_null(False)), p_f+5),
                               ),
                               p_f+9) ## should never get here


        # Notes on convertion from ls_fibermag to sdss_fiber2mag:
        # https://wiki.mpe.mpg.de/eRosita/EroAGN_eFEDS/SDSSIVSpecialPlates#Estimating_SDSS_fiber2mag_.2A_from_legacysurvey_photometry

        # Notes on converting from sdss_fiber2mag to sdss_psfmag
        # https://wiki.sdss.org/display/OPS/Contents+of+targetdb.magnitude#Contentsoftargetdb.magnitude-WhatmagnitudestoputintotheplPlugMapfilesforBOSSplatetargets?

        # A flux ratio of 0.6 (roughly what is seen in all three ls bands) is a magnitude difference of
        # fiber2mag(SDSS)-fibermag(LS) = 0.55mags
        flux_ratio = {'g' : 0.60, 'r' : 0.60, 'i' : 0.60, 'z' : 0.60 }
        # Then, also add the correction from sdss_fiber2mag to sdss_psfmag: psfmag_minus_fiber2mag('filter')

        # legacysurvey fibermags - derived from fiberfluxes - with limits to avoid divide by zero errors
        magnitude_g = (psfmag_minus_fiber2mag('g') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_g*flux_ratio['g']))).cast('float')
        magnitude_r = (psfmag_minus_fiber2mag('r') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_r*flux_ratio['r']))).cast('float')
        magnitude_z = (psfmag_minus_fiber2mag('z') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_z*flux_ratio['z']))).cast('float')
        # the simplest possible interpolation
        # TODO - we can do this better
        magnitude_i = (psfmag_minus_fiber2mag('i') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,
                                                     0.5*(ls.fiberflux_r+
                                                          ls.fiberflux_z)*flux_ratio['i']))).cast('float')


        query = (
            c
            .select(c.catalogid,
                    priority.alias('priority'),
                    value,
                    magnitude_g.alias("g"),
                    magnitude_r.alias("r"),
                    magnitude_i.alias("i"),
                    magnitude_z.alias("z"),
            )
            .join(c2ls)
            .join(ls)
            .join(x)
            .join(v, JOIN.LEFT_OUTER,
                  on=(fn.q3c_join(c.ra,c.dec,
                                 v.plug_ra,v.plug_dec,
                                 match_radius_spectro) &
                      (v.sn_median_all >= self.parameters['spec_sn_thresh']) &
                      (v.zwarning == 0) &
                      (v.z_err <= self.parameters['spec_z_err_thresh']) &
                      (v.z_err > 0.0)
                      )
                  )
            .join(s, JOIN.LEFT_OUTER,
                  on=(fn.q3c_join(c.ra,c.dec,
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
                   c2ls.version_id == version_id,
                   c2ls.best == True)
            .where(
                (x.ero_version == self.parameters['ero_version'] ),
                (v.pk.is_null()),
                (ls.fibertotflux_r < fiberflux_r_max),
                (
                    (ls.fiberflux_r >= fiberflux_r_min) |
                    (ls.fiberflux_z >= fiberflux_z_min)
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
            .distinct([ls.ls_id])   # avoid duplicates - trust the ls_id
        )

        return query










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




#deferred#
#deferred#
#deferred#
#deferred#
#deferred# class BhmSpidersBaseCarton(BaseCarton):
#deferred#     ''' Parent class that provides the mask selections for any SPIDERs catalogue'''
#deferred#
#deferred#     name = 'bhm_spiders'
#deferred#     category = 'science'
#deferred#     mapper = 'BHM'
#deferred#     program = 'BHM-SPIDERS'
#deferred#     tile = False
#deferred#
#deferred#     # list of skymasks - move to the config file?
#deferred#     skymasks = [
#deferred#     ]
#deferred#
#deferred#
#deferred#     def post_process(self, model, **kwargs):
#deferred#         # this is where we select on mask location
#deferred#         # get the ids and coords of all of the objects in the temp table
#deferred#         cat_id = np.array(model.catalog_id[:])
#deferred#         ra = np.array(model.ra[:])
#deferred#         dec = np.array(model.dec[:])
#deferred#
#deferred#         flags = np.ones(len(cat_id), np.bool)
#deferred#
#deferred#         for sm in self.skymasks:
#deferred#             sm.apply(lon=ra, lat=dec, flags=flags)
#deferred#
#deferred#         # not sure what to return at this point - a list of tuples?
#deferred#         result = [(i,flag) for i,flag in zip(cat_id,flags)]
#deferred#
#deferred#         return result
#deferred#
#deferred#
#deferred#
#deferred# class BhmSpidersWideBaseCarton(BhmSpidersBaseCarton):
#deferred#     ''' Parent class that provides the mask slections for any SPIDER-wide catalogue'''
#deferred#
#deferred#     name = 'bhm_spiders_wide'
#deferred#     category = 'science'
#deferred#     mapper = 'BHM'
#deferred#     program = 'BHM-SPIDERS'
#deferred#     tile = False
#deferred#
#deferred#     # list of skymasks
#deferred#     skymasks = [
#deferred#         SkyMask(filename=pkg_resources.resource_filename(
#deferred#             __name__,
#deferred#             'masks/eROSITA-DE_exgal_lsdr8_or_psdr2_proc.ply'),
#deferred#                 name="spiders_wide",
#deferred#                 masktype="mangle",
#deferred#                 sense="include",
#deferred#         ),
#deferred#         SkyMask(filename=pkg_resources.resource_filename(
#deferred#             __name__,
#deferred#             'masks/rsFields-annotated-lco-deep_proc.ply'),
#deferred#                 name="spiders_deep",
#deferred#                 masktype="mangle",
#deferred#                 sense="exclude",
#deferred#         ),
#deferred#     ]
#deferred#
#deferred#
#deferred# class BhmSpidersAgnWideLsCarton(BhmSpidersWideBaseCarton):
#deferred#
#deferred#     '''
#deferred#     spiders_agn_wide_ls:
#deferred#
#deferred#     SELECT * from bhm_spiders_agn_superset AS x
#deferred#     INNER JOIN legacy_survey_dr8 AS ls ON x.ls_id = ls.ls_id
#deferred#     WHERE x.ero_version = "version_code_TBD"
#deferred#     AND WHERE x.ero_det_like > X.X
#deferred#     AND WHERE x.xmatch_metric > 0.x
#deferred#     AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
#deferred#     AND WHERE ls.fibertotflux_r < CCC.C
#deferred#     '''
#deferred#
#deferred#     name = 'bhm_spiders_agn_wide_ls'
#deferred#     cadence = 'bhm_spiders_1x4'
#deferred#
#deferred#     def build_query(self, version_id, query_region=None):
#deferred#         c = Catalog.alias()
#deferred#         x = BHM_Spiders_AGN_Superset.alias()
#deferred#         ls = Legacy_Survey_DR8.alias()
#deferred#         c2ls = CatalogToLegacy_Survey_DR8.alias()
#deferred#
#deferred#         flux_r_max =  AB2nMgy(self.parameters['mag_r_min'])
#deferred#         flux_r_min =  AB2nMgy(self.parameters['mag_r_max'])
#deferred#         flux_z_min =  AB2nMgy(self.parameters['mag_z_max'])
#deferred#         target_value = peewee.Value(self.parameters.get('value', 1.0)).alias('value')
#deferred#
#deferred#         query = (
#deferred#             c
#deferred#             .select(c.catalogid,
#deferred#                     c.ra,
#deferred#                     c.dec,
#deferred#                     c.pmra,
#deferred#                     c.pmdec,
#deferred#                     x.target_priority.alias('priority'),
#deferred#                     ls.fiberflux_g.alias('lsfiberflux_g'),
#deferred#                     ls.fiberflux_r.alias('lsfiberflux_r'),
#deferred#                     ls.fiberflux_z.alias('lsfiberflux_z'),
#deferred#             )
#deferred#             .join(c2ls)
#deferred#             .join(ls)
#deferred#             .join(x)
#deferred#             .where(c.version_id == version_id,
#deferred#                    c2ls.version_id == version_id)
#deferred#             .where(
#deferred#                 (x.ero_version == self.parameters['ero_version'] ),
#deferred#                 (ls.fibertotflux_r < flux_r_max),
#deferred#                 ((ls.fiberflux_r   > flux_r_min) |
#deferred#                  (ls.fiberflux_z > flux_z_min) ),
#deferred#                 (x.ero_det_like > self.parameters['det_like_min']),
#deferred#                 (x.xmatch_metric > self.parameters['p_any_min']),
#deferred#             )
#deferred#         )
#deferred#
#deferred#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")
#deferred#
#deferred#         return query
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred# #waiting_for_psdr2# class BhmSpidersAgnWidePsCarton(BhmSpidersWideBaseCarton):
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#     '''
#deferred# #waiting_for_psdr2#     spiders_agn_wide_ps:
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#     SELECT * from bhm_spiders_agn_superset AS x
#deferred# #waiting_for_psdr2#     INNER JOIN panstarrs_dr2 AS ps ON x.ps1_dr2_objid = ls.ls_id
#deferred# #waiting_for_psdr2#     WHERE x.ero_version = "version_code_TBD"
#deferred# #waiting_for_psdr2#     AND WHERE x.ero_det_like > X.X
#deferred# #waiting_for_psdr2#     AND WHERE x.xmatch_metric > 0.x
#deferred# #waiting_for_psdr2#     AND WHERE (ls.fiberflux_r > A.A OR ls.fiberflux_z > B.B)
#deferred# #waiting_for_psdr2#     AND WHERE ls.fiberflux_r < CCC.C
#deferred# #waiting_for_psdr2#     '''
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#     name = 'bhm_spiders_agn_wide_ls'
#deferred# #waiting_for_psdr2#     cadence = 'bhm_spiders_1x4'
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#     def build_query(self, version_id, query_region=None):
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#         c = Catalog.alias()
#deferred# #waiting_for_psdr2#         x = BHM_Spiders_AGN_Superset.alias()
#deferred# #waiting_for_psdr2#         ps = PanStarrsDr2.alias()
#deferred# #waiting_for_psdr2#         c2ps = CatalogToPanStarrsDr2.alias()
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#         flux_r_max =  AB2Jy(self.parameters['r_mag_min'])
#deferred# #waiting_for_psdr2#         flux_r_min =  AB2Jy(self.parameters['r_mag_max'])
#deferred# #waiting_for_psdr2#         flux_z_min =  AB2Jy(self.parameters['z_mag_max'])
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#         query = (
#deferred# #waiting_for_psdr2#             c
#deferred# #waiting_for_psdr2#             .select(c.catalogid,
#deferred# #waiting_for_psdr2#                     c.ra,
#deferred# #waiting_for_psdr2#                     c.dec,
#deferred# #waiting_for_psdr2#                     c.pmra,
#deferred# #waiting_for_psdr2#                     c.pmdec,
#deferred# #waiting_for_psdr2#                     x.target_priority.alias('priority'),
#deferred# #waiting_for_psdr2#                     ps.g_stk_aper_flux.alias('psaperflux_g'),
#deferred# #waiting_for_psdr2#                     ps.r_stk_aper_flux.alias('psaperflux_r'),
#deferred# #waiting_for_psdr2#                     ps.i_stk_aper_flux.alias('psaperflux_i'),
#deferred# #waiting_for_psdr2#                     ps.z_stk_aper_flux.alias('psaperflux_z'))
#deferred# #waiting_for_psdr2#             .join(ps)
#deferred# #waiting_for_psdr2#             .join(c2ps)
#deferred# #waiting_for_psdr2#             .where(
#deferred# #waiting_for_psdr2#                 (ps.r_stk_aper_flux < flux_r_max) &
#deferred# #waiting_for_psdr2#                 ((ps.r_stk_aper_flux > flux_r_min) | (ps.z_stk_aper_flux > flux_z_min) ) &
#deferred# #waiting_for_psdr2#                 (tab.ero_det_like > self.parameters['det_like_min']) &
#deferred# #waiting_for_psdr2#                 (tab.xmatch_metric > self.parameters['p_any_min'])
#deferred# #waiting_for_psdr2#             )
#deferred# #waiting_for_psdr2#         )
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")
#deferred# #waiting_for_psdr2#
#deferred# #waiting_for_psdr2#         return query
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred#
#deferred# #
#deferred# # ##############################################################
#deferred# # ##############################################################
#deferred# # ##############################################################
#deferred# # ## MOCK TARGETS ##############################################
#deferred# # class BhmSpidersAgnWideMockCarton(BaseCarton):
#deferred# #
#deferred# #     name = 'bhm_spiders_agn_wide_mock'
#deferred# #     category = 'science'
#deferred# #     survey = 'BHM'
#deferred# #     cadence = 'bhm_spiders_1x4'
#deferred# #     tile = False
#deferred# #
#deferred# #     # list of masks - remove to config file?
#deferred# #     masks = [
#deferred# #         {"name": "wide",
#deferred# #          "type": "mangle",
#deferred# #          "polarity": "include",
#deferred# #          "filename": "eROSITA-DE_exgal_lsdr8_or_psdr2_proc.ply",
#deferred# #         },
#deferred# #         {"name": "deep",
#deferred# #          "type": "mangle",
#deferred# #          "polarity": "exclude",
#deferred# #          "filename" : "rsFields-annotated-lco-deep_proc.ply",
#deferred# #         },
#deferred# #     ]
#deferred# #
#deferred# #     def build_query(self, version_id, query_region=None):
#deferred# #         '''
#deferred# #         Pure database level query - generates a super-set of potential targets
#deferred# #         '''
#deferred# #         # get the table name from the config - maybe replace this with a list of options
#deferred# #         exec(f'tab = catalogdb.{params["catalogdb_table"]}')
#deferred# #         assert tab is not None, 'Failed to locate catalogdb_table'
#deferred# #
#deferred# #         query = (tab.select(tab.gaia_dr2_source_id.alias('catalog_id'),
#deferred# #                             tab.target_ra.alias('ra'),
#deferred# #                             tab.target_dec.alias('dec'),
#deferred# #                             tab.target_pmra.alias('pmra'),
#deferred# #                             tab.target_pmdec.alias('pmdec'),
#deferred# #                             tab.target_epoch.alias('epoch'),
#deferred# #                             tab.target_mag_r.alias('magnitude_g'),
#deferred# #                             tab.target_mag_r.alias('magnitude_r'),
#deferred# #                             tab.target_mag_r.alias('magnitude_z'))
#deferred# #                  .where((tab.target_mag_r > self.parameters['r_mag_min']) &
#deferred# #                         (tab.target_mag_r < self.parameters['r_mag_max']) &
#deferred# #                         (tab.ero_det_like_0 > self.parameters['det_like_0_min'])))
#deferred# #
#deferred# #         print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.parameters['req_ntargets']})")
#deferred# #
#deferred# #         return query
#deferred# #
#deferred# #
#deferred# #     def post_process(self, model, **kwargs):
#deferred# #         # this is where we select on mask location
#deferred# #         # get the coords of all of the objects oin the temp table
#deferred# #         ra = np.array(model.ra[:])
#deferred# #         dec = np.array(model.dec[:])
#deferred# #
#deferred# #         return True
#deferred# #
#deferred# # ## MOCK TARGETS ##############################################
#deferred# # ##############################################################
#deferred# # ##############################################################
#deferred# # ##############################################################
#deferred#
#deferred#
#deferred#
#deferred#
#deferred# def _test_xmatch_stuff():
#deferred#     # check numbers of targets in a test patch
#deferred#
#deferred#     version_id = 11
#deferred#     search_radius_deg = 0.01
#deferred#     ra0 = 135.0
#deferred#     dec0 = 1.0
#deferred#
#deferred#     c = Catalog.alias()
#deferred# #    x = BHM_Spiders_AGN_Superset.alias()
#deferred#     ls = Legacy_Survey_DR8.alias()
#deferred#     c2ls = CatalogToLegacy_Survey_DR8.alias()
#deferred#
#deferred#
#deferred#     query = (
#deferred#         c
#deferred#         .select(c.catalogid,
#deferred#                 c.ra,
#deferred#                 c.dec,
#deferred#                 c.lead,
#deferred#                 c.version,
#deferred#                 ls.ls_id,
#deferred#                 ls.ra,
#deferred#                 ls.dec,
#deferred#         )
#deferred#         .join(c2ls)
#deferred#         .join(ls)
#deferred#         .where(c.version_id == version_id,
#deferred#                c2ls.version_id == version_id)
#deferred#         .where(
#deferred#             peewee.fn.q3c_radial_query(c.ra,c.dec,
#deferred#                                        ra0, dec0,
#deferred#                                        search_radius_deg)
#deferred#         )
#deferred#     )
#deferred#
#deferred#
#deferred#     query.select().limit(1000).count()
#deferred#


'''
Exporting from the temp table

\copy (SELECT * FROM sandbox.temp_bhm_spiders_agn_efeds)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_agn_efeds.csv' with csv header
stilts tpipe in=~/scratch/targetdb/bhm_spiders_agn_efeds.csv out=~/scratch/targetdb/bhm_spiders_agn_efeds.fits ifmt=csv ofmt=fits-basic


'''
