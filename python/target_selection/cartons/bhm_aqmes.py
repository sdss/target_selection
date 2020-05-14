#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_spiders_agn.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

import peewee
#from peewee import JOIN
import sdssdb


from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             SDSS_DR14_QSO,)
#                                             SDSS_DR16_QSO )

from sdssdb.peewee.sdss5db.catalogdb import (CatalogToSDSS_DR14_QSO,)
#                                             CatalogToSDSS_DR16_QSO,)




from target_selection.cartons.base import BaseCarton
from target_selection.cartons.skymask import SkyMask
from target_selection.mag_flux import *

import pkg_resources



## This provides the following BHM cartons:

# BHM_AQMES_MED
# BHM_AQMES_MED_FAINT
# BHM_AQMES_WIDE2
# BHM_AQMES_WIDE3
# BHM_AQMES_WIDE2_FAINT
# BHM_AQMES_WIDE3_FAIINT
# BHM_AQMES_BONUS_DARK
# BHM_AQMES_BONUS_BRIGHT




class BhmAqmesBaseCarton(BaseCarton):
    ''' Parent class that provides the mask selections for all AQMES cartons'''

    name = 'bhm_aqmes'
    category = 'science'
    survey = 'BHM'
    tile = False

    # list of skymasks
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




class BhmAqmesMedCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr14_qso WHERE  psfmag_i BETWEEN 16.x AND 19.1
    '''

    name = 'bhm_aqmes_med'
    category = 'science'
    survey = 'BHM'
    tile = False
    cadence = 'bhm_aqmes_med_12x4'

    # list of skymasks
    skymasks = [
        SkyMask(filename=pkg_resources.resource_filename(
            __name__,
            'masks/candidate_target_fields_bhm_aqmes_med_v0.2.1.fits'),
                name="aqmes_med",
                masktype="circles",
                sense="include",
        ),
    ]

    def build_query(self):

        c = Catalog.alias()
        q = SDSS_DR14_QSO.alias()
        c2q = CatalogToSDSS_QSO_DR14.alias()

        query = (
            c
            .select(c.catalogid,
                    c.ra,
                    c.dec,
                    q.psfmag.alias('sdss_psfmag'),
            )
            .join(c2q)
            .join(q)
            .where(c.version_id == version_id,
                   c2q.version_id == version_id)
            .where(
                (q.psfmag[3] >= self.config['mag_i_min'] ),
                (q.psfmag[3] <= self.config['mag_i_max'] )
            )
        )

        print(f"This query will return nrows={query.count()}  (c.f. req_ntargets={self.config['req_ntargets']})")

        return query
