#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import peewee
import sdssdb
#from astropy.io import fits
#import pkg_resources

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             BHM_CSC,
                                             CatalogToBHM_CSC)


from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ChandraSourceCatalogue(CSC)-inprogress

## This module provides the following BHM cartons:
# bhm_csc_boss_dark
# bhm_csc_boss_bright
# bhm_csc_apogee



class BhmCscBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for all CSC cartons
    To be sub-classed, not typically to be called directly.
    '''

    name = 'bhm_csc_base'
    base_name = 'bhm_csc_base'
    category = 'science'
    mapper = 'BHM'
    program = 'BHM-CSC'
    tile = False
    priority = None
    alias_c = None
    alias_t = None
    instrument_code = None

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2t = CatalogToBHM_CSC.alias()
        t = BHM_CSC.alias()
        self.alias_c = c
        self.alias_t = t

        # set the Carton priority+values here - read from yaml
        target_priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')
        target_value = peewee.Value(float(self.parameters.get('value', 1.0))).alias('value')

        query = (
            c
            .select(c.catalogid,
                    target_priority,
                    target_value,
                    t.mag_g.alias('magnitude_g'),
                    t.mag_r.alias('magnitude_r'),
                    t.mag_i.alias('magnitude_i'),
                    t.mag_z.alias('magnitude_z'),
                    t.mag_h.alias('magnitude_h'),
            )
            .join(c2t)
            .join(t)
            .where(c.version_id == version_id,
                   c2t.version_id == version_id)
            .distinct([c2t.target_id])  # avoid duplicates - initially trust the CSC parent sample, also trust the Catalog as final arbiter
#            .distinct([c.catalogid])  # avoid duplicates - initially trust the CSC parent sample, also trust the Catalog as final arbiter
            .where
            (
                c2t.best == True
            )
        )

        return query


#-------bhm_csc_boss_dark ------ #

class BhmCscBossDarkCarton(BhmCscBaseCarton):
    '''
        SELECT * from bhm_csc AS c WHERE c.spectrograph = "BOSS" AND WHERE c.mag_i BETWEEN 17.x AND 21.x
    '''
    name = 'bhm_csc_boss_dark'
    cadence = 'bhm_csc_boss_1x4'
    instrument = 'BOSS'


    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i <  self.parameters['mag_i_max']) &
            (t.spectrograph == self.instrument)
        )
        return query

#-------bhm_csc_boss_bright ------ #

class BhmCscBossBrightCarton(BhmCscBaseCarton):
    '''
        SELECT * from bhm_csc AS c WHERE c.spectrograph = "BOSS" AND WHERE c.mag_i BETWEEN 1x.0 AND 18.x
    '''
    name = 'bhm_csc_boss_bright'
    cadence = 'bhm_csc_boss_1x1'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i <  self.parameters['mag_i_max']) &
            (t.spectrograph == self.instrument)
        )
        return query

#-------bhm_csc_apogee ------ #

class BhmCscApogeeCarton(BhmCscBaseCarton):
    '''
        SELECT * from bhm_csc AS c WHERE c.spectrograph = "APOGEE" AND WHERE c.mag_H BETWEEN 1x.x AND 1x.x
    '''
    name = 'bhm_csc_apogee'
    cadence = 'bhm_csc_apogee_3x1'
    instrument = 'APOGEE'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_h >= self.parameters['mag_h_min']) &
            (t.mag_h <  self.parameters['mag_h_max']) &
            (t.spectrograph == self.instrument)
        )
        return query
