#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# flake8: noqa
# isort: skip_file

import peewee
import sdssdb

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             BHM_CSC,
                                             CatalogToBHM_CSC)


from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ChandraSourceCatalogue(CSC)-inprogress

## This module provides the following BHM cartons:
# bhm_csc_boss-dark
# bhm_csc_boss-bright
# bhm_csc_apogee



class BhmCscBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for all CSC cartons
    To be sub-classed, not typically to be called directly.
    '''

    name = 'bhm_csc-base'
    base_name = 'bhm_csc-base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_csc'
    tile = False
    priority = None
    alias_t = None
    instrument = None

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2t = CatalogToBHM_CSC.alias()
        t = BHM_CSC.alias()
        self.alias_t = t

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float').alias('value')

        query = (
            c
            .select(c.catalogid,
                    priority,
                    value,
                    t.mag_g.alias('g'),
                    t.mag_r.alias('r'),
                    t.mag_i.alias('i'),
                    t.mag_z.alias('z'),
                    t.mag_h.alias('h'),
            )
            .join(c2t)
            .join(t)
            .where(c.version_id == version_id,
                   c2t.version_id == version_id,
                   c2t.best == True)
            .distinct([c2t.target_id])  # avoid duplicates - initially trust the CSC parent sample,
            .where
            (
                t.spectrograph == self.instrument
            )
        )

        return query


#-------bhm_csc_boss_dark ------ #

class BhmCscBossDarkCarton(BhmCscBaseCarton):
    '''
        SELECT * from bhm_csc AS c WHERE c.spectrograph = "BOSS" AND WHERE c.mag_i BETWEEN 17.x AND 21.x
    '''
    name = 'bhm_csc_boss-dark'
    cadence = 'bhm_csc_boss_1x4'
    instrument = 'BOSS'


    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i <  self.parameters['mag_i_max'])
        )
        return query

#-------bhm_csc_boss_bright ------ #

class BhmCscBossBrightCarton(BhmCscBaseCarton):
    '''
        SELECT * from bhm_csc AS c WHERE c.spectrograph = "BOSS" AND WHERE c.mag_i BETWEEN 1x.0 AND 18.x
    '''
    name = 'bhm_csc_boss-bright'
    cadence = 'bhm_csc_boss_1x1'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i <  self.parameters['mag_i_max'])
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
            (t.mag_h <  self.parameters['mag_h_max'])
        )
        return query

'''
Exporting from the temp table

\copy (SELECT * FROM sandbox.temp_bhm_csc_boss_dark)  TO '/home/tdwelly/scratch/targetdb/bhm_csc_boss_dark.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_csc_boss_bright)  TO '/home/tdwelly/scratch/targetdb/bhm_csc_boss_bright.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_csc_apogee)  TO '/home/tdwelly/scratch/targetdb/bhm_csc_apogee.csv' with csv header

for F in bhm_csc_*.csv; do   stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic; done
m
'''
