#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-03
# @Filename: bhm_aqmes.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import peewee
import sdssdb
from astropy.io import fits
import pkg_resources

#debug
#import os


from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             SDSS_DR13_PhotoObj,
                                             SDSS_DR16_SpecObj,
                                             SDSS_DR14_QSO,
                                             CatalogToSDSS_DR13_PhotoObj)

# when DR16Q is available
#from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
#                                             SDSS_DR16_QSO )
#                                             CatalogToSDSS_DR16_QSO,)


from target_selection.cartons.base import BaseCarton

radius_apo = 1.49 # degrees


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-AQMES-medium-inprogress

## This provides the following BHM cartons:
# bhm_aqmes_med
# bhm_aqmes_med_faint
# bhm_aqmes_wide2
# bhm_aqmes_wide3
# bhm_aqmes_wide2_faint
# bhm_aqmes_wide3_faint
# bhm_aqmes_bonus_dark
# bhm_aqmes_bonus_bright



class BhmAqmesBaseCarton(BaseCarton):
    ''' Parent class that provides the underlying selections for all AQMES cartons'''

    name = 'bhm_aqmes_base'
    category = 'science'
    mapper = 'BHM'
    program = 'BHM-AQMES'
    tile = False
    priority = None
    alias_c = None
    alias_t = None

    # read the AQMES field centres from a fits file and convert to a list of dicts
    def get_fieldlist(self):
        stub = self.parameters.get('fieldlist', None)
        if stub is None or stub == '' or stub == 'None':
            return None

        #filename = pkg_resources.resource_filename( __name__, stub)
        filename = pkg_resources.resource_filename('target_selection', stub)
        assert len(filename) > 0

        try:
            hdul = fits.open(filename)
        except:
            raise Exception(f"Failed to find/open fieldlist file: {filename}")

        assert len(hdul[1].data) > 0

        # choose the correct subset of fields based on the cadence name and form a list of dicts
        try:
            fieldlist = [ {'racen': r['RACEN'],
                           'deccen': r['DECCEN'],
                           'radius': radius_apo, }
                          for r in hdul[1].data
                          if r['CADENCE'] == self.cadence
            ]
        except:
            raise Exception(f"Error interpreting contents of fieldlist file: {filename}")

        assert len(fieldlist) > 0

        return fieldlist

    def append_spatial_query(self, query, fieldlist):
        '''extend the peewee query using a list of field centres'''
        if fieldlist is None :
            return query
        elif len(fieldlist) == 0 :
            return query

        q = False
        for f in fieldlist:
            q = ( q | peewee.fn.q3c_radial_query(self.alias_c.ra,
                                                 self.alias_c.dec,
                                                 f['racen'],
                                                 f['deccen'],
                                                 f['radius']))
        return query.where(q)


    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2p = CatalogToSDSS_DR13_PhotoObj.alias()
        p = SDSS_DR13_PhotoObj.alias()
        s = SDSS_DR16_SpecObj.alias()
        t = SDSS_DR14_QSO.alias()
        self.alias_c = c
        self.alias_t = t

        # set the Carton priority+values here - read from yaml
        target_priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')
        target_value = peewee.Value(float(self.parameters.get('value', 1.0))).alias('value')
        pmra =  peewee.Value(float(0.0)).alias('pmra')
        pmdec = peewee.Value(float(0.0)).alias('pmdec')

        query = (
            c
            .select(c.catalogid,
                    target_priority,
                    target_value,
                    pmra,
                    pmdec,
                    t.psfmag[1].alias('magnitude_g'),
                    t.psfmag[2].alias('magnitude_r'),
                    t.psfmag[3].alias('magnitude_i'),
                    t.psfmag[4].alias('magnitude_z'),
            )
            .join(c2p)
            .join(p)
            .join(s, on=(p.objid == s.bestobjid))
            .join(t, on=((s.plate == t.plate) & (s.mjd == t.mjd) & (s.fiberid == t.fiberid)))
            .where(c.version_id == version_id,
                   c2p.version_id == version_id,
                   c2p.best == True)
            .distinct([t.pk])   # avoid duplicates - trust the QSO parent sample
            .where
            (
                (t.psfmag[3] >= self.parameters['mag_i_min']),
                (t.psfmag[3] <  self.parameters['mag_i_max']),
                #            (t.z >= self.parameters['redshift_min']),
                #            (t.z <= self.parameters['redshift_max']),
            )
        )
        query = self.append_spatial_query(query, self.get_fieldlist())


        return query



#-------AQMES medium ------ #

class BhmAqmesMedCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE  psfmag_i BETWEEN 16.x AND 19.1
    '''
    name = 'bhm_aqmes_med'
    cadence = 'bhm_aqmes_medium_12x4'
    #cadence = 'dummy_cadence'

# add something like the following if want to add carton-specific selections
#    def build_query(self, version_id, query_region=None):
#        query = super().build_query(version_id, query_region)
#        query = query.where( # .... add extra terms here
#        )
#        return query


class BhmAqmesMedFaintCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE  psfmag_i BETWEEN 19.1 AND 21.0
    '''
    name = 'bhm_aqmes_med_faint'
    cadence = 'bhm_aqmes_medium_12x4'

#-------AQMES medium ------ #


#-------AQMES wide ------ #

class BhmAqmesWide3Carton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE psfmag_i BETWEEN 16.x AND 19.1
    '''
    name = 'bhm_aqmes_wide3'
    cadence = 'bhm_aqmes_wide_3x4'


class BhmAqmesWide3FaintCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE psfmag_i BETWEEN 19.1 AND 21.0
    '''
    name = 'bhm_aqmes_wide3_faint'
    cadence = 'bhm_aqmes_wide_3x4'


class BhmAqmesWide2Carton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE psfmag_i BETWEEN 16.x AND 19.1
    '''
    name = 'bhm_aqmes_wide2'
    cadence = 'bhm_aqmes_wide_2x4'


class BhmAqmesWide2FaintCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE psfmag_i BETWEEN 19.1 AND 21.0
    '''
    name = 'bhm_aqmes_wide2_faint'
    cadence = 'bhm_aqmes_wide_2x4'

#-------AQMES wide ------ #


#-------AQMES bonus ------ #

class BhmAqmesBonusDarkCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE psfmag_i BETWEEN 16.x AND 21.5
    '''
    name = 'bhm_aqmes_bonus_dark'
    cadence = 'bhm_spiders_1x4'


class BhmAqmesBonusBrightCarton(BhmAqmesBaseCarton):
    '''
    SELECT * FROM sdss_dr1x_qso WHERE psfmag_i BETWEEN 14.0 AND 18.0
    '''
    name = 'bhm_aqmes_bonus_bright'
    cadence = 'csc_bright_apogee_3x1'   # could add a new cadence name for this

#-------AQMES bonus ------ #








##################################################################################
### Debug stuff

'''
SQL:
select c.*,s.specobjid,t.class,t.psfmag[2] as psfmag_g,t.psfmag[4] AS psfmag_i,t.z,t.zwarning,t.plate,t.mjd,t.fiberid,q3c_dist(c.ra,c.dec,t.ra,t.dec)*3600.::FLOAT as sep from catalog AS C INNER JOIN catalog_to_sdss_dr13_photoobj as c2p ON c.catalogid = c2p.catalogid INNER JOIN sdss_dr13_photoobj AS p ON c2p.target_id = p.objid INNER JOIN sdss_dr16_specobj as s on p.objid = s.bestobjid INNER JOIN sdss_dr14_qso AS t ON (s.plate = t.plate AND s.mjd = t.mjd AND s.fiberid = t.fiberid ) WHERE (q3c_radial_query(c.ra, c.dec, 10.0, 25.0, 0.1) AND c.version_id = 13 AND (p.resolvestatus & 256) != 0 AND t.psfmag[4] < 19.1);

#python:
t = SDSS_DR14_QSO.alias()
for f in t._meta.fields:
     print (f)

'''

'''
# for testing do domething like this

import peewee
import sdssdb
from sdssdb.peewee.sdss5db.catalogdb import database
database.set_profile('tunnel_operations')
from target_selection.cartons.bhm_aqmes import *
c = BhmAqmesMedCarton(targeting_plan='0.1.0-beta.1')
q = c.build_query(version_id=13)
for r in q.limit(5).namedtuples():
    print(r)

'''


###################################################################################
