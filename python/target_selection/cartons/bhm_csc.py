#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# ## flake8: noqa
# isort: skip_file

import peewee
# import sdssdb

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             BHM_CSC,
                                             CatalogToBHM_CSC)

from target_selection.cartons.base import BaseCarton


# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ChandraSourceCatalogue(CSC)-inprogress  # noqa: E501

# This module provides the following BHM cartons:
#   bhm_csc_boss_dark
#   bhm_csc_boss_bright
#   bhm_csc_apogee


class BhmCscBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for all CSC cartons
    To be sub-classed, not typically to be called directly.

    SELECT DISTINCT ON (t.pk)
      c.catalogid, c.lead, c.ra, c.dec,
      q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
      t.mag_g, t.mag_r, t.mag_i
    FROM catalog AS c
    JOIN catalog_to_bhm_csc AS c2t
      ON c.catalogid = c2t.catalogid
    JOIN bhm_csc AS t
      ON c2t.target_id = t.pk
    WHERE
          c.version_id = 21
      AND c2t.version_id = 21
      AND c2t.best IS TRUE
      AND t.spectrograph = 'BOSS'
      AND t.mag_i BETWEEN 16.0 AND 22.0
    ORDER BY t.pk
    LIMIT 20;


    # check how many csc targets have positions defined by the  bhm_csc.oir_ra,oir_dec fields
    # these may actually be Chandra X-ray coords and not optical positions!

    SELECT
        count(*) as count,
        AVG(sep) as mean_sep,
        MIN(sep) as min_sep,
        MAX(sep) as max_sep,
        AVG(mag_r) as mean_mag_r,
        MIN(mag_r) as min_mag_r,
        MAX(mag_r) as max_mag_r
    FROM (
        SELECT DISTINCT ON (t.pk)
           c.catalogid, c.lead, c.ra, c.dec,
          q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
          t.mag_g, t.mag_r, t.mag_i
        FROM catalog AS c
        JOIN catalog_to_bhm_csc AS c2t
          ON c.catalogid = c2t.catalogid
        JOIN bhm_csc AS t
          ON c2t.target_id = t.pk
        WHERE
              c.version_id = 21
          AND c2t.version_id = 21
          AND c2t.best IS TRUE
          AND t.spectrograph = 'BOSS'
          AND t.mag_i BETWEEN 16.0 AND 22.0
        ORDER BY t.pk
        ) AS subq
      WHERE subq.lead = 'bhm_csc';

      SELECT DISTINCT ON (t.pk)
          c.catalogid, c.lead, c.ra, c.dec,
          q3c_dist(c.ra, c.dec,t.oir_ra, t.oir_dec)*3600.0 as sep,
          t.mag_g, t.mag_r, t.mag_i
        INTO sandbox.temp_td_csc_subset
        FROM catalog AS c
        JOIN catalog_to_bhm_csc AS c2t
          ON c.catalogid = c2t.catalogid
        JOIN bhm_csc AS t
          ON c2t.target_id = t.pk
        WHERE
              c.version_id = 21
          AND c2t.version_id = 21
          AND c2t.best IS TRUE
          AND t.spectrograph = 'BOSS'
          AND t.mag_i BETWEEN 16.0 AND 22.0
          AND c.lead = 'bhm_csc'
        ORDER BY t.pk;



    '''

    name = 'bhm_csc_base'
    base_name = 'bhm_csc_base'
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
        instrument = peewee.Value(self.instrument).alias('instrument')
        cadence = peewee.Value(self.cadence).cast('string').alias('cadence')

        # TODO need to process the bhm_csc.[g,r,i,z,h] magnitudes to deal with zeros/nulls

        query = (
            c.select(
                c.catalogid,
                c.ra,  # extra
                c.dec,  # extra
                priority,
                value,
                cadence,
                instrument,
                t.mag_g.alias('g'),
                t.mag_r.alias('r'),
                t.mag_i.alias('i'),
                t.mag_z.alias('z'),
                t.mag_h.alias('h'),
                t.oir_ra.alias('csc_ra'),
                t.oir_dec.alias('csc_dec'),
            )
            .join(c2t)
            .join(t)
            .where(
                c.version_id == version_id,
                c2t.version_id == version_id,
                c2t.best >> True)
            .distinct([c2t.target_id])  # avoid duplicates - initially trust the CSC parent sample,
            .where
            (
                t.spectrograph == self.instrument
            )
        )

        return query


class BhmCscBossDarkCarton(BhmCscBaseCarton):
    '''
    SELECT * from bhm_csc AS c
    WHERE c.spectrograph = "BOSS"
       ND c.mag_i BETWEEN 17.x AND 21.x
    '''
    name = 'bhm_csc_boss_dark'
    cadence = 'dark_1x4'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i < self.parameters['mag_i_max'])
        )
        return query


class BhmCscBossBrightCarton(BhmCscBaseCarton):
    '''
    SELECT * from bhm_csc AS c
    WHERE c.spectrograph = "BOSS"
      AND c.mag_i BETWEEN 1x.0 AND 18.x
    '''
    name = 'bhm_csc_boss_bright'
    cadence = 'bright_1x1'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_i >= self.parameters['mag_i_min']) &
            (t.mag_i < self.parameters['mag_i_max'])
        )
        return query


class BhmCscApogeeCarton(BhmCscBaseCarton):
    '''
    SELECT * from bhm_csc AS c
    WHERE c.spectrograph = "APOGEE"
      AND c.mag_H BETWEEN 1x.x AND 1x.x
    '''
    name = 'bhm_csc_apogee'
    cadence = 'bright_3x1'
    instrument = 'APOGEE'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.mag_h >= self.parameters['mag_h_min']) &
            (t.mag_h < self.parameters['mag_h_max'])
        )
        return query

#
# '''
# Exporting from the temp table
#
# \copy (SELECT * FROM sandbox.temp_bhm_csc_boss_dark)  TO 'bhm_csc_boss_dark.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_csc_boss_bright)
#      TO 'bhm_csc_boss_bright.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_csc_apogee)  TO 'bhm_csc_apogee.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_td_csc_subset)  TO 'temp_td_csc_subset.csv' with csv header
# for F in bhm_csc_*.csv; do
#    stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic;
# done
#
# '''
