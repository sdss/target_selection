#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_gua.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn
# import sdssdb

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToTIC_v8,
    TIC_v8,
    Gaia_DR2,
    # CatalogToGaia_unWISE_AGN,  # <-- this is old and does not work
    # CatalogToSDSS_DR16_SpecObj, # <-- rely on a q3c match in v0.5
    SDSS_DR16_SpecObj,
    Gaia_unWISE_AGN,
)
from target_selection.cartons.base import BaseCarton

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms
#
# This module provides the following BHM cartons:
# bhm_gua_dark
# bhm_gua_bright
#
# Updated by TD on 27/01/2021 to add v0.5 improvements
# switch to a q3c join to sdss_dr16_specobj

'''
  pseudo-SQL

    both Cartons:
        SELECT catalogid,{derived fields}
        FROM catalog AS c
        JOIN catalog_to_tic AS c2tic ON ...
        JOIN tic_v8 AS tic ON ...
        JOIN gaia_dr2_source AS g ON ...
        JOIN gaia_unwise_agn AS t ON g.source_id = t.gaia_sourceid
        LEFT OUTER JOIN sdss_dr16_specobj AS s
          ON ( q3c_join(s.ra,s.dec,c.ra,c.dec,{match_radius_spectro})
               AND s.zwarning = 0
               AND s.snmedianl > 2.0
               AND s.zerr < 0.01
               AND s.scienceprimary > 0
             )
        WHERE
                c.version_id = {version_id}
            AND c2tic.version_id = {version_id}
            AND c2tic.best = True
            AND t.prob_rf > 0.8
            AND s.specobjid = NULL

    # bhm_gua_dark carton:
            AND ( t.g > 16.5 AND
                  t.rp > 16.5 AND
                  (t.g < 21.2 OR t.rp < 21.0 )
                )

    # bhm_gua_bright carton:
            AND ( t.g > 13.0 AND
                  t.rp > 13.5 AND
                  (t.g < 18.5 OR t.rp < 18.5)
                )
'''


class BhmGuaBaseCarton(BaseCarton):
    '''
    Parent class that provides the basic selections for both Gaia UnWISE AGN cartons
    To be sub-classed, not to be called directly.

    To get from Catalog to GUA we join tables via :
    Catalog -> CatalogToTIC_v8 -> Gaia_DR2 -> Gaia_unWISE_AGN

    Cadence+priority per target depends on brightness
    Keep the two bhm_gua_* cartons separate to enable overlap in brightness ranges
    '''

    name = 'bhm_gua_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_filler'
    tile = False
    priority = None
    cadence = None
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        # ## c2t = CatalogToGaia_unWISE_AGN.alias() - deprecated - but leave this as a reminder
        c2tic = CatalogToTIC_v8.alias()
        tic = TIC_v8.alias()
        # c2s = CatalogToSDSS_DR16_SpecObj.alias()
        g = Gaia_DR2.alias()
        t = Gaia_unWISE_AGN.alias()

        # Keep things simple by trusting the specprimary labels in the
        # SDSS specObj file - consider only the 'best' spectrum per sky position
        ss = SDSS_DR16_SpecObj.alias()
        s = ss.select(
            ss.specobjid.alias('specobjid'),
            ss.ra.alias('ra'),
            ss.dec.alias('dec'),
        ).where
        (
            (ss.snmedian >= self.parameters['spec_sn_thresh']),
            (ss.zwarning == 0),
            (ss.zerr <= self.parameters['spec_z_err_thresh']),
            (ss.zerr > 0.0),
            (ss.scienceprimary > 0),
        ).alias('s')

        # s = SDSS_DR16_SpecObj.alias()

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get('priority', 10000)))
        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float')
        inertial = peewee.Value(True)
        cadence = peewee.Value(self.parameters['cadence'])
        instrument = peewee.Value(self.instrument)

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0

        query = (
            c.select(
                c.catalogid,
                c.ra,
                c.dec,
                t.gaia_sourceid,
                priority.alias('priority'),
                value.alias('value'),
                inertial.alias('inertial'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                # rely on the centralised magnitude routines
            )
            .join(c2tic)
            .join(tic)
            .join(g)
            .join(t, on=(g.source_id == t.gaia_sourceid))
            # sky join to 'good' primary spectra in SDSS-DR16
            .join(
                s,
                JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(s.c.ra, s.c.dec,
                                c.ra, c.dec,
                                match_radius_spectro) &

                )
            )
            # then reject any GUA targets with existing good DR16 spectroscopy
            .where(s.c.specobjid.is_null(True))
            # standard selection that chooses correct catalogdb version etc
            .where(c.version_id == version_id,
                   c2tic.version_id == version_id,
                   c2tic.best >> True)
            .where(
                (t.prob_rf >= self.parameters['prob_rf_min']),
                (t.g >= self.parameters['mag_g_min']),
                (t.rp >= self.parameters['mag_rp_min']),
                (
                    (t.g < self.parameters['mag_g_max']) |
                    (t.rp < self.parameters['mag_rp_max'])
                ),
            )
            # avoid duplicates - trust the gaia ids in the GUA parent sample
            .distinct([t.gaia_sourceid])
            # .switch(c)
            # .join(c2s, JOIN.LEFT_OUTER)
            # .join(s, JOIN.LEFT_OUTER)
            # .where(
            #     (s.specobjid.is_null()) |
            #     (s.zwarning != 0 ) |
            #     (s.snmedian < self.parameters['spec_sn_thresh']) |
            #     (s.zerr >  self.parameters['spec_z_err_thresh'])
            # )
        )

        return query


class BhmGuaDarkCarton(BhmGuaBaseCarton):
    '''
    -------  bhm_gua_dark   ------
    SQL as above plus:

        AND ( gua.g > 16.x AND gua.rp > 16.x)
    '''
    name = 'bhm_gua_dark'


class BhmGuaBrightCarton(BhmGuaBaseCarton):
    '''
    ------  bhm_gua_bright   ------
    SQL as above plus:

        AND ( gua.g < 18.x OR gua.rp < 18.x)

    '''
    name = 'bhm_gua_bright'


#
# Exporting from the temp table - from within psql
#
# psql> \copy (SELECT * FROM sandbox.temp_bhm_gaia_unwise_agn_dark)
# psql>   TO 'bhm_gaia_unwise_agn_dark.csv' with csv header
# psql> \copy (SELECT * FROM sandbox.temp_bhm_gaia_unwise_agn_bright)
# psql>   TO 'bhm_gaia_unwise_agn_bright.csv' with csv header
#
# Back in bash
# #> for F in bhm_gaia_unwise_agn_*.csv; do
# #>     stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic;
# #> done
#
