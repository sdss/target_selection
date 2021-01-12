#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2021-01-12
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (AllWise, Catalog,
                                             CatalogToTIC_v8,
                                             EROSITASupersetAGN,
                                             EROSITASupersetClusters,
                                             EROSITASupersetStars,
                                             EROSITASupersetCompactobjects,
                                             Gaia_DR2,
                                             TIC_v8, TwoMassPSC
                                             )

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# For example
# peewee Model name ---> postgres table name
# AllWise ---> 'catalogdb.allwise'
# EROSITASupersetAGN --->'catalogdb.erosita_superset_agn'
# EROSITASupersetClusters ---> 'catalogdb.erosita_superset_clusters'
# EROSITASupersetStars ---> 'catalogdb.erosita_superset_stars'
# EROSITASupersetCompactobjects ---> 'catalogdb.erosita_superset_compactobjects'
# Gaia_DR2 --->'catalogdb.gaia_dr2_source'
# TwoMassPSC --->'catalogdb.twomass_psc'
# MWM_ROSITA_Stars_Carton
# MWM_ROSITA_Compact_Gen_Carton
# MWM_ROSITA_Compact_Var_Carton

class MWM_ROSITA_Stars_Carton(BaseCarton):
    """MWM eROSITA Stars
    Owner: Lead by MWM (with assistance from BHM?)?? Jennifer Johnson, Tom Dwelly

    Shorthand name: erosita_stars

    What is it?: Optical counterparts to eROSITA sources with
    a high likelihood to be stars. Stellar identification are based on Bayesian scheme using the Gaia dr2 and Tycho catalogs. Those objects shall be observed to increase the legacy value of eROSITA / SDSS; possible science questions relate to Halpha emission although EW are expected to be around 1 Angstrom or less.

    Simplified Description of selection criteria: "Select all objects with xmatch_metric>0.5 from catalogdb.erosita_superset_stars that are in a magnitude range amenable for SDSS-V BOSS or APOGEE spectroscopy (target_priority==1; suitability estimated from lo < Gmag < hi)"

    Wiki page: N/A

    Additional source catalogs needed: None

    Additional cross-matching needed: No

    cadence options for these targets (list all options, even though no single target will receive more than one): boss_bright_1xN or boss_dark_1xN or apogee_1xN (spectrograph, moon phase and number of exposures depends on target magnitude). No detailed time spacing required.

    Pseudo SQL (optional):  

    SELECT c.catalogid, c.ra, c.dec, twomass.j_m, twomass.h_m, twomass.k_m,
     gaia.phot_g_mean_mag, gaia.phot_bp_mean_mag, gaia.phot_rp_mean_mag, gaia.parallax, gaia.pmra, gaia.pmdec, 
     estars.target_priority, estars.xmatch_metric, estars.ero_flux, estars.ero_ra, estars.ero_dec, estars.opt_ra, estars.opt_dec 
     FROM catalog c
     INNER JOIN catalog_to_tic_v8 ctic USING (catalogid)
     INNER JOIN tic_v8 tic ON tic.id = ctic.target_id
     INNER JOIN gaia_dr2_source gaia ON gaia.source_id = tic.gaia_int
     INNER JOIN twomass_psc twomass ON twomass.designation = tic.twomass_psc
     LEFT JOIN erosita_superset_stars estars ON estars.gaia_dr2_id = gaia.source_id
    WHERE (ctic.version_id = 21) AND /* control version! */
     (ctic.best is true) AND /* and enforce unique-ish crossmatch */
     estars.target_priority = 1 AND estars.xmatch_metric > 0.5
    ;

    The results of this query can then be sorted to assign cadences using the following logic:
    bright_bright_limit = 13   # (available for modification later) 
    ir_faint_limit = 13 # (available for modification later)
    - if bright_bright_limit > gaia.phot_g_mean_mag  & twomass.h_m < ir_faint_limit:
                 cadence = bright_apogee_1x1      &&    priority = 2400
    - if bright_bright_limit < gaia.phot_g_mean_mag < 17:
                 cadence = bright_boss_1x1          &&    priority = 2400
    - if 17 < gaia.phot_g_mean_mag < 19:
                  cadence = dark_boss_1x2           &&    priority = 1920
    - if 19 < gaia.phot_g_mean_mag:
                  cadence = dark_boss_1x3          &&    priority = 1920
    """

    name = 'mwm_yso_disk_apogee'
    category = 'science'
    cadence = None  # 'apogee_bright_3x1'
    program = 'mwm_yso'
    mapper = 'MWM'
    priority = 2700

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid, Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_dr2_ra'),
                         Gaia_DR2.dec.alias('gaia_dr2_dec'),
                         TwoMassPSC.pts_key,
                         TwoMassPSC.designation.alias('twomass_psc_designation'),
                         AllWise.designation.alias('allwise_designation'),
                         Gaia_DR2.phot_g_mean_mag, Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag.alias('gaia_dr2_rp'),
                         TwoMassPSC.j_m, TwoMassPSC.h_m,
                         TwoMassPSC.k_m,
                         Gaia_DR2.parallax)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(AllWise, on=(TIC_v8.allwise == AllWise.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 13,
                        (AllWise.w1mpro - AllWise.w2mpro) > 0.25,
                        (AllWise.w2mpro - AllWise.w3mpro) > 0.50,
                        (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                        Gaia_DR2.parallax > 0.3))

        # Gaia_DR2 pweewee model class corresponds to
        # table catalogdb.gaia_dr2_source.
        #
        # All values of TIC_v8.plx (for non-null entries) are not the same as
        # values of Gaia_DR2.parallax.
        # Hence, in the above query, we cannot use TIC_v8.plx instead
        # of Gaia_DR2.parallax.

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """
        cadence options for these targets:
        boss_bright_3x1 if RP<14.76 |
        boss_bright_4x1 if RP<15.075 |
        boss_bright_5x1 if RP<15.29 |
        boss_bright_6x1 if RP<15.5
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr2_rp from " +
            " sandbox.temp_mwm_yso_disk_boss ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_rp = output[i][1]

            if(current_rp < 14.76):
                current_cadence = 'boss_bright_3x1'
            elif(current_rp < 15.075):
                current_cadence = 'boss_bright_4x1'
            elif(current_rp < 15.29):
                current_cadence = 'boss_bright_5x1'
            elif(current_rp < 15.5):
                current_cadence = 'boss_bright_6x1'
            else:
                current_cadence = None

            self.database.execute_sql(
                " update sandbox.temp_mwm_yso_disk_boss " +
                " set cadence = '" + current_cadence + "'"
                " where catalogid = " + str(current_catalogid) + ";")


