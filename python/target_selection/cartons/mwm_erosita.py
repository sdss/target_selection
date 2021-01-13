#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2021-01-12
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
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

    Shorthand name: mwm_erosita_stars

    What is it?: Optical counterparts to eROSITA sources with
    a high likelihood to be stars. Stellar identification are
    based on Bayesian scheme using the Gaia dr2 and Tycho catalogs.
    Those objects shall be observed to increase the legacy value
    of eROSITA / SDSS; possible science questions relate to
    Halpha emission although EW are expected to be around 1 Angstrom or less.

    Simplified Description of selection criteria:
    "Select all objects with xmatch_metric>0.5
    from catalogdb.erosita_superset_stars that are
    in a magnitude range amenable for SDSS-V BOSS or APOGEE spectroscopy
    (target_priority==1; suitability estimated from lo < Gmag < hi)"

    Wiki page: N/A

    Additional source catalogs needed: None

    Additional cross-matching needed: No

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    boss_bright_1xN or boss_dark_1xN or apogee_1xN
    (spectrograph, moon phase and number of exposures depends on target magnitude).
    No detailed time spacing required.

    Pseudo SQL (optional):

    SELECT c.catalogid, c.ra, c.dec, twomass.j_m, twomass.h_m, twomass.k_m,
    gaia.phot_g_mean_mag, gaia.phot_bp_mean_mag,
    gaia.phot_rp_mean_mag, gaia.parallax, gaia.pmra, gaia.pmdec,
    estars.target_priority, estars.xmatch_metric, estars.ero_flux,
    estars.ero_ra, estars.ero_dec, estars.opt_ra, estars.opt_dec
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

    Due to below, above LEFT JOIN should be just INNER JOIN.
    i.e. every row of catalogdb.erosita_superset_stars has a gaia_dr2_id.
    sdss5db=# select count(1) from catalogdb.erosita_superset_stars where gaia_dr2_id is null;
    count
    -------
        0
    (1 row)


    See cadence logic in post_process() below.
    """

    name = 'mwm_erosita_stars'
    category = 'science'
    cadence = None  # assigned in post_process()
    program = 'mwm_erosita'
    mapper = 'MWM'
    priority = None  # assigned in post_processs()

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         TwoMassPSC.j_m,
                         TwoMassPSC.h_m.alias('twomass_psc_h_m'),
                         TwoMassPSC.k_m,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_dr2_ra'),
                         Gaia_DR2.dec.alias('gaia_dr2_dec'),
                         TwoMassPSC.pts_key,
                         TwoMassPSC.designation.alias('twomass_psc_designation'),
                         Gaia_DR2.phot_g_mean_mag.alias('gaia_dr2_g'),
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         Gaia_DR2.parallax, Gaia_DR2.pmra, Gaia_DR2.pmdec,
                         EROSITASupersetStars.target_priority,
                         EROSITASupersetStars.xmatch_metric,
                         EROSITASupersetStars.ero_flux,
                         EROSITASupersetStars.ero_ra,
                         EROSITASupersetStars.ero_dec,
                         EROSITASupersetStars.opt_ra,
                         EROSITASupersetStars.opt_dec)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(Gaia_DR2)
                 .join(EROSITASupersetStars,
                       on=(Gaia_DR2.source_id == EROSITASupersetStars.gaia_dr2_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        EROSITASupersetStars.target_priority == 1,
                        EROSITASupersetStars.xmatch_metric > 0.5))

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
        This post_process() method does two steps:
        (a) The source with the largest xmatch_metric is selected.
        
        (b) The results of the above query can then be sorted
        to assign cadences using the following logic:
        bright_bright_limit = 13   # (available for modification later)
        ir_faint_limit = 13 # (available for modification later)
        - if bright_bright_limit > gaia.phot_g_mean_mag  &
             twomass.h_m < ir_faint_limit:
                 cadence = bright_apogee_1x1      &&    priority = 2400
        - if bright_bright_limit < gaia.phot_g_mean_mag < 17:
                 cadence = bright_boss_1x1          &&    priority = 2400
        - if 17 < gaia.phot_g_mean_mag < 19:
                  cadence = dark_boss_1x2           &&    priority = 1920
        - if 19 < gaia.phot_g_mean_mag:
                  cadence = dark_boss_1x3          &&    priority = 1920
        """

        bright_bright_limit = 13
        ir_faint_limit = 13

        # select source with largest xmatch_metric
        self.database.execute_sql("update sandbox.temp_mwm_erosita_stars " +
                                  "set selected = false")

        cursor = self.database.execute_sql(
            "select catalogid, max(xmatch_metric) from " +
            " sandbox.temp_mwm_erosita_stars " +
            " group by catalogid ;")
       # TODO
        output = cursor.fetchall()

        list_of_catalog_id = [0] * len(output)
        current_count = 0
        current_target = 0
        for i in range(len(output)):
            current_xmatch_metric = output[i][1]
            if(count_count < 1):
                curent_count = current_count + 1
                list_of_catalog_id[current_target] = output[i][0]
                current_target = current_target + 1

        max_target = current_target
        for k in range(max_target + 1):
            self.database.execute_sql(
                " update sandbox.temp_mwm_erosita_stars set selected = true " +
                " where catalogid = " + str(list_of_catalog_id[k]) + ";")


        # Set cadence and priority

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr2_g, twomass_psc_h_m from " +
            " sandbox.temp_mwm_erosita_stars ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_g = output[i][1]
            current_h = output[i][2]

            if((current_g < bright_bright_limit) and
               (current_h < ir_faint_limit)):
                current_cadence = 'bright_apogee_1x1'
                current_priority = 2400
            elif((bright_bright_limit < current_g) and
                 (current_g < 17)):
                current_cadence = 'bright_boss_1x1'
                current_priority = 2400
            elif((17 < current_g) and
                 (current_g < 19)):
                current_cadence = 'dark_boss_1x2'
                current_priority = 1920
            elif(19 < current_g):
                current_cadence = 'dark_boss_1x3'
                current_priority = 1920
            else:
                current_cadence = None

            self.database.execute_sql(
                " update sandbox.temp_mwm_erosita_stars " +
                " set cadence = '" + current_cadence + "'"
                " where catalogid = " + str(current_catalogid) + ";")

            self.database.execute_sql(
                " update sandbox.temp_mwm_erosita_stars " +
                " set priority = '" + current_priority + "'"
                " where catalogid = " + str(current_catalogid) + ";")
