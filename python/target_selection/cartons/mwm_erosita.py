#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2021-01-12
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             EROSITASupersetCompactobjects,
                                             EROSITASupersetStars, Gaia_DR2,
                                             TIC_v8, TwoMassPSC)

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


class MWM_EROSITA_Stars_Carton(BaseCarton):
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
                 .distinct(CatalogToTIC_v8.catalogid)
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
        The results of the above query can then be sorted
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
               (current_h is not None) and
               (current_h < ir_faint_limit)):
                current_cadence = 'bright_apogee_1x1'
                current_priority = 2400
            elif((bright_bright_limit < current_g) and (current_g < 17)):
                current_cadence = 'bright_boss_1x1'
                current_priority = 2400
            elif((17 < current_g) and (current_g < 19)):
                current_cadence = 'dark_boss_1x2'
                current_priority = 1920
            elif(19 < current_g):
                current_cadence = 'dark_boss_1x3'
                current_priority = 1920
            else:
                current_cadence = None
                current_priority = None

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_erosita_stars " +
                    " set cadence = '" + current_cadence + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_erosita_stars " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_EROSITA_Compact_Gen_Carton(BaseCarton):
    """MWM eROSITA Compact-Objects (carton *_gen)

    Owner: eROSITA Compact objects WG Axel Schwope

    Shorthand name: mwm_erosita_compact_gen

    What is it?: eROSITA selected pointlike X-ray sources with
    likely Gaia counterpart (currently drawn from DR2),
    with high likelihood of being a galactic object
    according to proper motion and parallax,
    and not being an ordinary coronal emitter,
    also being (optically) reasonably faint,
    taking together being a candidate of a compact binary,
    i.e. accreting compact object.
    The original selection criteria that were applied
    to all X-ray/Gaia pairs are described on page Compact Binaries.
    These form the erosita_superset_compactobjects.

    Simplified Description of selection criteria:
    "Select all targets from erosita_superset_compactobjects
    that have minimal distance to eROSITA eRASS1 position BOSS spectroscopy"

    Wiki page: link or N/A

    Additional source catalogs needed: None

    Additional cross-matching needed: None

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    minimum exposure is 1x1 for bright 16<g < 17 mag),
    1x2 for medium (17<g<18), 1xn for the fainter targets (n as large as possible)

    Pseudo SQL (optional):

    select * from erosita_superset_compactobjects
    where xmatch_version = 'ASJK_0212020_select1uniq'

    reveals 78.252 targets
    Should one of these targets compete with one chosen
    from *_var (below), give preference to the *_var object

    Implementation:

    SELECT c.catalogid, c.ra, c.dec, gaia.phot_g_mean_mag,
    gaia.phot_bp_mean_mag, gaia.phot_rp_mean_mag, gaia.parallax,
    gaia.pmra, gaia.pmdec,
    estars.target_priority, estars.xmatch_metric,
    estars.ero_flux, estars.ero_ra, estars.ero_dec,
    estars.opt_ra, estars.opt_dec
    FROM catalog c
    INNER JOIN catalog_to_tic_v8 ctic USING (catalogid)
    INNER JOIN tic_v8 tic ON tic.id = ctic.target_id
    INNER JOIN gaia_dr2_source gaia ON gaia.source_id = tic.gaia_int
    LEFT JOIN erosita_superset_compactobjects estars ON estars.gaia_dr2_id = gaia.source_id
    WHERE (ctic.version_id = 21) AND /* control version! */
    (ctic.best is true) AND /* and enforce unique-ish crossmatch */
    estars.xmatch_version = 'ASJK_0212020_select1uniq'

    See cadence logic in post_process() below.
    """

    name = 'mwm_erosita_compact_gen'
    category = 'science'
    cadence = None  # assigned in post_process()
    program = 'mwm_erosita'
    mapper = 'MWM'
    priority = None  # assigned in post_processs()

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_dr2_ra'),
                         Gaia_DR2.dec.alias('gaia_dr2_dec'),
                         Gaia_DR2.phot_g_mean_mag.alias('gaia_dr2_g'),
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         Gaia_DR2.parallax, Gaia_DR2.pmra, Gaia_DR2.pmdec,
                         EROSITASupersetCompactobjects.target_priority,
                         EROSITASupersetCompactobjects.xmatch_metric,
                         EROSITASupersetCompactobjects.ero_flux,
                         EROSITASupersetCompactobjects.ero_ra,
                         EROSITASupersetCompactobjects.ero_dec,
                         EROSITASupersetCompactobjects.opt_ra,
                         EROSITASupersetCompactobjects.opt_dec)
                 .distinct(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .join(EROSITASupersetCompactobjects,
                       on=(Gaia_DR2.source_id == EROSITASupersetCompactobjects.gaia_dr2_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        EROSITASupersetCompactobjects.xmatch_version ==
                        'ASJK_0212020_select1uniq'))

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
        The results of the above query can then be sorted
        to assign cadences using the following logic:
         #------ then use following logic to assign cadences

        bright_bright_limit = 13   # (available for modification later)

        - if bright_bright_limit < gaia.phot_g_mean_mag < 17:
               cadence = bright_boss_1x1   &&   priority = 2400

        - if 17 < gaia.phot_g_mean_mag < 19:
               cadence = dark_boss_1x2     &&   priority = 1910

        - if 19 < gaia.phot_g_mean_mag:
        cadence = dark_boss_1x3     &&   priority = 1910

        Note: For the case gaia.phot_g_mean_mag < bright_bright_limit
              the cadence is None
        """

        bright_bright_limit = 13

        # Set cadence and priority

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr2_g from " +
            " sandbox.temp_mwm_erosita_compact_gen ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_g = output[i][1]

            if((current_g < bright_bright_limit)):
                current_cadence = None
                current_priority = None
            elif((bright_bright_limit < current_g) and (current_g < 17)):
                current_cadence = 'bright_boss_1x1'
                current_priority = 2400
            elif((17 < current_g) and (current_g < 19)):
                current_cadence = 'dark_boss_1x2'
                current_priority = 1910
            elif(19 < current_g):
                current_cadence = 'dark_boss_1x3'
                current_priority = 1910
            else:
                current_cadence = None
                current_priority = None

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_erosita_compact_gen " +
                    " set cadence = '" + current_cadence + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_erosita_compact_gen " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_EROSITA_Compact_Var_Carton(BaseCarton):
    """MWM eROSITA Compact-Objects (carton *_var)

    Owner: eROSITA Compact objects WG, Axel Schwope

    Shorthand name: mwm_erosita_compact_var

    What is it?: eROSITA selected pointlike X-ray sources with
    likely Gaia counterpart (currently drawn from DR2),
    with high likelihood of being a galactic object
    according to proper motion, parallax, and variability,
    and being reasonably faint, taking together
    being a candidate of a compact binary,
    i.e. accreting compact object.
    The original selection criteria that were applied to
    all X-ray/Gaia pairs are described on page Compact Binaries.
    These form the erosita_superset_compactobjects

    Simplified Description of selection criteria:
    "Select all targets from erosita_superset_compactobjects
    that have maximum varproxy for BOSS spectroscopy"

    Wiki page: link or N/A

    Additional source catalogs needed:

    Additional cross-matching needed:

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    1x1 for bright 16<g < 17 mag), 1x2 for medium (17<g<18),
    1xn for the fainter targets (n as large as possible)

    Pseudo SQL (optional):

    select * from erosita_superset_compactobjects
    where xmatch_version = 'ASJK_0212020_select2univ'

        reveals 17.058 targets
        Should one of these targets compete with
        one chosen from *_gen (above), give preference to the *_var object

    Implementation:

    SELECT c.catalogid, c.ra, c.dec, gaia.phot_g_mean_mag,
    gaia.phot_bp_mean_mag, gaia.phot_rp_mean_mag,
    gaia.parallax, gaia.pmra, gaia.pmdec,
    estars.target_priority, estars.xmatch_metric,
    estars.ero_flux, estars.ero_ra, estars.ero_dec,
    estars.opt_ra, estars.opt_dec
    FROM catalog c
    INNER JOIN catalog_to_tic_v8 ctic USING (catalogid)
    INNER JOIN tic_v8 tic ON tic.id = ctic.target_id
    INNER JOIN gaia_dr2_source gaia ON gaia.source_id = tic.gaia_int
    LEFT JOIN erosita_superset_compactobjects estars ON estars.gaia_dr2_id = gaia.source_id
    WHERE (ctic.version_id = 21) AND /* control version! */
    (ctic.best is true) AND /* and enforce unique-ish crossmatch */
     estars.xmatch_version = 'ASJK_0212020_select2univ'

    See cadence logic in post_process() below.
    """

    name = 'mwm_erosita_compact_var'
    category = 'science'
    cadence = None  # assigned in post_process()
    program = 'mwm_erosita'
    mapper = 'MWM'
    priority = None  # assigned in post_processs()

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_dr2_ra'),
                         Gaia_DR2.dec.alias('gaia_dr2_dec'),
                         Gaia_DR2.phot_g_mean_mag.alias('gaia_dr2_g'),
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         Gaia_DR2.parallax, Gaia_DR2.pmra, Gaia_DR2.pmdec,
                         EROSITASupersetCompactobjects.target_priority,
                         EROSITASupersetCompactobjects.xmatch_metric,
                         EROSITASupersetCompactobjects.ero_flux,
                         EROSITASupersetCompactobjects.ero_ra,
                         EROSITASupersetCompactobjects.ero_dec,
                         EROSITASupersetCompactobjects.opt_ra,
                         EROSITASupersetCompactobjects.opt_dec)
                 .distinct(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .join(EROSITASupersetCompactobjects,
                       on=(Gaia_DR2.source_id == EROSITASupersetCompactobjects.gaia_dr2_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        EROSITASupersetCompactobjects.xmatch_version ==
                        'ASJK_0212020_select2univ'))

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
        The results of the above query can then be sorted
        to assign cadences using the following logic:
         #------ then use following logic to assign cadences

        bright_bright_limit = 13   # (available for modification later)

        - if bright_bright_limit < gaia.phot_g_mean_mag < 17:
               cadence = bright_boss_1x1   &&   priority = 2400

        - if 17 < gaia.phot_g_mean_mag < 19:
               cadence = dark_boss_1x2     &&   priority = 1900

        - if 19 < gaia.phot_g_mean_mag:
               cadence = dark_boss_1x3     &&   priority = 1900

        Note: For the case gaia.phot_g_mean_mag < bright_bright_limit
              the cadence is None

        """

        bright_bright_limit = 13

        # Set cadence and priority

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr2_g from " +
            " sandbox.temp_mwm_erosita_compact_var ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_g = output[i][1]

            if((current_g < bright_bright_limit)):
                current_cadence = None
                current_priority = None
            elif((bright_bright_limit < current_g) and (current_g < 17)):
                current_cadence = 'bright_boss_1x1'
                current_priority = 2400
            elif((17 < current_g) and (current_g < 19)):
                current_cadence = 'dark_boss_1x2'
                current_priority = 1900
            elif(19 < current_g):
                current_cadence = 'dark_boss_1x3'
                current_priority = 1900
            else:
                current_cadence = None
                current_priority = None

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_erosita_compact_var " +
                    " set cadence = '" + current_cadence + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_erosita_compact_var " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")
