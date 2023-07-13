#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-22
# @Filename: mwm_halo_gaia_dr3.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import math

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (AllWise, Catalog,
                                             CatalogToAllWise,
                                             CatalogToGaia_DR3, Gaia_DR3,
                                             Gaia_dr3_vari_rrlyrae,
                                             Xpfeh_gaia_dr3)

from target_selection.cartons import BaseCarton
# from target_selection.exceptions import TargetSelectionError


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class MWM_halo_distant_rrl_Carton(BaseCarton):
    """ 5.1.30. mwm_halo_distant_rrl
Shorthand name: mwm_halo_distant_rrl
Existing carton code: replaces open fiber 2020 28c
Simplified Description of selection criteria:
all Gaia DR3 RRLs with BP < 18.8.
SELECT b.*, a.best_classification, a.pf, a.p1_o, a.zp_mag_g,
 a.average_rv, a.num_clean_epochs_g
FROM gaiadr3.vari_rrlyrae as a
INNER JOIN gaiadr3.gaia_source_lite as b ON a.source_id = b.source_id
WHERE b.phot_bp_mean_mag < 18.8
Return columns:
Metadata:
Priority: 3050 (may change after A/B tests)
Cadence: bright_1x1
Instrument: BOSS
can_offset = True
Lead contact: Alexander Ji
    """

    name = 'mwm_halo_distant_rrl'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 3050
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_dr3_vari_rrlyrae.best_classification,
                         Gaia_dr3_vari_rrlyrae.pf,
                         Gaia_dr3_vari_rrlyrae.p1_o,
                         Gaia_dr3_vari_rrlyrae.zp_mag_g,
                         Gaia_dr3_vari_rrlyrae.average_rv,
                         Gaia_dr3_vari_rrlyrae.num_clean_epochs_g)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Gaia_dr3_vari_rrlyrae,
                       on=(Gaia_DR3.source_id == Gaia_dr3_vari_rrlyrae.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 18.8,))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_halo_distant_rrl_dark_Carton(BaseCarton):
    """  5.1.31. mwm_halo_distant_rrl_dark
Shorthand name: mwm_halo_distant_rrl_dark
Existing carton code:
Simplified Description of selection criteria:
all Gaia DR3 RRLs with BP < 18.8.
SELECT b.*, a.best_classification, a.pf, a.p1_o, a.zp_mag_g,
a.average_rv, a.num_clean_epochs_g
FROM gaiadr3.vari_rrlyrae as a
INNER JOIN gaiadr3.gaia_source_lite as b ON a.source_id = b.source_id
WHERE b.phot_bp_mean_mag < 18.8
Return columns:
Metadata:
Priority: 3049 (may change after A/B tests)
Cadence: dark_1x2
Instrument: BOSS
can_offset = True
Lead contact: Alexander Ji
    """

    name = 'mwm_halo_distant_rrl_dark'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'dark_1x2'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 3049
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_dr3_vari_rrlyrae.best_classification,
                         Gaia_dr3_vari_rrlyrae.pf,
                         Gaia_dr3_vari_rrlyrae.p1_o,
                         Gaia_dr3_vari_rrlyrae.zp_mag_g,
                         Gaia_dr3_vari_rrlyrae.average_rv,
                         Gaia_dr3_vari_rrlyrae.num_clean_epochs_g)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Gaia_dr3_vari_rrlyrae,
                       on=(Gaia_DR3.source_id == Gaia_dr3_vari_rrlyrae.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 18.8,))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_halo_mp_xp_Carton(BaseCarton):
    """  5.1.32.  mwm_halo_mp_xp
Shorthand name: mwm_halo_mp_xp
Existing carton code: N/A
Simplified Description of selection criteria:
XG Boost on XP spectra + WISE photometry to determine red giant metallicities.
Link to paper: https://arxiv.org/abs/2302.02611 and
Zenodo: https://doi.org/10.5281/zenodo.7599789
Use catalogdb.xpfeh_gaia_dr3.

Cut on BP<17, teff_xgboost < 5500, logg_xgboost < 4, mh_xgboost < -1.0
w1mpro (from AllWise) absolute magnitude.

Specifically:
(Note that below there are two cuts for M_W1
with teff_xgboost on right hand side)

BP < 17
logg_xgboost < 4.0
teff_xgboost < 5500
mh_xgboost < -1.0
M_W1 > -0.3 - 0.006 * (5500 - teff_xgboost)
M_W1 > -0.01 * (5300 - teff_xgboost)

where M_W1 = w1mpro + 5 log10(parallax/100)
(note: solve these equations so that this is a cut on parallax).

Then there are three levels of priority based on selecting
mh_xgboost <= -2.0, -2.0 < mh_xgboost <= -1.5, -1.5 < mh_xgboost <= -1.0
at three different priorities.

Return columns:
Metadata:
Priority: three levels
2100 if mh_xgboost <= -2.0 (TBA pending A/B test)
2970 if -2.0 < mh_xgboost <= -1.5
6090 if -1.5 < mh_xgboost <= -1.0
Cadence: bright_1x1
Instrument: BOSS for G>13, APOGEE for G<13
can_offset = True
Lead contact: Alexander Ji, Rene Andrae
    """

    name = 'mwm_halo_mp_xp'
    category = 'science'
    instrument = None  # instrument set in post_process()
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = None  # priority set in post_process()
    can_offset = True

    def build_query(self, version_id, query_region=None):

        # We use the lower case variable name m_w1 so that the
        # corresponding column name is lower case.
        m_w1 = AllWise.w1mpro + 5 * peewee.fn.log10(Gaia_DR3.parallax / 100)

        # w1mpro is of PostgreSQL numeric(5,3) type
        # Hence, below we cast it to real so that astropy can handle it.
        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_g_mean_mag,
                         Gaia_DR3.parallax,
                         Xpfeh_gaia_dr3.logg_xgboost,
                         Xpfeh_gaia_dr3.teff_xgboost,
                         Xpfeh_gaia_dr3.mh_xgboost,
                         AllWise.w1mpro.cast('real'),
                         m_w1.alias('m_w1'))
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Xpfeh_gaia_dr3,
                       on=(Gaia_DR3.source_id == Xpfeh_gaia_dr3.source_id))
                 .switch(CatalogToGaia_DR3)
                 .join(CatalogToAllWise,
                       on=(CatalogToGaia_DR3.catalogid == CatalogToAllWise.catalogid))
                 .join(AllWise,
                       on=(CatalogToAllWise.target_id == AllWise.cntr))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        CatalogToAllWise.version_id == version_id,
                        CatalogToAllWise.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 17,
                        Gaia_DR3.parallax > 0,
                        Xpfeh_gaia_dr3.logg_xgboost < 4.0,
                        Xpfeh_gaia_dr3.teff_xgboost < 5500,
                        Xpfeh_gaia_dr3.mh_xgboost < -1.0,
                        m_w1 > -0.3 - 0.006 * (5500 - Xpfeh_gaia_dr3.teff_xgboost),
                        m_w1 > -0.01 * (5300 - Xpfeh_gaia_dr3.teff_xgboost)))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """
        Priority: three levels
        2100 if mh_xgboost <= -2.0 (TBA pending A/B test)
        2970 if -2.0 < mh_xgboost <= -1.5
        6090 if -1.5 < mh_xgboost <= -1.0

        Instrument: BOSS for G>13, APOGEE for G<13
        """

        cursor = self.database.execute_sql(
            "select catalogid, mh_xgboost, phot_g_mean_mag from " +
            " sandbox.temp_mwm_halo_mp_xp ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_mh_xgboost = output[i][1]
            current_phot_g_mean_mag = output[i][2]

            # From query in build_query(), we have mh_xgboost < -1.0
            if (current_mh_xgboost <= -2.0):
                current_priority = 2100
            elif ((current_mh_xgboost > -2.0) and
                  (current_mh_xgboost <= -1.5)):
                current_priority = 2970
            else:
                current_priority = 6090

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_mp_xp " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if (current_phot_g_mean_mag < 13):
                current_instrument = 'APOGEE'
            else:
                current_instrument = 'BOSS'

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_mp_xp " +
                    " set instrument = '" + current_instrument + "'"
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_mp_xp_dark_Carton(BaseCarton):
    """   5.1.33. mwm_halo_mp_xp_dark
Shorthand name: mwm_halo_mp_xp_dark
Existing carton code: N/A dark time cadence for mwm_halo_mp_xp_dark
Simplified Description of selection criteria:
XG Boost on XP spectra + WISE photometry to determine red giant metallicities.
Link to paper: https://arxiv.org/abs/2302.02611 and
Zenodo: https://doi.org/10.5281/zenodo.7599789
Use catalogdb.xpfeh_gaia_dr3.

Cut on BP<17, teff_xgboost < 5500, logg_xgboost < 4, mh_xgboost < -1.0
w1mpro (from AllWise) absolute magnitude.

Specifically:
(Note that below there are two cuts for M_W1
with teff_xgboost on right hand side)

BP < 17
logg_xgboost < 4.0
teff_xgboost < 5500
mh_xgboost < -1.0
M_W1 > -0.3 - 0.006 * (5500 - teff_xgboost)
M_W1 > -0.01 * (5300 - teff_xgboost)

where M_W1 = w1mpro + 5 log10(parallax/100)
(note: solve these equations so that this is a cut on parallax).

Then there are three levels of priority based on selecting
mh_xgboost <= -2.0, -2.0 < mh_xgboost <= -1.5, -1.5 < mh_xgboost <= -1.0
at three different priorities.
Return columns:
Metadata:
Priority:
2099 if mh_xgboost <= -2.0 (TBA pending A/B test)
2969 if -2.0 < mh_xgboost <= -1.5
6090 if -1.5 < mh_xgboost <= -1.0
Cadence: dark_1x2
Instrument: BOSS for G>13, APOGEE for G<13
can_offset = True
Lead contact: Alexander Ji, Rene Andrae
    """

    name = 'mwm_halo_mp_xp_dark'
    category = 'science'
    instrument = None  # instrument set in post_process()
    cadence = 'dark_1x2'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = None  # priority set in post_process()
    can_offset = True

    def build_query(self, version_id, query_region=None):

        # We use the lower case variable name m_w1 so that the
        # corresponding column name is lower case.
        m_w1 = AllWise.w1mpro + 5 * peewee.fn.log10(Gaia_DR3.parallax / 100)

        # w1mpro is of PostgreSQL numeric(5,3) type
        # Hence, below we cast it to real so that astropy can handle it.
        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_g_mean_mag,
                         Gaia_DR3.parallax,
                         Xpfeh_gaia_dr3.logg_xgboost,
                         Xpfeh_gaia_dr3.teff_xgboost,
                         Xpfeh_gaia_dr3.mh_xgboost,
                         AllWise.w1mpro.cast('real'),
                         m_w1.alias('m_w1'))
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Xpfeh_gaia_dr3,
                       on=(Gaia_DR3.source_id == Xpfeh_gaia_dr3.source_id))
                 .switch(CatalogToGaia_DR3)
                 .join(CatalogToAllWise,
                       on=(CatalogToGaia_DR3.catalogid == CatalogToAllWise.catalogid))
                 .join(AllWise,
                       on=(CatalogToAllWise.target_id == AllWise.cntr))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        CatalogToAllWise.version_id == version_id,
                        CatalogToAllWise.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 17,
                        Gaia_DR3.parallax > 0,
                        Xpfeh_gaia_dr3.logg_xgboost < 4.0,
                        Xpfeh_gaia_dr3.teff_xgboost < 5500,
                        Xpfeh_gaia_dr3.mh_xgboost < -1.0,
                        m_w1 > -0.3 - 0.006 * (5500 - Xpfeh_gaia_dr3.teff_xgboost),
                        m_w1 > -0.01 * (5300 - Xpfeh_gaia_dr3.teff_xgboost)))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """
        Priority: three levels
        2099 if mh_xgboost <= -2.0 (TBA pending A/B test)
        2969 if -2.0 < mh_xgboost <= -1.5
        6090 if -1.5 < mh_xgboost <= -1.0

        Instrument: BOSS for G>13, APOGEE for G<13
        """

        cursor = self.database.execute_sql(
            "select catalogid, mh_xgboost, phot_g_mean_mag from " +
            " sandbox.temp_mwm_halo_mp_xp_dark ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_mh_xgboost = output[i][1]
            current_phot_g_mean_mag = output[i][2]

            # From query in build_query(), we have mh_xgboost < -1.0
            if (current_mh_xgboost <= -2.0):
                current_priority = 2099
            elif ((current_mh_xgboost > -2.0) and
                  (current_mh_xgboost <= -1.5)):
                current_priority = 2969
            else:
                current_priority = 6090

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_mp_xp_dark " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if (current_phot_g_mean_mag < 13):
                current_instrument = 'APOGEE'
            else:
                current_instrument = 'BOSS'

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_mp_xp_dark " +
                    " set instrument = '" + current_instrument + "'"
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_Base_Carton(BaseCarton):
    """
    This base carton contains the first part of the below query.
    which is common to all the mmw_halo_local* cartons
    5.1.36.  mwm_halo_local
Shorthand name: mwm_halo_local

Simplified Description of selection criteria: Gaia DR3 vtan > 150 km/s,
parallax_over_error > 5, parallax > 0.5, RP < 20, RUWE < 1.4,
with priority levels based on vtan, absolute G mag (prioritizing MSTO stars),
and parallax_over_error.
Return columns:
Metadata:

mwm_halo_local_high_apogee_single
mwm_halo_local_high_boss_single
mwm_halo_local_low_boss_single
mwm_halo_local_low_apogee_single
mwm_halo_local_high_apogee
mwm_halo_local_high_boss
mwm_halo_local_low_boss
mwm_halo_local_low_apogee
Existing carton code: mwm_halo_local ← this forms the base catalog.
The individual cartons below just change priority and cadence and instrument.
Note: The mwm_halo_local_single and mwm_halo_local cartons differ only in cadence
(bright_1x1 vs dark_1x2) and the numerical value of the priority.
The structure of the if statements for the priority is the same for both.
For mwm_halo_local, the priority is reduced by 1 (i.e. a bit higher priority)
in all cases.
Simplified Description of selection criteria:
Gaia DR3 vtan > 150 km/s, parallax_over_error > 5, parallax > 0.5, RP < 20,
RUWE < 1.4, with priority levels based on vtan,
absolute G mag (prioritizing MSTO stars), and parallax_over_error.

Details here: Development of Pseudocode for v1.0 cartons#LocalHalo
(note: consolidated the 11 priorities into 2)

mwm_halo_local_high_apogee_single
vtan >= 150 and parallax_over_error >=10 OR
vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
G < 13
PRIORITY: 2981
CADENCE: bright_1x1
INSTRUMENT: APOGEE
can_offset = True

mwm_halo_local_high_boss_single
vtan >= 150 and parallax_over_error >=10 OR
vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
G >= 13
PRIORITY: 2981
CADENCE: bright_1x1
INSTRUMENT: BOSS
can_offset = True

mwm_halo_local_low_apogee_single
vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
G < 13
PRIORITY: 6501
CADENCE: bright_1x1
INSTRUMENT: APOGEE
can_offset = True

mwm_halo_local_low_boss_single
vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
G >= 13
PRIORITY: 6501
CADENCE: bright_1x1
INSTRUMENT: BOSS
can_offset = True

mwm_halo_local_high_apogee
vtan >= 150 and parallax_over_error >=10 OR
vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
G < 13
PRIORITY: 2980
CADENCE: dark_flexible_2x1
INSTRUMENT: APOGEE
can_offset = True

mwm_halo_local_high_boss
vtan >= 150 and parallax_over_error >=10 OR
vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
G >= 13
PRIORITY: 2980
CADENCE: dark_flexible_2x1
INSTRUMENT: BOSS
can_offset = True

mwm_halo_local_low_apogee
vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
G < 13
PRIORITY: 6501
CADENCE: dark_flexible_2x1
INSTRUMENT: APOGEE
can_offset = True

mwm_halo_local_low_apogee
vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
G < 13
PRIORITY: 6500
CADENCE: dark_flexible_2x1
INSTRUMENT: APOGEE
can_offset = True

mwm_halo_local_low_boss
vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
G >= 13
PRIORITY: 6500
CADENCE: dark_flexible_2x1
INSTRUMENT: BOSS
can_offset = True

Lead contact: Alexander Ji, Sebastien Lepine

Details:
Target Selection Algorithm, purely on Gaia DR3 source catalog.
The overall catalog is vtan > 150 km/s, parallax_over_error > 5, within 2 kpc.
We increase priority for vtan > 200 km/s,
main sequence turnoff (MSTO) stars (defined by 3 < M_G < 5),
and higher parallax_over_error

The base catalog is something like this:
SELECT source_id,ra,dec,phot_g_mean_mag,phot_rp_mean_mag,parallax,
parallax_over_error,pmra,pmdec,RUWE,
4.74 * sqrt(pmra*pmra + pmdec*pmdec)/parallax as vtan,
phot_g_mean_mag - (10 - 5*log10(parallax)) as M_G,
FROM gaiadr3.gaia_source
WHERE parallax > 0.5 and parallax_over_error > 5.
and
phot_rp_mean_mag < 20. and RUWE < 1.4
and
vtan > 150


    """

    def build_query(self, version_id, query_region=None):

        # Below we must use power(Gaia_DR3.pmra, 2) instead of
        # Gaia_DR3.pmra*Gaia_DR3.pmra. This is because
        # PostgreSQL does not allow specifying the same column twice.
        vtan = 4.74 * peewee.fn.sqrt(peewee.fn.power(Gaia_DR3.pmra, 2) +
                                     peewee.fn.power(Gaia_DR3.pmdec, 2)) / Gaia_DR3.parallax
        # We compute m_g in in post_process() instead of below since
        # PostgreSQL does not allow
        # specifying Gaia_DR3.phot_g_mean_mag in m_g and in select clause.
        # m_g = Gaia_DR3.phot_g_mean_mag - (10 - 5 * peewee.fn.log10(Gaia_DR3.parallax))

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_g_mean_mag,
                         Gaia_DR3.phot_rp_mean_mag,
                         Gaia_DR3.parallax,
                         Gaia_DR3.parallax_over_error,
                         Gaia_DR3.pmra,
                         Gaia_DR3.pmdec,
                         Gaia_DR3.ruwe,
                         vtan.alias('vtan'))  # use alias since vtan is an expression
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.parallax > 0.5,
                        Gaia_DR3.parallax_over_error > 5,
                        Gaia_DR3.phot_rp_mean_mag < 20,
                        Gaia_DR3.phot_g_mean_mag.is_null(False),
                        Gaia_DR3.ruwe < 1.4,
                        vtan > 150))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_halo_local_high_apogee_single_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_high_apogee_single
    vtan >= 150 and parallax_over_error >=10 OR
    vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
    G < 13
    PRIORITY: 2981
    CADENCE: bright_1x1
    INSTRUMENT: APOGEE
    can_offset = True
    """
    name = 'mwm_halo_local_high_apogee_single'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 2981
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag < 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_high_apogee_single " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_high_apogee_single ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((is_m_g_3to5 and is_parallax_over_error_5to10) or
               (current_parallax_over_error >= 10)):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_high_apogee_single " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_high_boss_single_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_high_boss_single
    vtan >= 150 and parallax_over_error >=10 OR
    vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
    G >= 13
    PRIORITY: 2981
    CADENCE: bright_1x1
    INSTRUMENT: BOSS
    can_offset = True
    """
    name = 'mwm_halo_local_high_boss_single'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 2981
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag >= 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_high_boss_single " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_high_boss_single ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((is_m_g_3to5 and is_parallax_over_error_5to10) or
               (current_parallax_over_error >= 10)):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_high_boss_single " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_low_apogee_single_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_low_apogee_single
    vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
    G < 13
    PRIORITY: 6501
    CADENCE: bright_1x1
    INSTRUMENT: APOGEE
    can_offset = True
    """
    name = 'mwm_halo_local_low_apogee_single'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 6501
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag < 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_low_apogee_single " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_low_apogee_single ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((not is_m_g_3to5) and is_parallax_over_error_5to10):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_low_apogee_single " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_low_boss_single_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_low_boss_single
    vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
    G >= 13
    PRIORITY: 6501
    CADENCE: bright_1x1
    INSTRUMENT: BOSS
    can_offset = True
    """
    name = 'mwm_halo_local_low_boss_single'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 6501
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag >= 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_low_boss_single " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_low_boss_single ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((not is_m_g_3to5) and is_parallax_over_error_5to10):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_low_boss_single " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_high_apogee_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_high_apogee
    vtan >= 150 and parallax_over_error >=10 OR
    vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
    G < 13
    PRIORITY: 2980
    CADENCE: dark_flexible_2x1
    INSTRUMENT: APOGEE
    can_offset = True
    """
    name = 'mwm_halo_local_high_apogee'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'dark_flexible_2x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 2980
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag < 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_high_apogee " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_high_apogee ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((is_m_g_3to5 and is_parallax_over_error_5to10) or
               (current_parallax_over_error >= 10)):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_high_apogee " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_high_boss_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_high_boss
    vtan >= 150 and parallax_over_error >=10 OR
    vtan >=150 and 3 < M_G < 5 and 10 > parallax_over_error > 5
    G >= 13
    PRIORITY: 2980
    CADENCE: dark_flexible_2x1
    INSTRUMENT: BOSS
    can_offset = True
    """
    name = 'mwm_halo_local_high_boss'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'dark_flexible_2x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 2980
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag >= 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_high_boss " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_high_boss ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((is_m_g_3to5 and is_parallax_over_error_5to10) or
               (current_parallax_over_error >= 10)):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_high_boss " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_low_apogee_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_low_apogee
    vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
    G < 13
    PRIORITY: 6501
    CADENCE: dark_flexible_2x1
    INSTRUMENT: APOGEE
    can_offset = True
    """
    name = 'mwm_halo_local_low_apogee'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'dark_flexible_2x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 6501
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag < 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_low_apogee " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_low_apogee ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((not is_m_g_3to5) and is_parallax_over_error_5to10):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_low_apogee " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_low_boss_Carton(MWM_halo_local_Base_Carton):
    """
    mwm_halo_local_low_boss
    vtan >=150 and not (3 < M_G < 5) and 10 > parallax_over_error > 5
    G >= 13
    PRIORITY: 6500
    CADENCE: dark_flexible_2x1
    INSTRUMENT: BOSS
    can_offset = True
    """
    name = 'mwm_halo_local_low_boss'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'dark_flexible_2x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 6500
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3.phot_g_mean_mag >= 13)
        return query

    def post_process(self, model):

        cursor = self.database.execute_sql(
            "update sandbox.temp_mwm_halo_local_low_boss " +
            "set selected = false;")

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_low_boss ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            # current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            # The 'where' clause in the query in build_query() ensures
            # that current_phot_g_mean_mag and current_parallax are not null.
            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            is_m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            is_parallax_over_error_5to10 = ((5 < current_parallax_over_error) and
                                            (current_parallax_over_error < 10))

            # Below we do not check for vtan > 150 since
            # the query in base class build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in base class build_query() has parallax_over_error > 5

            if ((not is_m_g_3to5) and is_parallax_over_error_5to10):
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_low_boss " +
                    " set selected = true " +
                    " where catalogid = " + str(current_catalogid) + ";")
