#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-22
# @Filename: mwm_halo_gaia_dr3.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             Gaia_DR3, Gaia_dr3_vari_rrlyrae,
                                             Xpfeh_gaia_dr3)

from target_selection.cartons import BaseCarton


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

Cut on BP<17, teff_xgboost < 5500, logg_xgboost < 4, W1 absolute magnitude.

Specifically:
BP < 17
logg_xgboost < 4.0
teff_xgboost < 5500
M_W1 > -0.3 - 0.006 * (5500 - teff_xgboost)
M_W1 > -0.01 * (5300 - teff_xgboost)

where M_W1 = W_1 + 5 log10(parallax/100)
(note: solve these equations so that this is a cut on parallax).

Then there are three levels of priority based on selecting
mh_xgboost <= -2.0, -2.0 < mh_xgboost <= -1.5, -1.5 < mh_xgboost <= -1.0
at three different priorities.

Return columns:
Metadata:
Priority: three levels
2100 if mh_xgboost <= -2.0 (TBA pending A/B test)
2970 if -2.0 < mh_xgboost <= -1.5
6090 if mh_xgboost > -1.5
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

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_g_mean_mag,
                         Xpfeh_gaia_dr3.logg_xgboost,
                         Xpfeh_gaia_dr3.teff_xgboost,
                         Xpfeh_gaia_dr3.mh_xgboost)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Xpfeh_gaia_dr3,
                       on=(Gaia_DR3.source_id == Xpfeh_gaia_dr3.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 17,
                        Xpfeh_gaia_dr3.logg_xgboost < 4.0,
                        Xpfeh_gaia_dr3.teff_xgboost < 5500))

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
        6090 if mh_xgboost > -1.5

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
                    " set priority = '" + current_priority + "'"
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

Cut on BP<17, teff_xgboost < 5500, logg_xgboost < 4, W1 absolute magnitude.

Specifically:
BP < 17
logg_xgboost < 4.0
teff_xgboost < 5500
M_W1 > -0.3 - 0.006 * (5500 - teff_xgboost)
M_W1 > -0.01 * (5300 - teff_xgboost)

where M_W1 = W_1 + 5 log10(parallax/100)
(note: solve these equations so that this is a cut on parallax).

Then there are three levels of priority based on selecting
mh_xgboost <= -2.0, -2.0 < mh_xgboost <= -1.5, -1.5 < mh_xgboost <= -1.0
at three different priorities.
Return columns:
Metadata:
Priority:
2099 if mh_xgboost <= -2.0 (TBA pending A/B test)
2969 if -2.0 < mh_xgboost <= -1.5
6090 if mh_xgboost > -1.5
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

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_g_mean_mag,
                         Xpfeh_gaia_dr3.logg_xgboost,
                         Xpfeh_gaia_dr3.teff_xgboost,
                         Xpfeh_gaia_dr3.mh_xgboost)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Xpfeh_gaia_dr3,
                       on=(Gaia_DR3.source_id == Xpfeh_gaia_dr3.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 17,
                        Xpfeh_gaia_dr3.logg_xgboost < 4.0,
                        Xpfeh_gaia_dr3.teff_xgboost < 5500))

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
        6090 if mh_xgboost > -1.5

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
                    " set priority = '" + current_priority + "'"
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


class MWM_halo_local_Carton(BaseCarton):
    """  5.1.36.  mwm_halo_local
Shorthand name: mwm_halo_local
Existing carton code: N/A
Simplified Description of selection criteria: Gaia DR3 vtan > 150 km/s,
parallax_over_error > 5, parallax > 0.5, RP < 20, RUWE < 1.4,
with priority levels based on vtan, absolute G mag (prioritizing MSTO stars),
and parallax_over_error.
Return columns:
Metadata:
Priority: 11 levels
vtan >= 200 and 3 < M_G < 5 and parallax_over_error >= 50: 2980

200 > vtan >= 150 and 3 < M_G < 5 and parallax_over_error >= 50: 2985

vtan >= 200 and not (3 < M_G < 5) and parallax_over_error >= 50: 2990

200 > vtan >= 150 and not (3 < M_G < 5) and parallax_over_error >= 50: 2995

vtan >= 200 and 3 < M_G < 5 and 50 > parallax_over_error >= 10: 3020

200 > vtan >= 150 and 3 < M_G < 5 and 50 > parallax_over_error >= 10: 3025

vtan >= 150 and 3 < M_G < 5 and 10 > parallax_over_error >= 5: 3030

vtan >= 200 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10: 3040

200 > vtan >= 150 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10: 3045

vtan >= 200 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5: 6091

200 > vtan >= 150 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5: 6092
Cadence: bright_1x1
Instrument: BOSS for G>13, APOGEE for G<13
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

    name = 'mwm_halo_local'
    category = 'science'
    instrument = None  # instrument set in post_process()
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = None  # priority set in post_process()
    can_offset = True

    def build_query(self, version_id, query_region=None):

        vtan = 4.74 * peewee.fn.sqrt(Gaia_DR3.pmra * Gaia_DR3.pmra +
                                     Gaia_DR3.pmdec * Gaia_DR3.pmdec) / Gaia_DR3.parallax
        m_g = Gaia_DR3.phot_g_mean_mag - (10 - 5 * peewee.fn.log10(Gaia_DR3.parallax))

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_g_mean_mag,
                         Gaia_DR3.phot_rp_mean_mag,
                         Gaia_DR3.parallax,
                         Gaia_DR3.paralax_over_error,
                         Gaia_DR3.pmra,
                         Gaia_DR3.pmdec,
                         Gaia_DR3.ruwe,
                         vtan,
                         m_g)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.parallax > 0.5,
                        Gaia_DR3.paralax_over_error > 5,
                        Gaia_DR3.phot_rp_mean_mag < 20,
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

    def post_process(self, model):
        """
Priority: 11 levels
vtan >= 200 and 3 < M_G < 5 and parallax_over_error >= 50: 2980
200 > vtan >= 150 and 3 < M_G < 5 and parallax_over_error >= 50: 2985
vtan >= 200 and not (3 < M_G < 5) and parallax_over_error >= 50: 2990
200 > vtan >= 150 and not (3 < M_G < 5) and parallax_over_error >= 50: 2995
vtan >= 200 and 3 < M_G < 5 and 50 > parallax_over_error >= 10: 3020
200 > vtan >= 150 and 3 < M_G < 5 and 50 > parallax_over_error >= 10: 3025
vtan >= 150 and 3 < M_G < 5 and 10 > parallax_over_error >= 5: 3030
vtan >= 200 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10: 3040
200 > vtan >= 150 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10: 3045
vtan >= 200 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5: 6091
200 > vtan >= 150 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5: 6092

        Instrument: BOSS for G>13, APOGEE for G<13
        """

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, m_g from " +
            " sandbox.temp_mwm_halo_local ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_vtan = output[i][1]
            # current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            # current_m_g = output[i][4]

            # TODO
            if (current_vtan <= -2.0):
                current_priority = 2100
            elif ((current_vtan > -2.0) and
                  (current_vtan <= -1.5)):
                current_priority = 2970
            else:
                current_priority = 6090

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local " +
                    " set priority = '" + current_priority + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            if (current_phot_g_mean_mag < 13):
                current_instrument = 'APOGEE'
            else:
                current_instrument = 'BOSS'

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local " +
                    " set instrument = '" + current_instrument + "'"
                    " where catalogid = " + str(current_catalogid) + ";")
