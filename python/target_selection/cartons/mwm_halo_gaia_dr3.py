#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-22
# @Filename: mwm_halo_gaia_dr3.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import math

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             Gaia_DR3, Gaia_dr3_vari_rrlyrae,
                                             Xpfeh_gaia_dr3)

from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError


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


class MWM_halo_local_Carton(BaseCarton):
    """  5.1.36.  mwm_halo_local
Shorthand name: mwm_halo_local
Existing carton code: N/A

Note: The mwm_halo_local and mwm_halo_local_dark  cartons differ
only in cadence (bright_1x1 vs dark_1x2) and
the numerical value of the priority.
The structure of the if statements for the priority is the same for both.
For mwm_halo_local_dark, the priority is reduced by 1 (i.e. a bit higher priority) in all cases.

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

Priority algorithm: (exact numbers are TBA, this is based on v0.5 priorities)
LH_MSTO1 priority 2980: vtan >= 200      and 3 < M_G < 5       and parallax_over_error >= 50 (N=14k)  # noqa: E501
LH_MSTO2 priority 2985: 200 > vtan >= 150 and 3 < M_G < 5       and parallax_over_error >= 50 (N=24k)  # noqa: E501
LH_ALL1  priority 2990: vtan >= 200      and not (3 < M_G < 5) and parallax_over_error >= 50 (N=28k)  # noqa: E501
LH_ALL2  priority 2995: 200 > vtan >= 150 and not (3 < M_G < 5) and parallax_over_error >= 50 (N=42k)  # noqa: E501
LH_MSTO3 priority 3020: vtan >= 200      and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=98k)  # noqa: E501
LH_MSTO4 priority 3025: 200 > vtan >= 150 and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=179k)  # noqa: E501
LH_MSTO5 priority 3030: vtan >= 150       and 3 < M_G < 5       and 10 > parallax_over_error >= 5 (N=12k)  # noqa: E501
LH_ALL3  priority 3040: vtan >= 200      and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=228k)  # noqa: E501
LH_ALL4  priority 3045: 200 > vtan >= 150 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=360k)  # noqa: E501
LH_ALL5  priority 6091: vtan >= 200      and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=274k)  # noqa: E501)
LH_ALL6  priority 6092: 200 > vtan >= 150 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=437k)  # noqa: E501

Note: LH_MSTO5 priority 3030 Covers both vtan>200 and vtan>150.

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

Priority algorithm: (exact numbers are TBA, this is based on v0.5 priorities)
LH_MSTO1 priority 2980: vtan >= 200      and 3 < M_G < 5       and parallax_over_error >= 50 (N=14k)  # noqa: E501
LH_MSTO2 priority 2985: 200 > vtan >= 150 and 3 < M_G < 5       and parallax_over_error >= 50 (N=24k)  # noqa: E501
LH_ALL1  priority 2990: vtan >= 200      and not (3 < M_G < 5) and parallax_over_error >= 50 (N=28k)  # noqa: E501
LH_ALL2  priority 2995: 200 > vtan >= 150 and not (3 < M_G < 5) and parallax_over_error >= 50 (N=42k)  # noqa: E501
LH_MSTO3 priority 3020: vtan >= 200      and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=98k)  # noqa: E501
LH_MSTO4 priority 3025: 200 > vtan >= 150 and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=179k)  # noqa: E501
LH_MSTO5 priority 3030: vtan >= 150       and 3 < M_G < 5       and 10 > parallax_over_error >= 5 (N=12k)  # noqa: E501
LH_ALL3  priority 3040: vtan >= 200      and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=228k)  # noqa: E501
LH_ALL4  priority 3045: 200 > vtan >= 150 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=360k)  # noqa: E501
LH_ALL5  priority 6091: vtan >= 200      and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=274k)  # noqa: E501)
LH_ALL6  priority 6092: 200 > vtan >= 150 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=437k)  # noqa: E501

Note: LH_MSTO5 priority 3030 Covers both vtan>200 and vtan>150.

        Instrument: BOSS for G>13, APOGEE for G<13
        """

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            parallax_over_error_10to50 = ((10 <= current_parallax_over_error) and
                                          (current_parallax_over_error < 50))

            current_priority = None
            # Below we do not check for vtan > 150 since
            # the query in build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in build_query() has parallax_over_error > 5
            #
            # The order below is the same order as in the comment above.
            if (m_g_3to5 and (current_parallax_over_error >= 50)):
                if (current_vtan >= 200):
                    current_priority = 2980
                else:
                    current_priority = 2985

            elif ((not m_g_3to5) and (current_parallax_over_error >= 50)):
                if (current_vtan >= 200):
                    current_priority = 2990
                else:
                    current_priority = 2995

            elif (m_g_3to5 and parallax_over_error_10to50):
                if (current_vtan >= 200):
                    current_priority = 3020
                else:
                    current_priority = 3025

            # See the case LH_MSTO5 in the above comments.
            # This is the odd one out since it sets the
            # same priority for both branches of if(current_vtan >= 200).
            # We could remove the inside if() but we are keeping it
            # so that the code is similar to the other elif sections.
            elif (m_g_3to5 and (10 > current_parallax_over_error)):
                if (current_vtan >= 200):
                    current_priority = 3030
                else:
                    current_priority = 3030

            elif ((not m_g_3to5) and parallax_over_error_10to50):
                if (current_vtan >= 200):
                    current_priority = 3040
                else:
                    current_priority = 3045

            elif ((not m_g_3to5) and (10 > current_parallax_over_error)):
                if (current_vtan >= 200):
                    current_priority = 6091
                else:
                    current_priority = 6092

            else:
                raise TargetSelectionError('error MWM_halo_local_Carton: ' +
                                           'post_process(): ' +
                                           ' no priority assigned for target: ' +
                                           'm_g=' + current_m_g + ', '
                                           'parallax_over_error=' +
                                           current_parallax_over_error + ',' +
                                           'vtan=' + current_vtan)

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            current_instrument = None
            if (current_phot_g_mean_mag < 13):
                current_instrument = 'APOGEE'
            else:
                current_instrument = 'BOSS'

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local " +
                    " set instrument = '" + current_instrument + "'"
                    " where catalogid = " + str(current_catalogid) + ";")


class MWM_halo_local_dark_Carton(BaseCarton):
    """  5.1.37.  mwm_halo_local_dark
Shorthand name: mwm_halo_local_dark
Existing carton code: N/A

Note: The mwm_halo_local and mwm_halo_local_dark  cartons differ
only in cadence (bright_1x1 vs dark_1x2) and
the numerical value of the priority.
The structure of the if statements for the priority is the same for both.
For mwm_halo_local_dark, the priority is reduced by 1 (i.e. a bit higher priority) in all cases.

Simplified Description of selection criteria: Gaia DR3 vtan > 150 km/s,
parallax_over_error > 5, parallax > 0.5, RP < 20, RUWE < 1.4,
with priority levels based on vtan, absolute G mag (prioritizing MSTO stars),
and parallax_over_error.
Return columns:
Metadata:

The priority numbers below are from mwm_halo_local.
The priority numbers for this carton (mwm_halo_local_dark) are one less
than the corresponding priority numbers for mwm_halo_local.

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

Cadence: dark_1x1
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

The priority numbers below are from mwm_halo_local.
The priority numbers for this carton (mwm_halo_local_dark) are one less
than the corresponding priority numbers for mwm_halo_local.

Priority algorithm: (exact numbers are TBA, this is based on v0.5 priorities)
LH_MSTO1 priority 2980: vtan >= 200      and 3 < M_G < 5       and parallax_over_error >= 50 (N=14k)  # noqa: E501
LH_MSTO2 priority 2985: 200 > vtan >= 150 and 3 < M_G < 5       and parallax_over_error >= 50 (N=24k)  # noqa: E501
LH_ALL1  priority 2990: vtan >= 200      and not (3 < M_G < 5) and parallax_over_error >= 50 (N=28k)  # noqa: E501
LH_ALL2  priority 2995: 200 > vtan >= 150 and not (3 < M_G < 5) and parallax_over_error >= 50 (N=42k)  # noqa: E501
LH_MSTO3 priority 3020: vtan >= 200      and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=98k)  # noqa: E501
LH_MSTO4 priority 3025: 200 > vtan >= 150 and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=179k)  # noqa: E501
LH_MSTO5 priority 3030: vtan >= 150       and 3 < M_G < 5       and 10 > parallax_over_error >= 5 (N=12k)  # noqa: E501
LH_ALL3  priority 3040: vtan >= 200      and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=228k)  # noqa: E501
LH_ALL4  priority 3045: 200 > vtan >= 150 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=360k)  # noqa: E501
LH_ALL5  priority 6091: vtan >= 200      and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=274k)  # noqa: E501)
LH_ALL6  priority 6092: 200 > vtan >= 150 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=437k)  # noqa: E501

Note: LH_MSTO5 priority 3030 Covers both vtan>200 and vtan>150.

    """

    name = 'mwm_halo_local_dark'
    category = 'science'
    instrument = None  # instrument set in post_process()
    cadence = 'dark_1x2'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = None  # priority set in post_process()
    can_offset = True

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
The priority numbers below are from mwm_halo_local.
The priority numbers for this carton (mwm_halo_local_dark) are one less
than the corresponding priority numbers for mwm_halo_local.

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

The priority numbers below are from mwm_halo_local.
The priority numbers for this carton (mwm_halo_local_dark) are one less
than the corresponding priority numbers for mwm_halo_local.

Priority algorithm: (exact numbers are TBA, this is based on v0.5 priorities)
LH_MSTO1 priority 2980: vtan >= 200      and 3 < M_G < 5       and parallax_over_error >= 50 (N=14k)  # noqa: E501
LH_MSTO2 priority 2985: 200 > vtan >= 150 and 3 < M_G < 5       and parallax_over_error >= 50 (N=24k)  # noqa: E501
LH_ALL1  priority 2990: vtan >= 200      and not (3 < M_G < 5) and parallax_over_error >= 50 (N=28k)  # noqa: E501
LH_ALL2  priority 2995: 200 > vtan >= 150 and not (3 < M_G < 5) and parallax_over_error >= 50 (N=42k)  # noqa: E501
LH_MSTO3 priority 3020: vtan >= 200      and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=98k)  # noqa: E501
LH_MSTO4 priority 3025: 200 > vtan >= 150 and 3 < M_G < 5       and 50 > parallax_over_error >= 10 (N=179k)  # noqa: E501
LH_MSTO5 priority 3030: vtan >= 150       and 3 < M_G < 5       and 10 > parallax_over_error >= 5 (N=12k)  # noqa: E501
LH_ALL3  priority 3040: vtan >= 200      and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=228k)  # noqa: E501
LH_ALL4  priority 3045: 200 > vtan >= 150 and not (3 < M_G < 5) and 50 > parallax_over_error >= 10 (N=360k)  # noqa: E501
LH_ALL5  priority 6091: vtan >= 200      and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=274k)  # noqa: E501)
LH_ALL6  priority 6092: 200 > vtan >= 150 and not (3 < M_G < 5) and 10 > parallax_over_error >= 5 (N=437k)  # noqa: E501

Note: LH_MSTO5 priority 3030 Covers both vtan>200 and vtan>150.

        Instrument: BOSS for G>13, APOGEE for G<13
        """

        cursor = self.database.execute_sql(
            "select catalogid, vtan, " +
            "parallax_over_error, phot_g_mean_mag, parallax from " +
            " sandbox.temp_mwm_halo_local_dark ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_vtan = output[i][1]
            current_parallax_over_error = output[i][2]
            current_phot_g_mean_mag = output[i][3]
            current_parallax = output[i][4]

            current_m_g = current_phot_g_mean_mag - (10 - 5 * math.log10(current_parallax))

            m_g_3to5 = ((3 < current_m_g) and (current_m_g < 5))
            parallax_over_error_10to50 = ((10 <= current_parallax_over_error) and
                                          (current_parallax_over_error < 50))

            current_priority = None
            # Below we do not check for vtan > 150 since
            # the query in build_query() has vtan > 150
            #
            # We do not check for parallax_over_error > 5 since
            # the query in build_query() has parallax_over_error > 5
            #
            # The order below is the same order as in the comment above.
            if (m_g_3to5 and (current_parallax_over_error >= 50)):
                if (current_vtan >= 200):
                    current_priority = 2979
                else:
                    current_priority = 2984

            elif ((not m_g_3to5) and (current_parallax_over_error >= 50)):
                if (current_vtan >= 200):
                    current_priority = 2989
                else:
                    current_priority = 2994

            elif (m_g_3to5 and parallax_over_error_10to50):
                if (current_vtan >= 200):
                    current_priority = 3019
                else:
                    current_priority = 3024

            # See the case LH_MSTO5 in the above comments.
            # This is the odd one out since it sets the
            # same priority for both branches of if(current_vtan >= 200).
            # We could remove the inside if() but we are keeping it
            # so that the code is similar to the other elif sections.
            elif (m_g_3to5 and (10 > current_parallax_over_error)):
                if (current_vtan >= 200):
                    current_priority = 3029
                else:
                    current_priority = 3029

            elif ((not m_g_3to5) and parallax_over_error_10to50):
                if (current_vtan >= 200):
                    current_priority = 3039
                else:
                    current_priority = 3044

            elif ((not m_g_3to5) and (10 > current_parallax_over_error)):
                if (current_vtan >= 200):
                    current_priority = 6090
                else:
                    current_priority = 6091

            else:
                raise TargetSelectionError('error MWM_halo_local_dark_Carton: ' +
                                           'post_process(): ' +
                                           ' no priority assigned for target: ' +
                                           'm_g=' + current_m_g + ', '
                                           'parallax_over_error=' +
                                           current_parallax_over_error + ',' +
                                           'vtan=' + current_vtan)

            if current_priority is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_dark " +
                    " set priority = '" + str(current_priority) + "'"
                    " where catalogid = " + str(current_catalogid) + ";")

            current_instrument = None
            if (current_phot_g_mean_mag < 13):
                current_instrument = 'APOGEE'
            else:
                current_instrument = 'BOSS'

            if current_instrument is not None:
                self.database.execute_sql(
                    " update sandbox.temp_mwm_halo_local_dark " +
                    " set instrument = '" + current_instrument + "'"
                    " where catalogid = " + str(current_catalogid) + ";")
