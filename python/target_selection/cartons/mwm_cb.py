#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_cb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (CataclysmicVariables, Catalog,
                                             CatalogToGUVCat, CatalogToTIC_v8,
                                             Gaia_DR2,
                                             GeometricDistances_Gaia_DR2,
                                             GUVCat, TIC_v8)

from target_selection.cartons import BaseCarton


class MWM_CB_300_Carton(BaseCarton):
    """MWM Compact Binaries 300pc.

    Definition: Cross-match Gaia & Bailer-Jones distances by source_id,
    cross-match with GALEX including proper motion corrections. Define a
    single linear relation in absolute FUV magnitude vs FUV - NUV,
    select all objects with distances less than 300pc.

    SQL:
        (FUVmag - 5 * log10(r_est/10)) < 14 * (FUVmag-NUVmag) - 46 & r_est <300

    Cadence: This target sample will be split into
        (1) H<11 to be observed with APOGEE.
        (2) phot_g_mean_mag < 16 to be observed with BOSS in bright time.
        (3) phot_g_mean_mag > 16 to be observed with BOSS in dark time.

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):

        GD = GeometricDistances_Gaia_DR2

        FUV = GUVCat.fuv_mag
        NUV = GUVCat.nuv_mag

        FUV_abs = FUV - 5 * peewee.fn.log(GD.r_est / 10)

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         GD.source_id,
                         GD.r_est,
                         FUV,
                         NUV,
                         TIC_v8.hmag.alias('h'),
                         Gaia_DR2.phot_g_mean_mag.alias('gaia_g'))
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(GD)
                 .join_from(CatalogToTIC_v8, CatalogToGUVCat,
                            on=(CatalogToGUVCat.catalogid == CatalogToTIC_v8.catalogid))
                 .join(GUVCat)
                 .where(GD.r_est < 300,
                        FUV_abs < 14 * (FUV - NUV) - 46)
                 .where(CatalogToGUVCat.version_id == version_id,
                        CatalogToGUVCat.best >> True,
                        CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_300_APOGEE_Carton(MWM_CB_300_Carton):
    """MWM CB 300pc targets to be observed with APOGEE."""

    name = 'mwm_cb_300pc_apogee'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 1400

    def build_query(self, version_id, **kwargs):

        query = super().build_query(version_id, **kwargs)
        query.where(TIC_v8.hmag < 11)

        return query


class MWM_CB_300_BOSS_Carton(MWM_CB_300_Carton):
    """MWM CB 300pc targets to be observed with BOSS."""

    name = 'mwm_cb_300pc_boss'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def post_process(self, model, **kwargs):

        # G > 16 => cadence = dark_2x1
        model.update(cadence='dark_2x1').where(model.gaia_g > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence='bright_2x1').where(model.gaia_g < 16).execute()

        return model


class MWM_CB_Gaia_Galex_Carton(BaseCarton):
    """MWM Compact Binaries 300pc.

    Definition: Cross-match GALEX & Gaia, taking into account proper motions.

    SQL:
        parallax / parallax_error > 3 & gaiag < 20 &
        fuvmag + 5 * log10(parallax / 1000) + 5 > 1.5 + 1.28 * (fuvmag - gaiag)

    Cadence: This target sample will be split into
        (1) H<11 to be observed with APOGEE.
        (2) phot_g_mean_mag < 16 to be observed with BOSS in bright time.
        (3) phot_g_mean_mag > 16 to be observed with BOSS in dark time.

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):

        FUV = GUVCat.fuv_mag
        gaiag = Gaia_DR2.phot_g_mean_mag

        FUV_abs = FUV + 5 * peewee.fn.log(Gaia_DR2.parallax / 1000.) + 5

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id.alias('gaia_source_id'),
                         Gaia_DR2.parallax,
                         Gaia_DR2.parallax_error,
                         FUV,
                         gaiag.alias('phot_g_mean_mag'),
                         TIC_v8.hmag.alias('h'))
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join_from(CatalogToTIC_v8, CatalogToGUVCat,
                            on=(CatalogToGUVCat.catalogid == CatalogToTIC_v8.catalogid))
                 .join(GUVCat)
                 .where(Gaia_DR2.parallax / Gaia_DR2.parallax_error > 3,
                        gaiag < 20,
                        FUV > -999.,
                        FUV_abs > (1.5 + 1.28 * (FUV - gaiag)))
                 .where(CatalogToGUVCat.version_id == version_id,
                        CatalogToGUVCat.best >> True,
                        CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_Gaia_Galex_APOGEE_Carton(MWM_CB_Gaia_Galex_Carton):
    """MWM CB GaiaGalex to be observed with APOGEE."""

    name = 'mwm_cb_gaiagalex_apogee'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 1400

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region=query_region)

        return query.where(TIC_v8.hmag < 11)


class MWM_CB_Gaia_Galex_BOSS_Carton(MWM_CB_Gaia_Galex_Carton):
    """MWM CB GaiaGalex to be observed with BOSS."""

    name = 'mwm_cb_gaiagalex_boss'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def post_process(self, model, **kwargs):

        # G > 16 => cadence = dark_2x1
        model.update(cadence='dark_2x1').where(model.phot_g_mean_mag > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence='bright_2x1').where(model.phot_g_mean_mag < 16).execute()

        return model


class MWM_CB_CV_Candidates_Carton(BaseCarton):
    """MWM Compact Binaries 300pc.

    Definition:
        List of cataclysmic variables (+candidates)
        compiled by the AAVSO.

    SQL:
        Select all records from table cataclysmic_variables.

    Cadence: This target sample will be split into
        (1) H<11 to be observed with APOGEE.
        (2) phot_g_mean_mag < 16 to be observed with BOSS in bright time.
        (3) phot_g_mean_mag > 16 to be observed with BOSS in dark time.

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         CataclysmicVariables.source_id,
                         CataclysmicVariables.phot_g_mean_mag,
                         TIC_v8.hmag.alias('h'))
                 .join(TIC_v8)
                 .join(CataclysmicVariables,
                       on=(CataclysmicVariables.source_id == TIC_v8.gaia_int))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_CV_Candidates_APOGEE_Carton(MWM_CB_CV_Candidates_Carton):
    """MWM CB CV to be observed with APOGEE."""

    name = 'mwm_cb_cvcandidates_apogee'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'APOGEE'
    cadence = 'bright_1x1'
    priority = 1400

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region=query_region)

        return query.where(TIC_v8.hmag < 11)


class MWM_CB_CV_Candidates_BOSS_Carton(MWM_CB_CV_Candidates_Carton):
    """MWM CB CV to be observed with BOSS."""

    name = 'mwm_cb_cvcandidates_boss'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def post_process(self, model, **kwargs):

        # G > 16 => cadence = dark_2x1
        model.update(cadence='dark_2x1').where(model.phot_g_mean_mag > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence='bright_2x1').where(model.phot_g_mean_mag < 16).execute()

        return model
