#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_cb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (HECATE_1_1, MILLIQUAS_7_7,
                                             XMM_OM_SUSS_5_0, BailerJonesEDR3,
                                             Catalog, CatalogToGaia_DR3,
                                             CatalogToMILLIQUAS_7_7,
                                             CatalogToXMM_OM_SUSS_5_0,
                                             Gaia_DR3, Galex_GR7_Gaia_DR3,
                                             GUVCat)

from target_selection.cartons import BaseCarton


class MWM_CB_GALEX_Mag(BaseCarton):
    """Magnitude-limited GALEX FUV / Gaia selection.

    Simplified Description of selection criteria:
        Uses Nicola Gentile Fusillo's Gaia/GALEX cross-match (which should be already ingested).
        Select objects that lie below the main-sequence in an FUV/optical absolute magnitude-color
        diagram. Remove AGN interlopers. Remove the LMC and SMC, where a lot of hot young stars
        mimic the properties of this target selection. For more details, see here.

    Pseudo-code:
        Use the following short-hands:
            pm = sqrt(pmra*pmra+pmdec*pmdec)
            pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm)+
                            (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            logpmdpm = log10(pm/pm_error)
            poe = parallax / parallax_error
        Then combine the following cuts with &&
            fuv_mag + 5*log10(parallax/1000) + 5 > 1.5 + 1.28*(fuv_mag - phot_g_mean_mag)
            fuv_mag + 5*log10(parallax/1000) + 5 > 2.5*(fuv_mag - phot_g_mean_mag) - 0.5
            fuv_magerr < 0.3
            !(logpmdpm < 0.301 && poe > -1.75*logpmdpm - 2.8 && poe <  1.75*logpmdpm + 2.8 ||
              (logpmdpm-0.301)*(logpmdpm-0.301) / (0.33*0.33) + (poe * poe) / (3.2* 3.2)<=1)
            !((sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8) ||
              sqrt(0.5*square(l-301) + 2*square(b+43.3) ) < 5)

    """

    name = 'mwm_cb_galex_mag'
    mapper = 'MWM'
    category = 'science'
    instrument = 'BOSS'
    cadence = None
    program = 'mwm_cb'
    priority = 1400
    can_offset = True

    def build_query(self, version_id, query_region=None):

        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt((G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm) +
                           (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm))
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        fuv_mag = GUVCat.fuv_mag
        Gm = G3.phot_g_mean_mag
        parallax = G3.parallax
        ll = G3.l
        bb = G3.b

        fuv_cut = (((fuv_mag + 5 * fn.log(parallax / 1000) + 5) > (1.5 + 1.28 * (fuv_mag - Gm))) &
                   ((fuv_mag + 5 * fn.log(parallax / 1000) + 5) > (2.5 * (fuv_mag - Gm) - 0.5)))

        astro_cut1 = ~((logpmdpm < 0.301) &
                       (poe > (-1.75 * logpmdpm - 2.8)) &
                       (poe < (1.75 * logpmdpm + 2.8)) |
                       (((logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33) +
                         (poe * poe) / (3.2 * 3.2)) <= 1))

        astro_cut2 = ~((fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8) |
                       (fn.sqrt(0.5 * fn.pow(ll - 301, 2) + 2 * fn.pow(bb + 43.3, 2)) < 5))

        priority = peewee.Case(None, ((Gm <= 16, 'bright_2x1'), (Gm > 16, 'dark_2x1')))

        query = (Gaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         G3.source_id,
                         ll,
                         bb,
                         parallax,
                         Gm,
                         fuv_mag,
                         priority.alias('priority'))
                 .join(Galex_GR7_Gaia_DR3,
                       on=(Gaia_DR3.source_id == Galex_GR7_Gaia_DR3.gaia_edr3_source_id))
                 .join(GUVCat,
                       on=(GUVCat.objid == Galex_GR7_Gaia_DR3.galex_objid))
                 .join_from(Gaia_DR3, CatalogToGaia_DR3)
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        fuv_cut,
                        G3.parallax > 0,
                        G3.parallax_error > 0,
                        GUVCat.fuv_magerr < 0.3,
                        astro_cut1,
                        astro_cut2)
                 .order_by(CatalogToGaia_DR3.catalogid, Galex_GR7_Gaia_DR3.galex_separation)
                 .distinct(CatalogToGaia_DR3.catalogid))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
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

        return query.where(TIC_v8.hmag < 11)


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

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    Definition: Cross-match GALEX & Gaia, taking into account proper motions.

    SQL:
        parallax / parallax_error > 3 & gaiag < 20 &
        fuvmag + 5 * log10(parallax / 1000) + 5 > 1.5 + 1.28 * (fuvmag - gaiag)

    Cadence: This target sample will be split into
        (1) H<11 to be observed with APOGEE.
        (2) phot_g_mean_mag < 16 to be observed with BOSS in bright time.
        (3) phot_g_mean_mag > 16 to be observed with BOSS in dark time.

        In this carton, we use FUV > -999.
        where FUV = GUVCat.fuv_mag.
        In mwm_cb_uvex.py, we use fuv_mag > -100.
        This difference is for historical reasons.

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

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    Definition:
        List of cataclysmic variables (+candidates)
        compiled by the AAVSO.

    SQL:
        Select all records from table cataclysmic_variables.

    Cadence: This target sample will be split into
        (1) H<11 to be observed with APOGEE.
        (2) phot_g_mean_mag < 16 to be observed with BOSS in bright time.
        (3) phot_g_mean_mag > 16 to be observed with BOSS in dark time.

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
