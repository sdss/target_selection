#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-05-22
# @Filename: mwm_cb_cvcandidates.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import CataclysmicVariables, Catalog, CatalogToTIC_v8, TIC_v8

from target_selection.cartons import BaseCarton


class MWM_CB_CV_Candidates_Carton(BaseCarton):
    """
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
        # The column CataclysmicVariables.source_id
        # corresponds to Gaia_DR2 source_id.
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid,
                CataclysmicVariables.source_id,
                CataclysmicVariables.phot_g_mean_mag,
                TIC_v8.hmag.alias("h"),
            )
            .join(TIC_v8)
            .join(CataclysmicVariables, on=(CataclysmicVariables.source_id == TIC_v8.gaia_int))
            .where(CatalogToTIC_v8.version_id == version_id, CatalogToTIC_v8.best >> True)
        )

        if query_region:
            query = query.join_from(CatalogToTIC_v8, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


class MWM_CB_CV_Candidates_APOGEE_Carton(MWM_CB_CV_Candidates_Carton):
    """MWM CB CV to be observed with APOGEE."""

    name = "mwm_cb_cvcandidates_apogee"
    mapper = "MWM"
    category = "science"
    program = "mwm_cb"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 1820
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region=query_region)

        return query.where(TIC_v8.hmag < 11)


class MWM_CB_CV_Candidates_BOSS_Carton(MWM_CB_CV_Candidates_Carton):
    """MWM CB CV to be observed with BOSS."""

    name = "mwm_cb_cvcandidates_boss"
    mapper = "MWM"
    category = "science"
    program = "mwm_cb"
    instrument = "BOSS"
    cadence = None  # cadence is set in post_process()
    priority = 1820
    can_offset = True

    def post_process(self, model, **kwargs):
        # G > 16 => cadence = dark_2x1
        model.update(cadence="dark_2x1").where(model.phot_g_mean_mag > 16).execute()

        # G < 16 => cadence = bright_2x1
        model.update(cadence="bright_2x1").where(model.phot_g_mean_mag < 16).execute()

        return model
