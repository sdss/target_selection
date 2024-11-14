#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_halo.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    BestBrightest,
    Catalog,
    CatalogToAllWise,
    CatalogToTIC_v8,
    Gaia_DR2,
    SkyMapperGaia,
    TIC_v8,
)

from target_selection.cartons import BaseCarton


class MWM_Halo_Best_Brightest_Base_Carton(BaseCarton):
    """MWM Halo Best & Brightest.

    Definition:
        Photometrically selected metal-poor stars to (hopefully) be
        observed with spare fibers. Selected using "v1" and "v2"
        from Schlaufman + Casey 2014.

        A sample of metal-poor giants selected based on optical and
        infrared photometry.

        Select all stars from B&B catalog in two rounds. First priority
        with version=2, second priority by version=1.

        Original selection criteria:

            0.45 < J-H < 0.6
            W3 > 8
            -0.04 < W1-W2 < 0.04
            J - W2 > 0.5
            e^z / (1 + e^z) > 0.13 where z is defined in Schlaufman+Casey
            (version=1): J-W2 > 0.5 (B-V -0.8) + 0.6
            (version=2): B-V < 1.2

            (250k sources for version=2, 396k sources for version=1)

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToAllWise.select(
                CatalogToAllWise.catalogid,
                BestBrightest.designation,
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.parallax,
                Gaia_DR2.parallax_error,
                BestBrightest.version,
            )
            .join(BestBrightest, on=(BestBrightest.cntr == CatalogToAllWise.target_id))
            .switch(CatalogToAllWise)
            .join(CatalogToTIC_v8, on=(CatalogToAllWise.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8)
            .join(Gaia_DR2)
            .where(
                CatalogToAllWise.version_id == version_id,
                CatalogToAllWise.best >> True,
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
            )
        )

        if query_region:
            query = query.join_from(CatalogToAllWise, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


class MWM_Halo_Best_Brightest_BOSS_Carton(MWM_Halo_Best_Brightest_Base_Carton):
    """MWM Halo B&B targets for BOSS."""

    name = "mwm_halo_bb_boss"
    mapper = "MWM"
    category = "science"
    program = "mwm_filler"
    instrument = "BOSS"
    cadence = "bright_1x1"
    priority = None

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region=query_region)

        return query

    def post_process(self, model, **kwargs):
        version = model.version
        g = model.phot_g_mean_mag
        parallax = model.parallax
        parallax_error = model.parallax_error
        abs_g = g - 10 + 5 * peewee.fn.log(parallax)

        # BB1: Select (Best & Brightest catalog version = 2) AND (Gaia G < 16) AND
        #   ((G - 10 + 5*log10(parallax) < 5) OR (parallax < 5*parallax_error))
        # Priority 6050
        # Expect 66808 sources
        (
            model.update({model.priority: 6050})
            .where(version == 2)
            .where(g < 16)
            .where(parallax > 0)
            .where((abs_g < 5) | (parallax < 5 * parallax_error))
        ).execute()

        # BB2: Select (Best & Brightest catalog version = 1) AND (Gaia G < 16) AND
        #   ((G - 10 + 5*log10(parallax) < 5) OR (parallax < 5*parallax_error))
        # Priority 6060
        # Expect 114667 sources
        (
            model.update({model.priority: 6060})
            .where(version == 1)
            .where(g < 16)
            .where(parallax > 0)
            .where((abs_g < 5) | (parallax < 5 * parallax_error))
        ).execute()

        # BB3: Select (Best & Brightest catalog version = 2) AND (Gaia G >= 16) AND
        #   ((G - 10 + 5*log10(parallax) < 5) OR (parallax < 5*parallax_error))
        # Priority 6065
        # Expect 25894 sources
        (
            model.update({model.priority: 6065})
            .where(version == 2)
            .where(g >= 16)
            .where(parallax > 0)
            .where((abs_g < 5) | (parallax < 5 * parallax_error))
        ).execute()

        # BB4: Select (Best & Brightest catalog version = 1) AND (Gaia G >= 16) AND
        #   ((G - 10 + 5*log10(parallax) < 5) OR (parallax < 5*parallax_error))
        # Priority 6075
        # Expect 25894 sources
        (
            model.update({model.priority: 6075})
            .where(version == 1)
            .where(g >= 16)
            .where(parallax > 0)
            .where((abs_g < 5) | (parallax < 5 * parallax_error))
        ).execute()

        # BB5: (Best & Brightest catalog version = 2) AND (Gaia G < 16) AND
        #   not in mwm_halo_bb1
        # Priority 6080
        # Expect 108439 sources
        (
            model.update({model.priority: 6080})
            .where(version == 2)
            .where(g < 16)
            .where((parallax < 0) | ~((abs_g < 5) | (parallax < 5 * parallax_error)))
        ).execute()

        # BB6: (Best & Brightest catalog version = 1) AND (Gaia G < 16) AND
        #   not in mwm_halo_bb2
        # Priority 6085
        # Expect 195816 sources
        (
            model.update({model.priority: 6085})
            .where(version == 1)
            .where(g < 16)
            .where((parallax < 0) | ~((abs_g < 5) | (parallax < 5 * parallax_error)))
        ).execute()

        # Assign a base priority to the remaining targets
        (model.update({model.priority: 6090}).where(model.priority.is_null())).execute()


class MWM_Halo_SkyMapper_Base_Carton(BaseCarton):
    """MWM Halo SkyMapper.

    Definition: Photometrically selected metal-poor stars with [Fe/H] < -1.7
    (determined using SkyMapper Ca K photometry).

    This is a base carton. Actual cartons are implemented as subclasses for the
    different magnitude cuts.

    """

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToTIC_v8.select(
                CatalogToTIC_v8.catalogid, SkyMapperGaia.gaia_source_id, Gaia_DR2.phot_g_mean_mag
            )
            .join(TIC_v8)
            .join(SkyMapperGaia, on=(TIC_v8.gaia_int == SkyMapperGaia.gaia_source_id))
            .join(Gaia_DR2, on=(Gaia_DR2.source_id == SkyMapperGaia.gaia_source_id))
            .where(
                SkyMapperGaia.feh < -1.7,
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
            )
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


class MWM_Halo_SkyMapper_BOSS_Carton(MWM_Halo_SkyMapper_Base_Carton):
    """MWM Halo SkyMapper BOSS targets."""

    name = "mwm_halo_sm_boss"
    mapper = "MWM"
    category = "science"
    program = "mwm_filler"
    instrument = "BOSS"
    cadence = "bright_1x1"
    priority = None

    def build_query(self, version_id, **kwargs):
        query = super().build_query(version_id, **kwargs)

        return query

    def post_process(self, model, **kwargs):
        # Select (SkyMapper Fe/H < -1.7) AND (Gaia G < 16)
        # Priority 6055
        # Expect 49457 sources
        model.update(priority=6055).where(model.phot_g_mean_mag < 16).execute()

        # Select (SkyMapper Fe/H < -1.7) AND (Gaia G >= 16)
        # Priority 6070
        # Expect 88841 sources
        model.update(priority=6070).where(model.phot_g_mean_mag >= 16).execute()
