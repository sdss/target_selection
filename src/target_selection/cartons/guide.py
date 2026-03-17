#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-26
# @Filename: guide.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToTIC_v8,
    Gaia_DR2,
    Gaia_DR2_Clean,
    TIC_v8,
)

from .base import BaseCarton


Gaia = Gaia_DR2
Gaia_Clean = Gaia_DR2_Clean


class GuideCarton(BaseCarton):
    name = "guide"
    category = "guide"
    cadence = None

    def build_query(self, version_id, query_region=None):
        sample = (
            Catalog.select(Catalog.catalogid, Catalog.ra, Catalog.dec, Gaia.phot_g_mean_mag)
            .join(CatalogToTIC_v8)
            .join(TIC_v8)
            .join(Gaia)
            .join(Gaia_Clean)
            .where(
                (Gaia.phot_g_mean_mag > self.parameters["g_min"])
                & (Gaia.phot_g_mean_mag < self.parameters["g_max"])
            )
            .where(Catalog.version_id == version_id, CatalogToTIC_v8.version_id == version_id)
        )

        if query_region:
            sample = sample.where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        sample = sample.cte("sample")

        # We should use q3c_join_pm and q3c_dist_pm here.
        subq = (
            Gaia.select(Gaia.source_id)
            .where(
                peewee.fn.q3c_join(
                    sample.c.ra, sample.c.dec, Gaia.ra, Gaia.dec, self.parameters["min_separation"]
                )
            )
            .order_by(peewee.fn.q3c_dist(sample.c.ra, sample.c.dec, Gaia.ra, Gaia.dec).asc())
            .limit(1)
            .offset(1)
        )

        decollided = (
            sample.select(sample.c.catalogid)
            .join(subq, peewee.JOIN.LEFT_LATERAL, on=peewee.SQL("true"))
            .where(subq.c.source_id.is_null())
            .with_cte(sample)
        )

        return decollided
