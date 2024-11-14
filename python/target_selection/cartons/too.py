#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-02-16
# @Filename: too.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

from sdssdb.peewee.sdss5db.catalogdb import CatalogToToO_Target, ToO_Metadata, ToO_Target
from sdssdb.peewee.sdss5db.targetdb import Carton, CartonToTarget, Target, Version

from .base import BaseCarton


__all__ = ["ToO_Carton"]


class ToO_Carton(BaseCarton):
    """Target of opportunity carton.

    Selects all the targets in ``catalogdb.too_target`` that don't yet exist in
    the carton.

    """

    name = "too"
    category = "science"
    cadence = "bright_1x1"
    priority = 3000
    program = "too"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        C2TT = CatalogToToO_Target

        too_in_carton = (
            Target.select(Target.catalogid)
            .join(CartonToTarget)
            .join(Carton)
            .join(Version)
            .where(Version.plan == self.plan, Carton.carton == self.name)
        ).alias("too_in_carton")

        query = (
            ToO_Target.select(
                C2TT.catalogid,
                ToO_Target.fiber_type.alias("instrument"),
                ToO_Metadata.g_mag.alias("g"),
                ToO_Metadata.r_mag.alias("r"),
                ToO_Metadata.i_mag.alias("i"),
                ToO_Metadata.z_mag.alias("z"),
                ToO_Metadata.optical_prov,
            )
            .join(C2TT, on=(ToO_Target.too_id == C2TT.target_id))
            .switch(ToO_Target)
            .join(ToO_Metadata, on=(ToO_Target.too_id == ToO_Metadata.too_id))
            .where(
                C2TT.version_id == version_id,
                C2TT.best >> True,
                C2TT.catalogid.not_in(too_in_carton),
            )
        )

        return query
