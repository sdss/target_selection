#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-31
# @Filename: bhm_rm.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn

from target_selection.cartons.base import BaseCarton
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    BHM_RM_v1_3,
    CatalogToLegacy_Survey_DR8,
    CatalogToLegacy_Survey_DR10,
    CatalogToGaia_DR2,
    CatalogToGaia_DR3,
    CatalogToPanstarrs1,
)


#  This module provides the following BHM cartons in v1.0:
#  bhm_rm_core
#  bhm_rm_known_spec
#  bhm_rm_var
#  bhm_rm_ancillary
#  bhm_rm_xrayqso


class BhmRmBaseCarton(BaseCarton):
    """
    This class provides common setting and the masking routines used by all RM cartons
    """

    name = "bhm_rm_base"
    base_name = "bhm_rm_base"
    category = "science"
    mapper = "BHM"
    program = "bhm_rm"
    instrument = "BOSS"
    tile = False
    priority = None
    inertial = True
    alias_c = None
    alias_t = None
    alias_tw = None
    can_offset = False

    def get_fieldlist(self):
        """Read the RM field centres from the yaml"""
        fieldlist = []
        base_parameters = self.config["parameters"].get(self.base_name, None)
        if base_parameters:
            fieldlist = base_parameters["fieldlist"]
        return fieldlist

    def append_spatial_query(self, query, fieldlist):
        """extend the peewee query using a list of field centres"""
        if fieldlist is None:
            return query
        elif len(fieldlist) == 0:
            return query

        q = False
        for f in fieldlist:
            q = q | peewee.fn.q3c_radial_query(
                self.alias_c.ra, self.alias_c.dec, f["racen"], f["deccen"], f["radius"]
            )
        return query.where(q)

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2ls8 = CatalogToLegacy_Survey_DR8.alias()
        c2ls10 = CatalogToLegacy_Survey_DR10.alias()
        c2g3 = CatalogToGaia_DR3.alias()
        c2g2 = CatalogToGaia_DR2.alias()
        c2ps = CatalogToPanstarrs1.alias()
        t = BHM_RM_v1_3.alias()
        self.alias_c = c
        self.alias_t = t

        fieldlist = self.get_fieldlist()

        # fold in tiers of magnitude-based priority
        priority_mag_step = 0.5
        priority_mag_bright = 17.0
        priority_mag_faint = 22.0
        priority_mag_bright_known_spec = 20.5
        priority_floor = self.parameters.get("priority", 10000)
        priority1 = peewee.Case(
            None,
            (
                ((t.mag_i <= priority_mag_bright), priority_floor + 0),
                (
                    (
                        (self.name == "bhm_rm_known_spec")
                        & ~(t.rm_field_name.contains("SDSS-RM"))
                        & (t.mag_i <= priority_mag_bright_known_spec)
                    ),
                    priority_floor + 0,
                ),
                (
                    (t.mag_i <= priority_mag_faint),
                    priority_floor
                    + 5
                    * (
                        1
                        + peewee.fn.floor(
                            (t.mag_i - priority_mag_bright) / priority_mag_step
                        ).cast("int")
                    ),
                ),
                ((t.mag_i > priority_mag_faint), priority_floor + 95),
            ),
            None,
        )

        # combine the priorities
        priority = priority1

        value = peewee.Value(self.parameters.get("value", 1.0)).cast("float")
        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast("bool")

        # This is the scheme used in v0
        cadence_v0 = peewee.Case(
            None, ((t.rm_field_name.contains("S-CVZ"), "bhm_rm_lite5_100x8"),), "bhm_rm_174x8"
        )

        # this gives the new names for the same cadences assumed in v0
        cadence_v0p5 = peewee.Case(
            None, ((t.rm_field_name.contains("S-CVZ"), "dark_100x8"),), "dark_174x8"
        )

        # the following will replace old generic cadences when relevant table has been populated
        # TODO - replace when correct cadences are loaded
        cadence_v1p0 = peewee.Case(
            None,
            (
                (t.rm_field_name.contains("SDSS-RM"), "bhm_rm_sdss-rm"),
                (t.rm_field_name.contains("COSMOS"), "bhm_rm_cosmos"),
                (t.rm_field_name.contains("XMM-LSS"), "bhm_rm_xmm-lss"),
                (t.rm_field_name.contains("S-CVZ"), "bhm_rm_cvz-s"),
                (t.rm_field_name.contains("CDFS"), "bhm_rm_cdfs"),
                (t.rm_field_name.contains("ELIAS-S1"), "bhm_rm_elias-s1"),
            ),
            "dark_174x8",
        )

        opt_prov = peewee.Value("psfmag")

        route_taken = peewee.Case(
            None,
            (
                (c2ls10.catalogid.is_null(False), "lsdr10"),
                (c2g3.catalogid.is_null(False), "gdr3"),
                (c2ls8.catalogid.is_null(False), "lsdr8"),
                (c2g2.catalogid.is_null(False), "gdr2"),
                (c2ps.catalogid.is_null(False), "ps1dr2"),
            ),
            "unknown",
        )
        query = (
            t.select(
                c.catalogid,
                c2ls10.catalogid.alias("c2ls10_catalogid"),  # extra
                c2g3.catalogid.alias("c2g3_catalogid"),  # extra
                c2ls8.catalogid.alias("c2ls8_catalogid"),  # extra
                c2g2.catalogid.alias("c2g2_catalogid"),  # extra
                c2ps.catalogid.alias("c2ps_catalogid"),  # extra
                c.ra,  # extra
                c.dec,  # extra
                t.rm_field_name.alias("rm_field_name"),  # extra
                t.pkey.alias("rm_pkey"),  # extra
                instrument.alias("instrument"),
                priority.alias("priority"),
                value.alias("value"),
                cadence_v0p5.alias("cadence"),
                cadence_v0.alias("cadence_v0"),  # extra
                cadence_v0p5.alias("cadence_v0p5"),  # extra
                cadence_v1p0.alias("cadence_v1p0_maybe"),  # extra
                t.mag_g.alias("g"),
                t.mag_r.alias("r"),
                t.mag_i.alias("i"),
                t.mag_z.alias("z"),
                t.gaia_g.alias("gaia_g"),
                t.gaia_bp.alias("bp"),
                t.gaia_rp.alias("rp"),
                opt_prov.alias("optical_prov"),
                inertial.alias("inertial"),
                c2ls10.best.alias("c2ls10_best"),  # extra
                c2g3.best.alias("c2g3_best"),  # extra
                c2ls8.best.alias("c2ls8_best"),  # extra
                c2g2.best.alias("c2g2_best"),  # extra
                c2ps.best.alias("c2ps_best"),  # extra
                t.catalogidv05.alias("rm_catalogidv05"),  # extra
                t.ra.alias("rm_ra"),  # extra
                t.dec.alias("rm_dec"),  # extra
                route_taken.alias("route_taken"),  # extra
            )
            .join(c2ls10, JOIN.LEFT_OUTER, on=(c2ls10.target_id == t.ls_id_dr10))
            .join(c2g3, JOIN.LEFT_OUTER, on=(c2g3.target_id == t.gaia_dr3_source_id))
            .join(c2ls8, JOIN.LEFT_OUTER, on=(c2ls8.target_id == t.ls_id_dr8))
            .join(c2g2, JOIN.LEFT_OUTER, on=(c2g2.target_id == t.gaia_dr2_source_id))
            .join(c2ps, JOIN.LEFT_OUTER, on=(c2ps.target_id == t.panstarrs1_catid_objid))
            .join(
                c,
                on=(
                    fn.coalesce(
                        c2ls10.catalogid,
                        c2g3.catalogid,
                        c2ls8.catalogid,
                        c2g2.catalogid,
                        c2ps.catalogid,
                    )
                    == c.catalogid
                ),
            )
            .where(
                c.version_id == version_id,
                (
                    ((route_taken == "lsdr10") & (c2ls10.version_id == version_id))
                    | ((route_taken == "gdr3") & (c2g3.version_id == version_id))
                    | ((route_taken == "lsdr8") & (c2ls8.version_id == version_id))
                    | ((route_taken == "gdr2") & (c2g2.version_id == version_id))
                    | ((route_taken == "ps1dr2") & (c2ps.version_id == version_id))
                ),
                # the method below was throwing out ~20 cases where the lsdr8 ls_id
                # was unexpectedly unmatched to a catalogid in the v1 crossmatch
                # fn.coalesce(c2ls10.version_id, version_id) == version_id,
                # fn.coalesce(c2g3.version_id, version_id) == version_id,
                # fn.coalesce(c2ls8.version_id, version_id) == version_id,
                # fn.coalesce(c2g2.version_id, version_id) == version_id,
                # fn.coalesce(c2ps.version_id, version_id) == version_id,
                # ## fn.coalesce(c2ls10.best,True) >> True   # TODO check if this is dropping RM
                # ##                                         # targets like it does for AQMES
            )
            .where(
                (
                    (t.mag_i >= self.parameters["mag_i_min"])
                    & (t.mag_i < self.parameters["mag_i_max"])
                )
                | (
                    # S-CVZ targets often have only Gaia photom
                    (t.rm_field_name.contains("S-CVZ"))
                    & (t.gaia_g >= self.parameters["mag_g_min_cvz_s"])
                    & (t.gaia_g < self.parameters["mag_g_max_cvz_s"])
                )
            )
            .where(
                # Reject any objects where the t.rm_unsuitable flag is set
                t.rm_unsuitable >> False,
            )
            .distinct([t.pkey])  # avoid duplicates - trust the RM parent sample
            # - only needed if NOT using c2t.best = True condition
        )
        query = self.append_spatial_query(query, fieldlist)

        return query


class BhmRmCoreCarton(BhmRmBaseCarton):
    name = "bhm_rm_core"

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.rm_core >> True),
            ~(t.rm_field_name.contains("SDSS-RM")),  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmKnownSpecCarton(BhmRmBaseCarton):
    """
    bhm_rm_known_spec:  select all spectroscopically confirmed QSOs where redshift is extragalactic
    """

    name = "bhm_rm_known_spec"

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            ((t.rm_known_spec >> True),),
            (
                ~(t.rm_field_name.contains("SDSS-RM"))
                |
                # include extra constraints on SDSS-RM targets
                (t.mag_i < self.parameters["mag_i_max_sdss_rm"])
            ),
            (
                ~(t.rm_field_name.contains("COSMOS"))
                |
                # include extra constraints on COSMOS targets
                (t.mag_i < self.parameters["mag_i_max_cosmos"])
            ),
            (
                ~(t.rm_field_name.contains("XMM-LSS"))
                |
                # include extra constraints on XMM-LSS targets
                (t.mag_i < self.parameters["mag_i_max_xmm_lss"])
            ),
        )

        return query


class BhmRmVarCarton(BhmRmBaseCarton):
    """bhm_rm_var: selected based on g-band variability > 0.05 mag
    and bright enough to be detected by Gaia (G<~21)
    """

    name = "bhm_rm_var"

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.rm_var >> True),
            ~(t.rm_field_name.contains("SDSS-RM")),  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmAncillaryCarton(BhmRmBaseCarton):
    """
    bhm_rm_ancillary: from the Gaia_unWISE AGN catalog or the XDQSO catalog,
                      but requiring no proper motion/parallax detection from Gaia DR2
    """

    name = "bhm_rm_ancillary"

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t

        query = query.where(
            (t.rm_ancillary >> True),
            ~(t.rm_field_name.contains("SDSS-RM")),  # ignore this carton in the SDSS-RM field
        )

        return query


class BhmRmXrayQsoCarton(BhmRmBaseCarton):
    """
    bhm_rm_xrayqso:
       selected based on X-ray and SED

    """

    name = "bhm_rm_xrayqso"

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.rm_xrayqso > 0),
            ~(t.rm_field_name.contains("SDSS-RM")),  # ignore this carton in the SDSS-RM field
        )

        return query
