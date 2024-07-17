#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-12
# @Filename: bhm_csc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import JOIN, fn

from sdssdb.peewee.sdss5db.catalogdb import (
    BHM_CSC_v3,
    CatalogFromSDSS_DR19p_Speclite,
    CatalogToGaia_DR3,
    CatalogToLegacy_Survey_DR10,
    CatalogToPanstarrs1,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    Legacy_Survey_DR10,
    Panstarrs1,
    SDSS_DR19p_Speclite,
    TwoMassPSC,
)

from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import AB2Jy, AB2nMgy


# DEBUG STUFF TO USE TEMP TABLE
# CatalogToSDSS_DR19p_Speclite._meta.table_name = 'temp_catalog_to_sdss_dr19p_speclite'
# CatalogToSDSS_DR19p_Speclite._meta._schema = 'sandbox'

# Details: Start here
# https://wiki.sdss.org/display/OPS/Cartons+for+v0.5#Cartonsforv0.5-BHMCSC # noqa: E501

# This module provides the following BHM cartons:
#   bhm_csc_boss
#   bhm_csc_boss_d3
#   bhm_csc_apogee


class BhmCscBossCarton(BaseCarton):
    """ """

    name = "bhm_csc_boss"
    category = "science"
    mapper = "BHM"
    program = "bhm_csc"
    instrument = "BOSS"
    tile = False
    can_offset = True
    only_faintest_cadence = False

    def build_query(self, version_id, query_region=None):
        x = BHM_CSC_v3.alias()
        ps = Panstarrs1.alias()
        c2ps = CatalogToPanstarrs1.alias()
        g3 = Gaia_DR3.alias()
        c2g3 = CatalogToGaia_DR3.alias()
        ls = Legacy_Survey_DR10.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()

        g_psf_flux_max = AB2Jy(self.parameters["g_psf_mag_min"])
        r_psf_flux_max = AB2Jy(self.parameters["r_psf_mag_min"])
        i_psf_flux_max = AB2Jy(self.parameters["i_psf_mag_min"])
        z_psf_flux_max = AB2Jy(self.parameters["z_psf_mag_min"])
        g_psf_flux_min = AB2Jy(self.parameters["g_psf_mag_max"])
        r_psf_flux_min = AB2Jy(self.parameters["r_psf_mag_max"])
        i_psf_flux_min = AB2Jy(self.parameters["i_psf_mag_max"])
        z_psf_flux_min = AB2Jy(self.parameters["z_psf_mag_max"])
        fiberflux_g_max = AB2nMgy(self.parameters["fibermag_g_min"])
        fiberflux_r_max = AB2nMgy(self.parameters["fibermag_r_min"])
        fiberflux_i_max = AB2nMgy(self.parameters["fibermag_i_min"])
        fiberflux_z_max = AB2nMgy(self.parameters["fibermag_z_min"])
        fiberflux_g_min = AB2nMgy(self.parameters["fibermag_g_max"])
        fiberflux_r_min = AB2nMgy(self.parameters["fibermag_r_max"])
        fiberflux_i_min = AB2nMgy(self.parameters["fibermag_i_max"])
        fiberflux_z_min = AB2nMgy(self.parameters["fibermag_z_max"])

        gaia_g_mag_min = self.parameters["gaia_g_mag_min"]
        gaia_g_mag_max = self.parameters["gaia_g_mag_max"]

        i_psf_flux_min_for_cadence1 = AB2Jy(self.parameters["i_psf_mag_max_for_cadence1"])
        fiberflux_r_min_for_cadence1 = AB2nMgy(self.parameters["fibermag_r_max_for_cadence1"])
        gaia_g_mag_max_for_cadence1 = self.parameters["gaia_g_mag_max_for_cadence1"]
        i_psf_flux_min_for_cadence2 = AB2Jy(self.parameters["i_psf_mag_max_for_cadence2"])
        fiberflux_r_min_for_cadence2 = AB2nMgy(self.parameters["fibermag_r_max_for_cadence2"])
        gaia_g_mag_max_for_cadence2 = self.parameters["gaia_g_mag_max_for_cadence2"]

        value = peewee.Value(self.parameters["value"]).cast("real")

        cadence1 = self.parameters["cadence1"]
        cadence2 = self.parameters["cadence2"]
        cadence3 = self.parameters["cadence3"]
        cadence = peewee.Case(
            None,
            (
                (
                    (
                        (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1)
                        | (ls.fiberflux_r > fiberflux_r_min_for_cadence1)
                        | (g3.phot_g_mean_mag < gaia_g_mag_max_for_cadence1)
                    ),
                    cadence1,
                ),
                (
                    (
                        (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence2)
                        | (ls.fiberflux_r > fiberflux_r_min_for_cadence2)
                        | (g3.phot_g_mean_mag < gaia_g_mag_max_for_cadence2)
                    ),
                    cadence2,
                ),
            ),
            cadence3,
        )

        # #########################################################################
        # prepare the existing spectroscopy catalogues - these affect only priority
        spec_sn_thresh = self.parameters["spec_sn_thresh"]
        spec_z_err_thresh = self.parameters["spec_z_err_thresh"]
        dpriority_has_spec = self.parameters["dpriority_has_spec"]

        # SDSS DR19p
        # downslect only 'good' spectra
        c2s19 = CatalogFromSDSS_DR19p_Speclite.alias()
        ss19 = SDSS_DR19p_Speclite.alias()
        s19 = (
            ss19.select(
                ss19.pk.alias("s19_pk"),
            )
            .where(
                ss19.sn_median_all >= spec_sn_thresh,
                ss19.zwarning == 0,
                ss19.z_err <= spec_z_err_thresh,
                ss19.z_err > 0.0,
                ss19.specprimary > 0,
            )
            .alias("s19")
        )

        # adjust priority if target aleady has a good SDSS spectrum
        priority_1 = peewee.Case(None, ((s19.c.s19_pk.is_null(False), 1),), 0)
        #
        # Compute net priority
        priority_floor_dark = peewee.Value(self.parameters["priority_floor_dark"])
        priority_floor_bright = peewee.Value(self.parameters["priority_floor_bright"])
        priority = peewee.Case(
            None,
            (
                (
                    (
                        (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1)
                        | (g3.phot_g_mean_mag < gaia_g_mag_max_for_cadence1)
                    ),
                    priority_floor_bright + (priority_1 * dpriority_has_spec) + x.xpriority - 1,
                ),
            ),
            priority_floor_dark + (priority_1 * dpriority_has_spec) + x.xpriority - 1,
        )

        query = (
            x.select(
                fn.coalesce(c2g3.catalogid, c2ls.catalogid, c2ps.catalogid).alias("catalogid"),
                priority.alias("priority"),
                cadence.alias("cadence"),
                value.alias("value"),
                x.pk.alias("csc_pk"),  # extra
                x.csc21p_id.alias("csc_csc21p_id"),  # extra
                x.ra.alias("csc_opt_ra"),  # extra
                x.dec.alias("csc_opt_dec"),  # extra
                x.best_mag.alias("csc_best_mag"),  # extra
                x.mag_type.alias("csc_mag_type"),  # extra
                x.best_oir_cat.alias("csc_best_oir_cat"),  # extra
                x.xpriority.alias("csc_xpriority"),  # extra
                x.xband.alias("csc_xband"),  # extra
                x.logfx.alias("csc_logfx"),  # extra
                x.gaia_dr3_srcid.alias("csc_gaia_dr3_srcid"),  # extra
                x.ls_dr10_lsid.alias("csc_ls_dr10_lsid"),  # extra
                x.ps21p_ippobjid.alias("csc_ps21p_ippobjid"),  # extra
                g3.phot_g_mean_mag.alias("g3_g_mag"),  # extra
                ls.fiberflux_r.alias("ls_fiberflux_r"),  # extra
                ls.fiberflux_i.alias("ls_fiberflux_i"),  # extra
                ps.r_stk_psf_flux.alias("ps_r_stk_psf_flux"),  # extra
                ps.i_stk_psf_flux.alias("ps_i_stk_psf_flux"),  # extra
                c2s19.target_id.alias("sdss_dr19p_speclite_pk"),  # extra
            )
            .join(g3, join_type=JOIN.LEFT_OUTER, on=(x.gaia_dr3_srcid == g3.source_id))
            .join(ls, join_type=JOIN.LEFT_OUTER, on=(x.ls_dr10_lsid == ls.ls_id))
            .join(ps, join_type=JOIN.LEFT_OUTER, on=(x.ps21p_ippobjid == ps.catid_objid))
            .join(c2ps, join_type=JOIN.LEFT_OUTER, on=(ps.catid_objid == c2ps.target_id))
            .join(c2ls, join_type=JOIN.LEFT_OUTER, on=(ls.ls_id == c2ls.target_id))
            .join(c2g3, join_type=JOIN.LEFT_OUTER, on=(g3.source_id == c2g3.target_id))
            .join(
                c2s19,
                join_type=JOIN.LEFT_OUTER,
                on=(
                    fn.coalesce(c2g3.catalogid, c2ls.catalogid, c2ps.catalogid) == c2s19.catalogid
                ),
            )
            .join(s19, join_type=JOIN.LEFT_OUTER, on=(s19.c.s19_pk == c2s19.target_id))
            .where(
                fn.coalesce(c2g3.version_id, c2ls.version_id, c2ps.version_id) == version_id,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
                fn.coalesce(c2g3.best, c2ls.best, c2ps.best) >> True,
            )
            .where(
                (
                    (x.best_oir_cat == "gdr3")
                    & (g3.phot_g_mean_mag.between(gaia_g_mag_min, gaia_g_mag_max))
                )
                | (
                    (x.best_oir_cat == "lsdr10")
                    & (ls.fiberflux_g < fiberflux_g_max)
                    & (ls.fiberflux_r < fiberflux_r_max)
                    & (ls.fiberflux_i < fiberflux_i_max)
                    & (ls.fiberflux_z < fiberflux_z_max)
                    & (
                        (ls.fiberflux_g > fiberflux_g_min)
                        | (ls.fiberflux_r > fiberflux_r_min)
                        | (ls.fiberflux_i > fiberflux_i_min)
                        | (ls.fiberflux_z > fiberflux_z_min)
                    )
                )
                | (
                    (x.best_oir_cat == "ps1dr2")
                    & (ps.g_stk_psf_flux < g_psf_flux_max)
                    & (ps.r_stk_psf_flux < r_psf_flux_max)
                    & (ps.i_stk_psf_flux < i_psf_flux_max)
                    & (ps.z_stk_psf_flux < z_psf_flux_max)
                    & (
                        (ps.g_stk_psf_flux > g_psf_flux_min)
                        | (ps.r_stk_psf_flux > r_psf_flux_min)
                        | (ps.i_stk_psf_flux > i_psf_flux_min)
                        | (ps.z_stk_psf_flux > z_psf_flux_min)
                    )
                )
            )
            # .distinct(fn.coalesce(fn.max(c2ps.catalogid), fn.max(c2g3.catalogid)))
            .distinct([fn.coalesce(c2g3.catalogid, c2ls.catalogid, c2ps.catalogid)])
        )

        if self.only_faintest_cadence:
            query = query.where(cadence == cadence3)

        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    x.ora, x.odec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


class BhmCscBossD3Carton(BhmCscBossCarton):
    name = "bhm_csc_boss_d3"
    only_faintest_cadence = True


# #######
class BhmCscApogeeCarton(BaseCarton):
    """
    SELECT * from bhm_csc_v2 AS x
    WHERE
      AND x.best_mag BETWEEN 7 AND 14
    """

    name = "bhm_csc_apogee"
    mapper = "BHM"
    program = "bhm_csc"
    category = "science"
    this_cadence = "bright_3x1"  # TODO Check this is still a valid choice
    instrument = "APOGEE"
    tile = False
    can_offset = True

    def build_query(self, version_id, query_region=None):
        x = BHM_CSC_v3.alias()
        # tic = TIC_v8.alias()
        # c2tic = CatalogToTIC_v8.alias()
        # change to now rely directly on twomass_psc table instead of TIC_v8
        tm = TwoMassPSC.alias()
        c2tm = CatalogToTwoMassPSC.alias()

        hmag_max_for_cadence1 = self.parameters["hmag_max_for_cadence1"]

        value = peewee.Value(self.parameters["value"]).cast("real")

        cadence = peewee.Case(
            None,
            (((x.best_mag < hmag_max_for_cadence1), self.parameters["cadence1"]),),
            self.parameters["cadence2"],
        )

        # Compute net priority
        priority = peewee.Value(self.parameters["priority_floor"]) + x.xpriority - 1

        query = (
            x.select(
                # c2tic.catalogid.alias('catalogid'),
                c2tm.catalogid.alias("catalogid"),
                x.csc21p_id.alias("csc_csc21p_id"),  # extra
                priority.alias("priority"),
                cadence.alias("cadence"),
                value.alias("value"),
                tm.j_m.alias("j"),
                tm.h_m.alias("h"),
                tm.k_m.alias("k"),
                tm.pts_key.alias("twomass_pts_key"),
                x.best_oir_cat.alias("csc_best_oir_cat"),  # extra
                x.ra.alias("csc_oir_ra"),  # extra
                x.dec.alias("csc_oir_dec"),  # extra
                x.best_mag.alias("csc_best_mag"),  # extra
                x.mag_type.alias("csc_mag_type"),  # extra
                x.xpriority.alias("csc_xpriority"),  # extra
                x.xband.alias("csc_xband"),  # extra
                x.logfx.alias("csc_logfx"),  # extra
                x.tmass_designation.alias("csc_tmass_designation"),  # extra
            )
            .join(tm, on=(x.tmass_designation == tm.designation))
            .join(c2tm, on=(tm.pts_key == c2tm.target_id))
            .where(
                c2tm.version_id == version_id,
                c2tm.best >> True,
                # x.best_oir_cat == '2mass',
                # x.best_mag >= self.parameters['hmag_min'],
                # x.best_mag < self.parameters['hmag_max'],
                # x.best_mag != 'NaN',
                tm.h_m.between(self.parameters["hmag_min"], self.parameters["hmag_max"]),
            )
            .distinct([c2tm.catalogid])
        )

        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    x.ra2m, x.dec2m, query_region[0], query_region[1], query_region[2]
                )
            )

        return query
