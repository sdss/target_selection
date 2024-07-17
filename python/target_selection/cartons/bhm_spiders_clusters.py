#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-17
# @Filename: bhm_spiders_clusters.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn


from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import AB2nMgy  # , AB2Jy

# general catalogdb imports
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    EROSITASupersetv1Clusters,
)

# imports of existing spectro catalogues
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogFromSDSS_DR19p_Speclite,
    SDSS_DR19p_Speclite,
)

# additional imports required by bhm_spiders_clusters_lsdr10
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToLegacy_Survey_DR10,
    Legacy_Survey_DR10,
)

# # additional imports required by bhm_spiders_clusters_ps1dr2
# from sdssdb.peewee.sdss5db.catalogdb import (
#     Panstarrs1,
#     CatalogToPanstarrs1,    # only exists after v0.5 cross-match
# )

# DEBUG STUFF TO USE TEMP TABLE
# CatalogToSDSS_DR19p_Speclite._meta.table_name = 'temp_catalog_to_sdss_dr19p_speclite'
# CatalogToSDSS_DR19p_Speclite._meta._schema = 'sandbox'

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms

# ############################################
# ############################################
# ############################################
# ############################################
# # This file provides the following BHM cartons in v1.0:
#
#   bhm_spiders_clusters_lsdr10
#   bhm_spiders_clusters_lsdr10_d3
#
# ############################################
# ############################################
# ############################################
# ############################################

# Notes on how many targets to expect:
# => SELECT ero_version,xmatch_method,xmatch_version,opt_cat,count(*)
#      FROM erosita_superset_v1_clusters GROUP BY ero_version,xmatch_method,xmatch_version,opt_cat;
#
#  ero_version | xmatch_method | xmatch_version  | opt_cat |  count
# -------------+---------------+-----------------+---------+---------
#  eRASS4_020  | eromapper     | lsdr10a_grz_z_v | lsdr10  | 2330083

#
#
# END PREAMBLE
# ##################################################################################


class BhmSpidersClustersLsdr10Carton(BaseCarton):
    name = "bhm_spiders_clusters_lsdr10"
    category = "science"
    mapper = "BHM"
    program = "bhm_spiders"
    tile = False
    instrument = "BOSS"
    inertial = True
    can_offset = True
    only_faintest_cadence = False

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        ls = Legacy_Survey_DR10.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()

        x = EROSITASupersetv1Clusters.alias()
        # xx = EROSITASupersetv1Clusters.alias()
        # x = (
        #     xx
        #     .select(
        #         fn.rank().over(partition_by=[xx.ero_detuid],
        #                        order_by=[xx.xmatch_metric.desc()]).alias('x_rank'),
        #         xx.ero_detuid.alias('ero_detuid'),
        #         xx.ls_id.alias('ls_id'),
        #         xx.target_has_spec.alias('target_has_spec'),
        #     )
        #     .where(
        #         (xx.ero_version == self.parameters['ero_version']),
        #         (xx.xmatch_method == self.parameters['xmatch_method']),
        #         (xx.xmatch_version == self.parameters['xmatch_version']),
        #         (xx.opt_cat == self.parameters['opt_cat']),
        #         (xx.xmatch_metric > self.parameters['xmatch_metric_min']),
        #         (xx.ero_det_like > self.parameters['det_like_min']),
        #     )
        #     .alias('x')
        # )

        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast("bool")

        fibertotflux_i_max = AB2nMgy(self.parameters["fibertotmag_i_min"])
        fibertotflux_i_min = AB2nMgy(self.parameters["fibertotmag_i_max"])
        fibertotflux_z_max = AB2nMgy(self.parameters["fibertotmag_z_min"])
        fibertotflux_z_min = AB2nMgy(self.parameters["fibertotmag_z_max"])
        fibertotflux_i_max_for_core = AB2nMgy(self.parameters["fibertotmag_i_min_for_core"])
        fibertotflux_i_min_for_core = AB2nMgy(self.parameters["fibertotmag_i_max_for_core"])
        fibertotflux_z_max_for_core = AB2nMgy(self.parameters["fibertotmag_z_min_for_core"])
        fibertotflux_z_min_for_core = AB2nMgy(self.parameters["fibertotmag_z_max_for_core"])

        fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters["fibertotmag_r_for_cadence1"])
        fibertotflux_z_min_for_cadence1 = AB2nMgy(self.parameters["fibertotmag_z_for_cadence1"])
        fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters["fibertotmag_r_for_cadence2"])
        gaia_g_max_for_cadence1 = self.parameters["gaia_g_max_for_cadence1"]
        gaia_rp_max_for_cadence1 = self.parameters["gaia_rp_max_for_cadence1"]

        # flux30 = AB2nMgy(30.00)

        # #########################################################################
        # prepare the spectroscopy catalogues
        spec_sn_thresh = self.parameters["spec_sn_thresh"]
        spec_z_err_thresh = self.parameters["spec_z_err_thresh"]

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
        # #########################################################################

        # logic is written this backwards way so that a failure to meet ant core
        # criterion results in non-core status
        is_core = peewee.Case(
            None,
            (
                (
                    ~(
                        (
                            ls.fiberflux_r.between(
                                fibertotflux_i_min_for_core, fibertotflux_i_max_for_core
                            )
                        )
                        | (
                            ls.fiberflux_i.between(
                                fibertotflux_z_min_for_core, fibertotflux_z_max_for_core
                            )
                        )
                    ),
                    False,
                ),
            ),
            True,
        )

        # priority is determined by target rank within cluster
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:

        # priority = peewee.Case(
        #     None,
        #     (
        #         (
        #             x.c.x_rank == 1,
        #             self.parameters['priority_floor_bcg']
        #         ),
        #         (
        #             x.c.x_rank > 1,
        #             self.parameters['priority_floor_member'] +
        #             fn.least(self.parameters['priority_levels'] - 2,
        #                      x.c.x_rank - 2)
        #         ),
        #     ),
        #     None)

        dpriority_non_core = peewee.Case(
            None, ((is_core, 0),), self.parameters["dpriority_non_core"]
        )

        priority_floor = peewee.Case(
            None,
            ((x.target_priority == 0, self.parameters["priority_floor_bcg"]),),
            self.parameters["priority_floor_member"],
        )

        priority = priority_floor + x.target_priority + dpriority_non_core

        value = peewee.Case(
            None,
            (
                (~is_core, 0),
                (x.target_priority == 0, self.parameters["value_bcg"]),
                (x.target_priority > 0, self.parameters["value_member"]),
            ),
            None,
        ).cast("float")

        # choose cadence based on fiber magnitude in r-band
        cadence1 = self.parameters["cadence1"]
        cadence2 = self.parameters["cadence2"]
        cadence3 = self.parameters["cadence3"]
        cadence4 = "unknown_cadence"  # catch failures
        cadence = peewee.Case(
            None,
            (
                (
                    (
                        (ls.fibertotflux_r > fibertotflux_r_min_for_cadence1)
                        | (ls.fibertotflux_z > fibertotflux_z_min_for_cadence1)
                        | (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1))
                        | (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))
                    ),
                    cadence1,
                ),
                (ls.fibertotflux_r > fibertotflux_r_min_for_cadence2, cadence2),
                (ls.fibertotflux_r <= fibertotflux_r_min_for_cadence2, cadence3),
            ),
            cadence4,
        )

        # compute transformed SDSS mags for pointlike and extended sources uniformly
        # transform the legacysurvey grz into sdss psfmag griz

        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_clusters_lsdr8/lsdr8_fibermag_to_sdss_fiber2mag_?_results.log   # noqa
        coeffs = {
            "g2_e": -0.897719,
            "g1_e": 2.298300,
            "g0_e": -1.019299,
            "i2_e": -0.950114,
            "i1_e": 0.981972,
            "i0_e": -0.261645,
            "r2_e": -0.201741,
            "r1_e": 0.697128,
            "r0_e": -0.120926,
            "z2_e": -1.424312,
            "z1_e": 2.415301,
            "z0_e": -0.677163,
        }

        nMgy_min = 1e-3  # equiv to AB=30
        # extended - start from ls8 fiberfluxes
        g0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g))
        r0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r))
        i0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_i))
        z0_e = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_z))

        g0_t = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fibertotflux_g))
        r0_t = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fibertotflux_r))
        i0_t = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fibertotflux_i))
        z0_t = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fibertotflux_z))

        g0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g))
        r0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r))
        i0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_i))
        z0_p = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z))

        g_r_e = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.fiberflux_g)
            / peewee.fn.greatest(nMgy_min, ls.fiberflux_r)
        )
        r_z_e = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.fiberflux_r)
            / peewee.fn.greatest(nMgy_min, ls.fiberflux_z)
        )

        g_e = g0_e + coeffs["g0_e"] + coeffs["g1_e"] * g_r_e + coeffs["g2_e"] * g_r_e * g_r_e
        r_e = r0_e + coeffs["r0_e"] + coeffs["r1_e"] * g_r_e + coeffs["r2_e"] * g_r_e * g_r_e
        i_e = r0_e + coeffs["i0_e"] + coeffs["i1_e"] * r_z_e + coeffs["i2_e"] * r_z_e * r_z_e
        z_e = z0_e + coeffs["z0_e"] + coeffs["z1_e"] * r_z_e + coeffs["z2_e"] * r_z_e * r_z_e

        # validity checks
        valid = g0_e.between(0.1, 29.9) & r0_e.between(0.1, 29.9) & z0_e.between(0.1, 29.9)

        opt_prov = peewee.Case(None, ((valid, "sdss_fiber2mag_from_lsdr10"),), "undefined")
        magnitude_g = peewee.Case(None, ((valid, g_e),), "NaN")
        magnitude_r = peewee.Case(None, ((valid, r_e),), "NaN")
        magnitude_i = peewee.Case(None, ((valid, i_e),), "NaN")
        magnitude_z = peewee.Case(None, ((valid, z_e),), "NaN")
        magnitude_gaia_g = peewee.Case(
            None, ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),), "NaN"
        )
        magnitude_gaia_bp = peewee.Case(
            None, ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),), "NaN"
        )
        magnitude_gaia_rp = peewee.Case(
            None, ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),), "NaN"
        )

        spec_sn_thresh = self.parameters["spec_sn_thresh"]
        spec_z_err_thresh = self.parameters["spec_z_err_thresh"]

        query = (
            c.select(
                c.catalogid.alias("catalogid"),
                ls.ls_id.alias("ls_id"),  # extra
                x.pkey.alias("ero_pkey"),  # extra
                x.ero_detuid.cast("text").alias("ero_detuid"),  # extra
                s19.c.s19_pk.alias("sdss_dr19p_speclite_pk"),  # extra
                c.ra.alias("ra"),  # extra
                c.dec.alias("dec"),  # extra
                priority.alias("priority"),
                value.alias("value"),
                cadence.alias("cadence"),
                instrument.alias("instrument"),
                opt_prov.alias("optical_prov"),
                magnitude_g.alias("g"),
                magnitude_r.alias("r"),
                magnitude_i.alias("i"),
                magnitude_z.alias("z"),
                magnitude_gaia_g.alias("gaia_g"),
                magnitude_gaia_bp.alias("bp"),
                magnitude_gaia_rp.alias("rp"),
                inertial.alias("inertial"),
                g0_p.alias("ls10_mag_g"),  # extra
                r0_p.alias("ls10_mag_r"),  # extra
                i0_p.alias("ls10_mag_i"),  # extra
                z0_p.alias("ls10_mag_z"),  # extra
                g0_e.alias("ls10_fibermag_g"),  # extra
                r0_e.alias("ls10_fibermag_r"),  # extra
                i0_e.alias("ls10_fibermag_i"),  # extra
                z0_e.alias("ls10_fibermag_z"),  # extra
                g0_t.alias("ls10_fibertotmag_g"),  # extra
                r0_t.alias("ls10_fibertotmag_r"),  # extra
                i0_t.alias("ls10_fibertotmag_i"),  # extra
                z0_t.alias("ls10_fibertotmag_z"),  # extra
                is_core.alias("is_core"),  # extra
                x.target_priority.alias("orig_target_priority"),
                x.eromapper_lambda.alias("eromapper_lambda"),
                x.eromapper_z_lambda.alias("eromapper_z_lambda"),
                x.xmatch_flags.alias("xmatch_flags"),
                x.target_has_spec.alias("target_has_spec"),
            )
            .join(c2ls)
            .join(ls)
            .join(x, on=(ls.ls_id == x.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s19, JOIN.LEFT_OUTER)
            .join(s19, JOIN.LEFT_OUTER, on=(s19.c.s19_pk == c2s19.target_id))
            # finished joining the spectroscopy
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True,
                fn.coalesce(c2s19.version_id, version_id) == version_id,
            )
            .where(
                s19.c.s19_pk.is_null(True),
            )
            .where(
                (x.ero_version == self.parameters["ero_version"]),
                (x.xmatch_method == self.parameters["xmatch_method"]),
                (x.xmatch_version == self.parameters["xmatch_version"]),
                (x.opt_cat == self.parameters["opt_cat"]),
                (x.xmatch_metric > self.parameters["xmatch_metric_min"]),
                (x.ero_det_like > self.parameters["det_like_min"]),
            )
            .where(
                (
                    (ls.fibertotflux_i.between(fibertotflux_i_min, fibertotflux_i_max))
                    | (ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max))
                ),
                (x.target_has_spec == 0),
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters["gaia_g_mag_limit"])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters["gaia_rp_mag_limit"])),
            )
            .distinct([ls.ls_id])
        )

        if self.only_faintest_cadence:
            query = query.where(cadence == cadence3)

        if query_region:
            query = query.where(
                peewee.fn.q3c_radial_query(
                    c.ra, c.dec, query_region[0], query_region[1], query_region[2]
                )
            )

        return query


#
# END BhmSpidersClustersLsdr10Carton
# ##################################################################################


class BhmSpidersClustersLsdr10D3Carton(BhmSpidersClustersLsdr10Carton):
    name = "bhm_spiders_clusters_lsdr10_d3"
    only_faintest_cadence = True
