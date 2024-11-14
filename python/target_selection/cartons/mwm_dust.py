#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-16
# @Filename: mwm_dust.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import pandas
import peewee

from sdssdb.peewee.sdss5db import targetdb
from sdssdb.peewee.sdss5db.catalogdb import (
    GLIMPSE,
    AllWise,
    Catalog,
    CatalogToAllWise,
    CatalogToGaia_DR3,
    CatalogToGLIMPSE,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


def lbp2xyz(ll, bb, pp):
    dist = 1.0 / pp  # kpc
    z = dist * numpy.sin(numpy.radians(bb))
    x = dist * numpy.cos(numpy.radians(bb)) * numpy.cos(numpy.radians(ll))
    y = dist * numpy.cos(numpy.radians(bb)) * numpy.sin(numpy.radians(ll))

    # z is definitely right; x and y depend on conventions.  As long
    # as this routine is used for x and y throughout this code, this should
    # not introduce any issues.

    return x, y, z


# subselect() is used in post_process() to randomly select a subset
# from the results of the  query in build_query().
# Hence, the query in build_query() uses order_by().
# The order_by() in the query ensures that the random selection in
# post_process() gives the same result every time we run this carton.


def subselect(data, othersel, downsampledby=1):
    """Turn x, y, z into pixel; bin up to see where we are missing."""

    ll = data["gallong"]
    bb = data["gallat"]
    xyz = lbp2xyz(ll, bb, data["parallax"])

    dustsel = numpy.ones(len(data), numpy.bool)

    ranges = [[-5, 5], [-5, 5], [-0.2, 0.2]]
    resolution = 0.1

    countshape = [numpy.round((r[1] - r[0]) / resolution).astype("i4") + 2 for r in ranges]

    coords = [
        numpy.clip(numpy.floor((c - r[0]) / resolution).astype("i4"), -1, s - 2) + 1
        for (c, r, s) in zip(xyz, ranges, countshape)
    ]
    centers = [
        r[0] + resolution * (numpy.arange((r[1] - r[0]) / resolution) + 0.5) for r in ranges
    ]

    flatcoord = numpy.ravel_multi_index(coords, countshape)

    countsd = numpy.bincount(
        flatcoord[dustsel & ~othersel],
        minlength=numpy.product(countshape),
    )
    countso = numpy.bincount(
        flatcoord[othersel & dustsel],
        minlength=numpy.product(countshape),
    )
    countsd = countsd.reshape(countshape)
    countso = countso.reshape(countshape)

    pnumdensd = map_coordinates_wrap(
        (countsd[1:-1, 1:-1, 1:-1], centers),
        [c[dustsel & ~othersel] for c in xyz],
        order=0,
    )
    pnumdenso = map_coordinates_wrap(
        (countso[1:-1, 1:-1, 1:-1], centers),
        [c[dustsel & ~othersel] for c in xyz],
        order=0,
    )

    ntopick = 10.0 / downsampledby
    pnumdensleft = ntopick - pnumdenso

    numpy.random.seed(0)
    randnum = numpy.random.rand(numpy.sum(dustsel & ~othersel))

    subsel = randnum < pnumdensleft / 1.0 / (pnumdensd + (pnumdensd == 0))

    dustsel[othersel] = 0
    dustsel[dustsel & ~othersel] = subsel

    return dustsel


def map_coordinates_wrap(grid, coord, **kw):
    gridpts = grid[0]
    gridcoord = grid[1]
    outputcoordnorm = [
        numpy.interp(c, gc, numpy.arange(len(gc))) for c, gc in zip(coord, gridcoord)
    ]
    from scipy.ndimage import interpolation

    return interpolation.map_coordinates(
        gridpts,
        outputcoordnorm,
        cval=0.0,
        mode="constant",
        **kw,
    )


# These were initial selection criteria, not used anymore. We use the
# real GG carton instead.


def jkselect(data):
    jk0 = data["j_ks_0"]
    return (jk0 > 0.7) & (data["hmag"] < 11)


def ghselect(data):
    gaiag = data["gaiamag"]
    return (data["hmag"] < 11) & numpy.where(
        numpy.isfinite(gaiag),
        gaiag - data["hmag"] > 3.5,
        True,
    )


# The below MWM_Dust_Core_apogee_Carton carton is obsolete.
# Do not run this carton for v1.0.
# However, do not delete this carton since
# this carton is the base carton for MWM_Dust_Core_Dist_apogee_Carton.
# This carton requires that the MWM_Galatic_Core_apogee_Carton
# has been run before it in the same target_selection_plan
# in target_selection.yml or in an earlier target_selection_plan.
# See target_selection.yml for example for mwm_dust_core_apogee.
class MWM_Dust_Core_apogee_Carton(BaseCarton):
    """MWM Dust Core apogee Carton.

    Definition:

        A sample of bright, nearby, midplane giants that, when combined
        with the more distant and/or heavily reddened sample in other cartons,
        produces 100 stars per (100pc)^3 volume in the near disk.

        - Use Gaia parallax and (l,b) to get distance
        - Use H and 4.5 micron to get A_Ks and E_JKs, a la APOGEE
          (using RJCE dereddening prescription from Majewski et al. 2011).
        - Use J, K, and E_JKs to get (J-Ks)_0
        - Use distance, K, and A_Ks to get absolute mag M_K
        - 0 < Gaia parallax error/parallax < 0.2
        - M_K < 2.6
        - H < 11.2
        - distance < 5 kpc
        - |z| < 0.2 kpc
        - (J-Ks)_0 > 0.5
        - Can be spatially subselected to complement GG sampling in order to
          obtain 100 stars per (100pc)^3 volume.
        - Galactic Genesis quality flags:
            (ph_qual[1] = 'A' OR ph_qual[1] = 'B')
                AND gal_contam=0
                AND cc_flg[1]=0
                AND (rd_flg[1] > 0 AND rd_flg[1] <= 3)

    Non-SQL implementation:

        https://faun.rc.fas.harvard.edu/eschlafly/sdss5/dustsel.py

    """

    name = "mwm_dust_core_apogee"
    mapper = "MWM"
    category = "science"
    program = "mwm_dust"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    priority = 2720
    can_offset = True

    mwm_galactic_carton = None
    mwm_galactic_plan = None

    def build_query(self, version_id, query_region=None):
        # Get the plan used to run mwm_galactic_core_apogee.
        if self.name == "mwm_dust_core_apogee":
            self.mwm_galactic_carton = "mwm_galactic_core_apogee"
        elif self.name == "mwm_dust_core_dist_apogee":
            self.mwm_galactic_carton = "mwm_galactic_core_dist_apogee"
        else:
            raise ValueError(f"Invalid carton {self.name}.")

        if self.parameters:
            self.mwm_galactic_plan = self.parameters.get(
                self.mwm_galactic_carton + "_plan", self.plan
            )
        else:
            self.mwm_galactic_plan = self.plan

        gg_exists = (
            targetdb.Carton.select()
            .join(targetdb.Version)
            .where(
                targetdb.Carton.carton == self.mwm_galactic_carton,
                targetdb.Version.plan == self.mwm_galactic_plan,
                targetdb.Version.target_selection >> True,
            )
            .exists()
        )
        if not gg_exists:
            raise RuntimeError(self.mwm_galactic_carton + " has not been loaded yet.")

        fn = peewee.fn

        # GG quality flags
        ph_qual = TwoMassPSC.ph_qual
        cc_flg = TwoMassPSC.cc_flg
        rd_flg = TwoMassPSC.rd_flg
        rd_flag_1 = peewee.fn.substr(rd_flg, 2, 1).cast("integer")
        gal_contam = TwoMassPSC.gal_contam

        gallong = Gaia_DR3.l
        gallat = Gaia_DR3.b

        ipar = 1.0 / Gaia_DR3.parallax  # kpc
        zz = ipar * fn.sin(fn.radians(gallat))
        xx = ipar * fn.cos(fn.radians(gallat)) * fn.cos(fn.radians(gallong))
        yy = ipar * fn.cos(fn.radians(gallat)) * fn.sin(fn.radians(gallong))

        dist = fn.sqrt(fn.pow(xx, 2) + fn.pow(yy, 2) + fn.pow(zz, 2))

        aks_glimpse = 0.918 * (GLIMPSE.mag_h - GLIMPSE.mag4_5 - 0.08)
        aks_allwise = 0.918 * (AllWise.h_m_2mass - AllWise.w2mpro - 0.08)
        aks = fn.coalesce(aks_glimpse, aks_allwise)

        Ej_ks = 1.5 * aks

        j_ks_0_glimpse = GLIMPSE.mag_j - GLIMPSE.mag_ks - Ej_ks
        j_ks_0_allwise = AllWise.j_m_2mass - AllWise.k_m_2mass - Ej_ks
        j_ks_0 = fn.coalesce(j_ks_0_glimpse, j_ks_0_allwise)

        plxfracunc = 1 / Gaia_DR3.parallax_over_error
        dm = 5 * fn.log(1000.0 / Gaia_DR3.parallax / 10)
        absmag = TwoMassPSC.k_m - aks - dm

        query = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                peewee.Value(False).alias("selected"),  # Set selected to False
                Gaia_DR3.source_id.alias("gaia_souce_id"),
                Gaia_DR3.l.alias("gallong"),
                Gaia_DR3.b.alias("gallat"),
                Gaia_DR3.parallax,
                Gaia_DR3.phot_g_mean_mag,
                TwoMassPSC.h_m,
                TwoMassPSC.k_m,
                aks.alias("a_ks"),
                j_ks_0.alias("j_ks_0"),
            )
            .join(CatalogToGaia_DR3)
            .join(
                CatalogToAllWise,
                peewee.JOIN.LEFT_OUTER,
                on=(CatalogToAllWise.catalogid == CatalogToGaia_DR3.catalogid),
            )
            .join(AllWise, peewee.JOIN.LEFT_OUTER)
            .join_from(
                CatalogToGaia_DR3,
                CatalogToTwoMassPSC,
                peewee.JOIN.LEFT_OUTER,
                on=(CatalogToTwoMassPSC.catalogid == CatalogToGaia_DR3.catalogid),
            )
            .join(TwoMassPSC, peewee.JOIN.LEFT_OUTER)
            .join_from(
                CatalogToAllWise,
                CatalogToGLIMPSE,
                peewee.JOIN.LEFT_OUTER,
                on=(CatalogToGLIMPSE.catalogid == CatalogToAllWise.catalogid),
            )
            .join(GLIMPSE, peewee.JOIN.LEFT_OUTER)
            .where(
                ((CatalogToAllWise.version_id == version_id) & (CatalogToAllWise.best >> True))
                | (CatalogToAllWise.catalogid >> None)
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
            )
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
            )
            .where(
                ((CatalogToGLIMPSE.version_id == version_id) & (CatalogToGLIMPSE.best >> True))
                | (CatalogToGLIMPSE.catalogid >> None)
            )
            .where(
                TwoMassPSC.h_m < 11.2,
                fn.abs(zz) < 0.2,
                j_ks_0.is_null(False),
                j_ks_0 > 0.5,
                dist < 5,
                plxfracunc < 0.2,
                plxfracunc > 0,
                absmag < 2.6,
            )
            .where(
                ph_qual.regexp(".(A|B)."),
                gal_contam == 0,
                peewee.fn.substr(cc_flg, 2, 1) == "0",
                rd_flag_1 > 0,
                rd_flag_1 <= 3,
            )
            .order_by(CatalogToGaia_DR3.catalogid)
        )

        # The above order_by() ensures that the random selection in
        # post_process() gives the same result every time we run this carton.

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

    def post_process(self, model):
        """Select samples based on GG using Eddie Schlafly code.

        Goal subsample and to select only targets that are not already part
        of the Galactic Genesis sample. We use the real GG sample from
        targetdb.

        """

        data = pandas.read_sql(
            f"SELECT * FROM {model._meta.schema}.{model._meta.table_name};",
            self.database,
        )

        mwm_galactic = (
            targetdb.Target.select(targetdb.Target.catalogid)
            .join(targetdb.CartonToTarget)
            .join(targetdb.Carton)
            .join(targetdb.Version)
            .where(
                targetdb.Carton.carton == self.mwm_galactic_carton,
                targetdb.Version.plan == self.mwm_galactic_plan,
                targetdb.Version.target_selection >> True,
            )
        )

        gg_catids = tuple(zip(*mwm_galactic.tuples()))[0]
        gg_mask = numpy.in1d(data.catalogid, gg_catids)

        # subselect() is used below to randomly select a subset
        # from the results of the  query in build_query().
        # Hence, the query in build_query() uses order_by().
        # The order_by() in the query ensures that the random selection in
        # subselect() gives the same result every time we run this carton.
        #
        # Subsample based on GG and update DB.
        dust_gg_subsel = subselect(data, gg_mask)
        dust_gg_cid = peewee.ValuesList(
            zip(data.catalogid[dust_gg_subsel]),
            columns=("catalogid",),
            alias="vl",
        )

        (
            model.update({model.selected: True})
            .from_(dust_gg_cid)
            .where(model.catalogid == dust_gg_cid.c.catalogid)
            .execute()
        )

        return


# The below MWM_Dust_Core_Dist_apogee_Carton carton requires
# that the MWM_Galatic_Core_Dist_apogee_Carton
# has been run before this carton in the same target_selection_plan
# in target_selection.yml or in an earlier target_selection_plan.
# See target_selection.yml for example for mwm_dust_core_dist_apogee.
class MWM_Dust_Core_Dist_apogee_Carton(MWM_Dust_Core_apogee_Carton):
    """MWM Dust Core Dist apogee Carton.

    Same definition as the MWM_Dust_Core_apogee_Carton carton but
     subsampling using the mwm_galactic_core_dist_apogee.

    """

    name = "mwm_dust_core_dist_apogee"
