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

from sdssdb.peewee.sdss5db.catalogdb import (GLIMPSE, AllWise, Catalog,
                                             CatalogToAllWise,
                                             CatalogToGLIMPSE, CatalogToTIC_v8,
                                             Gaia_DR2, TIC_v8, TwoMassPSC)

from . import BaseCarton


def lbp2xyz(ll, bb, pp):

    dist = 1. / pp  # kpc
    z = dist * numpy.sin(numpy.radians(bb))
    x = (dist * numpy.cos(numpy.radians(bb)) * numpy.cos(numpy.radians(ll)))
    y = (dist * numpy.cos(numpy.radians(bb)) * numpy.sin(numpy.radians(ll)))

    # z is definitely right; x and y depend on conventions.  As long
    # as this routine is used for x and y throughout this code, this should
    # not introduce any issues.

    return x, y, z


def subselect(data, othersel, downsampledby=1):
    """Turn x, y, z into pixel; bin up to see where we are missing."""

    ll = data['l']
    bb = data['b']
    xyz = lbp2xyz(ll, bb, data['parallax'])

    dustsel = numpy.ones(len(data), numpy.bool)

    ranges = [[-5, 5], [-5, 5], [-0.2, 0.2]]
    resolution = 0.1

    countshape = [numpy.round((r[1] - r[0]) / resolution).astype('i4') + 2
                  for r in ranges]

    coords = [numpy.clip(numpy.floor((c - r[0]) / resolution).astype('i4'),
                         -1, s - 2) + 1
              for (c, r, s) in zip(xyz, ranges, countshape)]
    centers = [r[0] + resolution * (numpy.arange((r[1] - r[0]) / resolution) + 0.5)
               for r in ranges]

    flatcoord = numpy.ravel_multi_index(coords, countshape)

    countsd = numpy.bincount(flatcoord[dustsel & ~othersel],
                             minlength=numpy.product(countshape))
    countso = numpy.bincount(flatcoord[othersel & dustsel],
                             minlength=numpy.product(countshape))
    countsd = countsd.reshape(countshape)
    countso = countso.reshape(countshape)

    pnumdensd = map_coordinates_wrap((countsd[1:-1, 1:-1, 1:-1], centers),
                                     [c[dustsel & ~othersel] for c in xyz],
                                     order=0)
    pnumdenso = map_coordinates_wrap((countso[1:-1, 1:-1, 1:-1], centers),
                                     [c[dustsel & ~othersel] for c in xyz],
                                     order=0)

    ntopick = 10. / downsampledby
    pnumdensleft = ntopick - pnumdenso

    numpy.random.seed(0)
    randnum = numpy.random.rand(numpy.sum(dustsel & ~othersel))

    subsel = (randnum < pnumdensleft / 1. / (pnumdensd + (pnumdensd == 0)))

    dustsel[othersel] = 0
    dustsel[dustsel & ~othersel] = subsel

    return dustsel


def map_coordinates_wrap(grid, coord, **kw):
    gridpts = grid[0]
    gridcoord = grid[1]
    outputcoordnorm = [numpy.interp(c, gc, numpy.arange(len(gc)))
                       for c, gc in zip(coord, gridcoord)]
    from scipy.ndimage import interpolation
    return interpolation.map_coordinates(gridpts, outputcoordnorm,
                                         cval=0., mode='constant', **kw)


def jkselect(data):
    jk0 = data['j_ks_0']
    return (jk0 > 0.7) & (data['h_m'] < 11)


def ghselect(data):
    gaiag = data['phot_g_mean_mag']
    return ((data['h_m'] < 11) & numpy.where(numpy.isfinite(gaiag),
                                             gaiag - data['h_m'] > 3.5,
                                             True))


class MWM_Dust_Carton(BaseCarton):
    """MWM Dust Carton.

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

    Non-SQL implementation:

        https://faun.rc.fas.harvard.edu/eschlafly/sdss5/dustsel.py

    """

    name = 'mwm_dust'
    mapper = 'MWM'
    category = 'science'
    program = 'Dust'

    def build_query(self, version_id, query_region=None):

        fn = peewee.fn

        twomass_cte = (TwoMassPSC
                       .select(TwoMassPSC.designation,
                               TwoMassPSC.k_m,
                               TwoMassPSC.h_m)
                       .where(TwoMassPSC.h_m < 11.2)
                       .cte('twomass_cte'))

        gallong = Gaia_DR2.l
        gallat = Gaia_DR2.b

        ipar = 1. / Gaia_DR2.parallax / 1000  # kpc
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

        plxfracunc = Gaia_DR2. parallax_error / Gaia_DR2. parallax
        dm = 5 * fn.log(1000. / Gaia_DR2.parallax / 10)
        absmag = twomass_cte.c.k_m - aks - dm

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR2.l, Gaia_DR2.b,
                         Gaia_DR2.parallax,
                         Gaia_DR2.phot_g_mean_mag,
                         twomass_cte.c.h_m,
                         aks.alias('a_ks'),
                         j_ks_0.alias('j_ks_0'),
                         peewee.Value(False).alias('dustghsubsel'),
                         peewee.Value(False).alias('dustjksubsel'))
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(twomass_cte,
                       on=(TIC_v8.twomass_psc == twomass_cte.c.designation))
                 .join_from(TIC_v8, Gaia_DR2)
                 .join_from(Catalog, CatalogToAllWise)
                 .join(AllWise)
                 .join_from(Catalog, CatalogToGLIMPSE)
                 .join(GLIMPSE)
                 .where(CatalogToAllWise.version_id == version_id,
                        CatalogToAllWise.best >> True)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(CatalogToGLIMPSE.version_id == version_id,
                        CatalogToGLIMPSE.best >> True)
                 .where(fn.abs(zz) < 0.2,
                        j_ks_0 > 0.5,
                        dist < 5,
                        plxfracunc < 0.2, plxfracunc > 0,
                        absmag < 2.6)
                 .with_cte(twomass_cte))

        if query_region:
            query = (query
                     .join_from(CatalogToAllWise, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """Select J-K and G-H samples using Eddie Schlafly code.

        Goal subsample and to select only targets that are not already part
        of the Galactic Genesis sample. Two possible conditions for GG are
        given based on (G-H) or (J-Ks) colours.

        """

        data = model.select().tuples()
        data = pandas.read_sql(
            f'SELECT * FROM {model._meta.schema}.{model._meta.table_name};',
            self.database)

        # Change selected to False for now.
        data.selected = False

        ghs = ghselect(data)
        jks = jkselect(data)

        dust_gh_subsel = subselect(data, ghs)
        dust_gh_cid = peewee.ValuesList(zip(data.catalogid[dust_gh_subsel]),
                                        columns=('catalogid',), alias='vl')

        # Subsample based on (G-H) and update DB.
        model.update({model.dustghsubsel: False}).execute()
        (model
         .update({model.dustghsubsel: True})
         .from_(dust_gh_cid)
         .where(model.catalogid == dust_gh_cid.c.catalogid)
         .execute())

        dust_jks_subsel = subselect(data, jks)
        dust_jks_cid = peewee.ValuesList(zip(data.catalogid[dust_jks_subsel]),
                                         columns=('catalogid',), alias='vl')

        # Subsample based on (J-Ks) and update DB.
        model.update({model.dustjksubsel: False}).execute()
        (model
         .update({model.dustjksubsel: True})
         .from_(dust_jks_cid)
         .where(model.catalogid == dust_jks_cid.c.catalogid)
         .execute())

        return
