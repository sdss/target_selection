#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-31
# @Filename: mwm_ob.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
from peewee import fn
from scipy.special import erf

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2_RUWE,
                                             Gaia_DR2_TwoMass_Best_Neighbour,
                                             TIC_v8)

from .. import log
from . import BaseCarton


class MWM_OB_Carton(BaseCarton):
    """Milky Waky OB stars.

    Definition: Select all the hot, young stars in Gaia and 2MASS  with M_K < 0
    mag (M ~ 4 M_Sun) with Gaia G < 16 mag in the Milky Way, then subsampled
    to reduce the sample size / observing burden.

    Query on the Gaia archive:

    SELECT  g.source_id
    FROM gaiadr2.gaia_source as g
    JOIN gaiadr2.ruwe AS r
        USING (source_id)
    INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch
        ON g.source_id = xmatch.source_id
    INNER JOIN gaiadr1.tmass_original_valid AS tm
        ON tm.tmass_oid = xmatch.tmass_oid
    WHERE parallax < power(10.,(10.-tm.ks_m-0.)/5.)
        AND tm.j_m - tm.ks_m - 0.25*(g.phot_g_mean_mag - tm.ks_m) < 0.10
        AND tm.j_m - tm.ks_m - 0.25*(g.phot_g_mean_mag - tm.ks_m) > -0.30
        AND tm.j_m - tm.h_m < 0.15*(g.phot_g_mean_mag -tm.ks_m) + 0.05
        AND tm.j_m - tm.h_m > 0.15*(g.phot_g_mean_mag -tm.ks_m) - 0.15
        AND tm.j_m - tm.ks_m < 0.23*(g.phot_g_mean_mag -tm.ks_m) + 0.03
        AND g.phot_g_mean_mag > 2*(g.phot_g_mean_mag -tm.ks_m) + 3.0
        AND g.phot_g_mean_mag < 2*(g.phot_g_mean_mag -tm.ks_m) + 11.
        AND xmatch.angular_distance < 1. AND r.ruwe <1.4;

    Notes:
        - ks_m is the same as twomass_psc.h_m.

    """

    name = 'mwm_ob'
    category = 'science'

    def build_query(self, version_id, query_region=None):

        km = TIC_v8.kmag
        hm = TIC_v8.hmag
        jm = TIC_v8.jmag
        Gm = TIC_v8.gaiamag

        TMBN = Gaia_DR2_TwoMass_Best_Neighbour

        query = (Catalog.select(Catalog.catalogid,
                                Catalog.parallax,
                                km.alias('ks_m'))
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(Gaia_DR2_RUWE,
                       on=(TIC_v8.gaia_int == Gaia_DR2_RUWE.source_id))
                 .join(TMBN, on=(TMBN.source_id == Gaia_DR2_RUWE.source_id))
                 .where(TIC_v8.plx < fn.pow(10, ((10. - km) / 5.)),
                        jm - km - 0.25 * (Gm - km) < 0.10,
                        jm - km - 0.25 * (Gm - km) > -0.30,
                        jm - hm < 0.15 * (Gm - km) + 0.05,
                        jm - hm > 0.15 * (Gm - km) - 0.15,
                        jm - km < 0.23 * (Gm - km) + 0.03,
                        Gm > 2 * (Gm - km) + 3.0,
                        Gm < 2 * (Gm - km) + 11.,
                        TMBN.angular_distance < 1.,
                        Gaia_DR2_RUWE.ruwe < 1.4))

        if query_region:
            query = query.where(fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                    query_region[0],
                                                    query_region[1],
                                                    query_region[2]))

        return query

    def post_process(self, Model):
        """Subsample of main sample.

        See files  https://keeper.mpdl.mpg.de/f/a84318bb5213431ea513/?dl=1 and
        https://keeper.mpdl.mpg.de/f/ee2d70d8eb2f463f9562/?dl=1.

        """

        data = numpy.array(Model.select(Model.catalogid,
                                        Model.ks_m,
                                        Model.parallax),
                           dtype=[('catalogid', numpy.int64),
                                  ('ks_m', float),
                                  ('parallax', float)])

        log.debug(f'Number of initial rows: {len(data)}.')

        M_K_star = data['ks_m'] + 5. * numpy.log10(data['parallax'] / 100.)
        data_sel_idx = numpy.isfinite(M_K_star)
        data_sel = data[data_sel_idx]
        M_K_star = M_K_star[data_sel_idx]

        alpha = 0.4
        M_K = numpy.arange(min(M_K_star), 0., 0.01)
        p = erf(-alpha * (M_K - 0.1))

        catalogid_new = []

        for i in range(len(p) - 1):

            w = numpy.where((M_K_star < M_K[i + 1]) & (M_K_star >= M_K[i]))
            p_ = numpy.mean(p[i:i + 1])

            N = len(data[w])
            N_new = numpy.int(numpy.round(p_ * N, 2))
            catalogid_ = numpy.random.choice(data_sel['catalogid'][w], N_new,
                                             replace=False)

            catalogid_new += catalogid_

        log.debug(f'Number of selected rows: {len(catalogid_new)}.')

        with self.database.atomic():
            (Model
             .update({Model.selected: False})
             .where(Model.catalogid.not_in(catalogid_new)))
