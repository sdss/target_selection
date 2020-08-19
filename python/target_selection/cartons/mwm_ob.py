#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-31
# @Filename: mwm_ob.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             GAIA_ASSAS_SN_Cepheids, Gaia_DR2,
                                             Gaia_DR2_RUWE,
                                             Gaia_DR2_TwoMass_Best_Neighbour,
                                             TIC_v8, TwoMassPSC)

from . import BaseCarton


TMBN = Gaia_DR2_TwoMass_Best_Neighbour


class MWM_OB_Carton(BaseCarton):
    """Milky Way OB stars.

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
        WHERE parallax < power(10., (10. - tm.ks_m) / 5.)
            AND g.phot_g_mean_mag < 16.
            AND tm.j_m - tm.ks_m - 0.25 * (g.phot_g_mean_mag - tm.ks_m) < 0.10
            AND tm.j_m - tm.ks_m - 0.25 * (g.phot_g_mean_mag - tm.ks_m) > -0.30
            AND tm.j_m - tm.h_m < 0.15 * (g.phot_g_mean_mag - tm.ks_m) + 0.05
            AND tm.j_m - tm.h_m > 0.15 * (g.phot_g_mean_mag - tm.ks_m) - 0.15
            AND tm.j_m - tm.ks_m < 0.23 * (g.phot_g_mean_mag - tm.ks_m) + 0.03
            AND g.phot_g_mean_mag > 2 * (g.phot_g_mean_mag - tm.ks_m) + 3.0
            AND g.phot_g_mean_mag < 2 * (g.phot_g_mean_mag - tm.ks_m) + 11.
            AND xmatch.angular_distance < 1.
            AND r.ruwe < 1.4;

    Notes:
        - ks_m is the same as twomass_psc.h_m.

    """

    name = 'mwm_ob_core'
    mapper = 'MWM'
    category = 'science'
    cadence = 'mwm_ob_3x1'
    program = 'mwm_ob'
    priority = 2910

    def build_query(self, version_id, query_region=None):

        km = TwoMassPSC.k_m
        hm = TwoMassPSC.h_m
        jm = TwoMassPSC.j_m
        Gm = Gaia_DR2.phot_g_mean_mag

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.parallax,
                         Gaia_DR2.source_id.alias('gaia_source_id'),
                         km.alias('ks_m'))
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(Gaia_DR2_RUWE,
                       on=(Gaia_DR2.source_id == Gaia_DR2_RUWE.source_id))
                 .join(TMBN,
                       on=(TMBN.source_id == Gaia_DR2_RUWE.source_id))
                 .join(TwoMassPSC,
                       on=(TMBN.tmass_pts_key == TwoMassPSC.pts_key))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(Gaia_DR2.parallax < fn.pow(10, ((10. - km - 0.61) / 5.)),
                        Gm < 16.,
                        Gm > 12,
                        jm - km - 0.25 * (Gm - km) < 0.10,
                        jm - km - 0.25 * (Gm - km) > -0.30,
                        jm - hm < 0.15 * (Gm - km) + 0.05,
                        jm - hm > 0.15 * (Gm - km) - 0.15,
                        jm - km < 0.23 * (Gm - km) + 0.03,
                        Gm > 2 * (Gm - km) + 3.0,
                        TMBN.angular_distance < 1.))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, Model):
        """Subsample of main sample.

        See files  https://keeper.mpdl.mpg.de/f/a84318bb5213431ea513/?dl=1 and
        https://keeper.mpdl.mpg.de/f/ee2d70d8eb2f463f9562/?dl=1.

        After the gamma-3 sims it was decided that the subsampling was not
        needed but we leave it here just in case.

        """

        return

        # data = numpy.array(Model.select(Model.catalogid,
        #                                 Model.ks_m,
        #                                 Model.parallax).tuples(),
        #                    dtype=[('catalogid', numpy.int64),
        #                           ('ks_m', float),
        #                           ('parallax', float)])

        # self.log.debug(f'Number of initial rows: {len(data)}.')

        # data = data[data['parallax'] > 0]
        # M_K_star = data['ks_m'] + 5. * numpy.log10(data['parallax'] / 100.)
        # data_sel_idx = numpy.isfinite(M_K_star)
        # data_sel = data[data_sel_idx]
        # M_K_star = M_K_star[data_sel_idx]

        # alpha = 0.4
        # M_K = numpy.arange(min(M_K_star), 0., 0.01)
        # p = erf(-alpha * (M_K - 0.1))

        # catalogid_new = []

        # # p ~ 3000 so this for loop is actually not that inefficient.
        # for i in range(len(p) - 1):

        #     w = numpy.where((M_K_star < M_K[i + 1]) & (M_K_star >= M_K[i]))
        #     p_ = numpy.mean(p[i:i + 1])

        #     N = len(data[w])
        #     N_new = numpy.int(numpy.round(p_ * N, 2))
        #     catalogid_ = numpy.random.choice(data_sel['catalogid'][w], N_new,
        #                                      replace=False)
        #     catalogid_new += catalogid_.tolist()

        # catalogid_new = set(catalogid_new)

        # self.log.debug(f'Number of selected rows: {len(catalogid_new)}.')

        # self.log.debug('Applying selected mask.')

        # values = ValuesList(zip(catalogid_new), columns=('catalogid',),
        #                     alias='vl')

        # with self.database.atomic():
        #     # Change everything to selected=False
        #     (Model.update({Model.selected: False}).execute())
        #     # Select the catalogids we calculated. If we tried to do a
        #     # .where(Model.catalogid != values.c.catalogid) that would take
        #     # forever, not sure why.
        #     (Model
        #      .update({Model.selected: True})
        #      .from_(values)
        #      .where(Model.catalogid == values.c.catalogid)
        #      .execute())


# class MWM_OB_MC_Carton(BaseCarton):
#     """Magellanic Clouds OB stars.

#     Definition: Select all the hot, young stars in Gaia and 2MASS with M_K < 0
#     mag (M ~ 4 M_Sun) with Gaia G < 16 mag in the Magellanic Clouds.

#     Pseudo-query:

#         phot_g_mean_ mag < 16.0
#         AND parallax < power(10., (10. -  ks_m) / 5.)
#         AND jm - ks_m - 0.25 * (phot_g_mean_mag - ks_m) < 0.10
#         AND j_m - h_m < 0.15 * (phot_g_mean_mag -ks_m) + 0.05
#         AND j_m - ks_m < 0.23 * (phot_g_mean_mag -ks_m) + 0.03
#         AND phot_g_mean_mag > 2 * (phot_g_mean_mag -ks_m) + 3.0
#         AND  265.0 < l < 310.0
#         AND -50 < b < -25.0
#         AND xmatch.angular distance < 1.0

#     Notes:
#         - ks_m is the same as twomass_psc.h_m.

#     """

#     name = 'mwm_ob_mc'
#     mapper = 'MWM'
#     category = 'science'
#     cadence = 'mwm_ob_3x1'
#     program = 'OB'

#     def build_query(self, version_id, query_region=None):

#         b = Gaia_DR2.b
#         l = Gaia_DR2.l  # noqa

#         km = TwoMassPSC.k_m
#         hm = TwoMassPSC.h_m
#         jm = TwoMassPSC.j_m
#         Gm = Gaia_DR2.phot_g_mean_mag

#         query = (CatalogToTIC_v8
#                  .select(CatalogToTIC_v8.catalogid,
#                          Gaia_DR2.source_id.alias('gaia_source_id'))
#                  .join(TIC_v8)
#                  .join(Gaia_DR2)
#                  .join(TMBN, on=(TMBN.source_id == Gaia_DR2.source_id))
#                  .join(TwoMassPSC,
#                        on=(TMBN.tmass_pts_key == TwoMassPSC.pts_key))
#                  .where(CatalogToTIC_v8.version_id == version_id,
#                         CatalogToTIC_v8.best >> True)
#                  .where(Gaia_DR2.parallax < fn.pow(10, ((10. - km) / 5.)),
#                         Gm < 16.,
#                         Gm > 12,
#                         jm - km - 0.25 * (Gm - km) < 0.10,
#                         jm - hm < 0.15 * (Gm - km) + 0.05,
#                         jm - km < 0.23 * (Gm - km) + 0.03,
#                         Gm > 2 * (Gm - km) + 3.0,
#                         TMBN.angular_distance < 1.,
#                         b > -50., b < -25., l > 265., l < 310.))

#         if query_region:
#             query = (query
#                      .join_from(CatalogToTIC_v8, Catalog)
#                      .where(peewee.fn.q3c_radial_query(Catalog.ra,
#                                                        Catalog.dec,
#                                                        query_region[0],
#                                                        query_region[1],
#                                                        query_region[2])))

#         return query


class MWM_OB_Cepheids_Carton(BaseCarton):
    """Milky Way Cepheids.

    Definition: List of Cepheids compiled by Inno et al. (in prep). The
    catalogue is obtained by using Gaia DR2 and ASAS-SN. Gaia parallaxes and
    variability information are used to select an initial sample, for which
    multi-epoch V-band photometry from the ASAS-SN survey is retrieved. From
    the ASAS-SN lightcurves, periods and Fourier parameters are derived. The
    Fourier parameters are used to identify classical Cepheids. The final
    catalogue is assembled by including several public classic Cepheid
    databases: OGLE-cCs, Gaia-cCs, ASASSN-cCs, VSX-cCs, Simbad-cCs,  and
    WISE-cCs. The final catalogue consists of ~3000 targets, between
    6 < G < 17 mag.

    """

    name = 'mwm_ob_cepheids'
    mapper = 'MWM'
    category = 'science'
    cadence = 'mwm_ob_3x1'
    program = 'mwm_ob'
    priority = 2910

    def build_query(self, version_id, query_region=None):

        query = (GAIA_ASSAS_SN_Cepheids
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id)
                 .join(Gaia_DR2)
                 .join(TIC_v8)
                 .join(CatalogToTIC_v8)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.phot_g_mean_mag > 12))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
