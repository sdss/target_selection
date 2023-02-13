#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-31
# @Filename: mwm_ob.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             CatalogToTIC_v8,
                                             CatalogToTwoMassPSC,
                                             GAIA_ASSAS_SN_Cepheids, Gaia_DR2,
                                             Gaia_DR2_TwoMass_Best_Neighbour,
                                             Gaia_DR3,
                                             Gaia_dr3_astrophysical_parameters,
                                             TIC_v8, TwoMassPSC)

from target_selection.cartons import BaseCarton


TMBN = Gaia_DR2_TwoMass_Best_Neighbour


class MWM_OB_Carton(BaseCarton):
    """Milky Way OB stars.

    Definition: Select all the hot, young stars in Gaia and 2MASS  with M_K < 0
    mag (M ~ 4 M_Sun) with Gaia G < 16 mag in the Milky Way, then subsampled
    to reduce the sample size / observing burden.

    Query on the Gaia archive (this is the v0.5 query, v1 is slightly modified):

        SELECT  g.source_id
        FROM gaiadr3.gaia_source as g
        JOIN gaiadr3.ruwe AS r
            USING (source_id)
        INNER JOIN gaiadr3.tmass_best_neighbour AS xmatch
            ON g.source_id = xmatch.source_id
        INNER JOIN gaiadr1.tmass_original_valid AS tm
            ON tm.tmass_oid = xmatch.tmass_oid
        WHERE parallax < fn.pow(10, ((10. - tm.ks_m - 0.61) / 5.))
            AND g.phot_g_mean_mag < 16.
            AND tm.j_m - tm.ks_m - 0.25 * (g.phot_g_mean_mag - tm.ks_m) < 0.10
            AND tm.j_m - tm.ks_m - 0.25 * (g.phot_g_mean_mag - tm.ks_m) > -0.30
            AND tm.j_m - tm.h_m < 0.15 * (g.phot_g_mean_mag - tm.ks_m) + 0.05
            AND tm.j_m - tm.h_m > 0.15 * (g.phot_g_mean_mag - tm.ks_m) - 0.15
            AND tm.j_m - tm.ks_m < 0.23 * (g.phot_g_mean_mag - tm.ks_m) + 0.03
            AND g.phot_g_mean_mag > 2 * (g.phot_g_mean_mag - tm.ks_m) + 3.0

    Notes:
        - ks_m is the same as twomass_psc.h_m.

    """

    name = 'mwm_ob_core'
    mapper = 'MWM'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'bright_3x1'
    program = 'mwm_ob'
    priority = None
    can_offset = True

    def build_query(self, version_id, query_region=None):

        km = TwoMassPSC.k_m
        hm = TwoMassPSC.h_m
        jm = TwoMassPSC.j_m
        Gm = Gaia_DR3.phot_g_mean_mag

        query = (Gaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.ra,
                         Gaia_DR3.dec,
                         km,
                         hm,
                         jm,
                         Gm,
                         Gaia_dr3_astrophysical_parameters.teff_esphs,
                         Gaia_DR3.parallax,
                         Gaia_DR3.source_id,
                         km.alias('ks_m'))
                 .join(CatalogToGaia_DR3)
                 .join(CatalogToTwoMassPSC,
                       on=(CatalogToTwoMassPSC.catalogid == CatalogToGaia_DR3.catalogid))
                 .join(TwoMassPSC)
                 .join_from(Gaia_DR3, Gaia_dr3_astrophysical_parameters,
                            on=(Gaia_dr3_astrophysical_parameters.source_id == Gaia_DR3.source_id))
                 .where(CatalogToTwoMassPSC.version_id == version_id,
                        CatalogToTwoMassPSC.best >> True,
                        CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True)
                 .where(Gaia_DR3.parallax < fn.pow(10, ((10. - km - 0.61) / 5.)),
                        Gm < 16.,
                        jm - km - 0.25 * (Gm - km) < 0.10,
                        jm - km - 0.25 * (Gm - km) > -0.30,
                        jm - hm < 0.15 * (Gm - km) + 0.05,
                        jm - hm > 0.15 * (Gm - km) - 0.15,
                        jm - km < 0.23 * (Gm - km) + 0.03,
                        Gm > 2 * (Gm - km) + 3.0))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """Adjust priorities based on ``Teff_esphs``."""

        (model.update({model.priority: 2700})
         .where(model.teff_esphs > 10000)).execute()

        (model.update({model.priority: 2910})
         .where((model.teff_esphs < 10000) | (model.teff_esphs.is_null()))).execute()


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
    instrument = 'BOSS'
    cadence = 'bright_3x1'
    program = 'mwm_ob'
    priority = 2910
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (GAIA_ASSAS_SN_Cepheids
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id)
                 .join(Gaia_DR2)
                 .join(TIC_v8)
                 .join(CatalogToTIC_v8)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
