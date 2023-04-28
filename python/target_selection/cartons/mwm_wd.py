#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2022-12-11
# @Filename: mwm_wd.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             Gaia_DR3, Gedr3spur_main,
                                             WD_gaia_dr3)

from target_selection.cartons import BaseCarton


class MWM_WD_PWD_Carton(BaseCarton):
    """MWM White Dwarfs.

    Definition:

    Reference: Gentile Fusillo et al. 2021
    A catalogue of white dwarfs in Gaia EDR3.

    select all targets from table wd_gaia_d3
    where Pwd > 0.5 and Gmag <= 20.

    (above Gmag means gmag_vega)

    """

    name = 'mwm_wd_pwd'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_wd'
    instrument = 'BOSS'
    cadence = 'dark_2x1'
    priority = 1400
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         WD_gaia_dr3.gaiaedr3,
                         WD_gaia_dr3.pwd,
                         WD_gaia_dr3.gmag_vega)
                 .join(WD_gaia_dr3,
                       on=(CatalogToGaia_DR3.target_id == WD_gaia_dr3.gaiaedr3))
                 .where(WD_gaia_dr3.pwd > 0.5,
                        WD_gaia_dr3.gmag_vega <= 20.0,
                        CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_WD_Gaia_Carton(BaseCarton):
    """MWM Gaia White Dwarfs.

    Definition:

    bp_rp=gaiadr3.phot_bp_mean_mag-gaiadr3.phot_rp_mean_mag
    par_log = log10(gaiadr3.parallax/100.)

    First select:
    gaiadr3.phot_g_mean_mag < 20.0
    and gaiadr3.parallax*power(10.,0.2*(22.-gaiadr3.phot_g_mean_mag))> 3.
    and phot_g_mean_mag+5*par_log-0.5-5.5*(bp_rp+0.5) > 0.
    and ((bp_rp<0.4) or ((gaiadr3.phot_g_mean_mag+5*par_log)> 11+2.2*(bp_rp-0.5) and (bp_rp >0.4))

    Then select:
    gaiadr3.phot_g_mean_mag+5*par_log - 6.7-2.2*bp_rp  > 0. or bp_rp < -0.6

    And finally
    gDR3spur_main.fidelity_v2 > 0.9

    """

    name = 'mwm_wd_gaia'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_wd'
    instrument = 'BOSS'
    cadence = 'dark_2x1'
    priority = 1400
    can_offset = True

    def build_query(self, version_id, query_region=None):

        bp_rp = Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag
        gaiaG = Gaia_DR3.phot_g_mean_mag
        par_log = fn.log10(Gaia_DR3.parallax / 100.)

        gaia_cut = ((gaiaG < 20.0) &
                    (Gaia_DR3.parallax * fn.power(10., 0.2 * (22. - gaiaG)) > 3.) &
                    (gaiaG + 5 * par_log - 0.5 - 5.5 * (bp_rp + 0.5) > 0.) &
                    ((bp_rp < 0.4) | ((bp_rp > 0.4) &
                                      ((gaiaG + 5 * par_log) > (11 + 2.2 * (bp_rp - 0.5))))))

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.ra,
                         Gaia_DR3.dec,
                         bp_rp,
                         gaiaG,
                         Gaia_DR3.parallax,
                         Gedr3spur_main.fidelity_v2)
                 .join(Gaia_DR3)
                 .join(Gedr3spur_main, on=(Gaia_DR3.source_id == Gedr3spur_main.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        gaia_cut,
                        ((gaiaG + 5 * par_log - 6.7 - 2.2 * bp_rp) > 0.) | (bp_rp < -0.6),
                        Gedr3spur_main.fidelity_v2 > 0.9))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
