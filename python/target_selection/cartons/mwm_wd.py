#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_wd.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2_WD, TIC_v8)

from . import BaseCarton


class MWM_WD_Carton(BaseCarton):
    """MWM White Dwarfs.

    Definition: all targets from Gentile Fusillo et al. 2019 (table
    gaia_dr2_wd) where Pwd > 0.5 and Gmag <= 20.

    """

    name = 'mwm_wd_core'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_wd'
    cadence = 'mwm_wd_2x1'
    priority = 1400

    def build_query(self, version_id, query_region=None):

        query = (Gaia_DR2_WD
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2_WD.source_id,
                         Gaia_DR2_WD.pwd,
                         Gaia_DR2_WD.gmag)
                 .join(TIC_v8, on=(TIC_v8.gaia_int == Gaia_DR2_WD.source_id))
                 .join(CatalogToTIC_v8)
                 .where(Gaia_DR2_WD.pwd > self.parameters['pwd'],
                        Gaia_DR2_WD.g_gaia_mag <= self.parameters['gaia_mag'],
                        CatalogToTIC_v8.version_id == version_id,
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
