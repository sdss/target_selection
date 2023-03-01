#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2022-12-11
# @Filename: mwm_wd.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
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
