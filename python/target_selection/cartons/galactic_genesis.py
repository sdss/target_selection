#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-29
# @Filename: galactic_genesis.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import Catalog, CatalogToTIC_v8, TIC_v8

from . import BaseCarton


class GalacticGenesisCarton(BaseCarton):
    """Galactic Genesis carton.

    Definition: Selection of all IR-bright, red stars – vast majority are red
    giants; follows the density distribution of the MW (concentrated in the
    plane and bulge). Select sources brighter than H<11 AND ( (G-H) > 3.5 OR
    Gaia non-detection). Approximately 5 million stars.

    """

    name = 'galactic_genesis'
    category = 'science'

    def build_query(self, version_id, query_region=None):

        gg = (TIC_v8
              .select(Catalog.catalogid)
              .join(CatalogToTIC_v8)
              .join(Catalog)
              .where(TIC_v8.hmag < self.parameters['h_max'],
                     (((TIC_v8.hmag - TIC_v8.gmag) > self.parameters['h_g']) |
                      TIC_v8.gaia >> None),
                     Catalog.version_id == version_id))

        if query_region:
            gg = gg.where(peewee.fn.q3c_radial_query(Catalog.ra, Catalog.dec,
                                                     query_region[0],
                                                     query_region[1],
                                                     query_region[2]))

        return gg
