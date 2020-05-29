#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-29
# @Filename: galactic_genesis.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             TIC_v8, TwoMassPSC)

from . import BaseCarton


class GalacticGenesisCarton(BaseCarton):

    name = 'galactic_genesis'
    category = 'science'

    def build_query(self, version_id):

        gg = (TwoMassPSC
              .select(Catalog.catalogid)
              .join(TIC_v8, 'LEFT OUTER')
              .join(CatalogToTIC_v8)
              .join(Catalog)
              .where(TwoMassPSC.h_m < self.parameters['h_max'],
                     ((TwoMassPSC.h_m - TIC_v8.gmag) > self.parameters['h_g'] |
                      TIC_v8.gaia >> None),
                     Catalog.version_id == version_id))

        return gg
