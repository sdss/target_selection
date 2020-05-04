#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-26
# @Filename: guide.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             Gaia_DR2, Gaia_DR2_Clean, TIC_v8)

from .base import BaseCarton


class GuideCarton(BaseCarton):

    name = 'guide'
    category = 'guide'

    tile = False
    tile_region = None
    tile_num = 501

    def build_query(self, version_id):

        sample_data = (Catalog
                       .select(Catalog.catalogid,
                               Catalog.ra,
                               Catalog.dec,
                               Gaia_DR2.phot_g_mean_mag)
                       .join(CatalogToTIC_v8)
                       .join(TIC_v8)
                       .join(Gaia_DR2)
                       .join(Gaia_DR2_Clean)
                       .where((Gaia_DR2.phot_g_mean_mag > self.config['g_min']) &
                              (Gaia_DR2.phot_g_mean_mag < self.config['g_max']))
                       .where(Catalog.version_id == version_id,
                              CatalogToTIC_v8.version_id == version_id)
                       .cte('sample_data'))

        # We should use q3c_join_pm and q3c_dist_pm here.
        subq = (Gaia_DR2
                .select(Gaia_DR2.source_id)
                .where(peewee.fn.q3c_join(sample_data.c.ra, sample_data.c.dec,
                                          Gaia_DR2.ra,
                                          Gaia_DR2.dec,
                                          self.config['min_separation']))
                .order_by(peewee.fn.q3c_dist(sample_data.c.ra,
                                             sample_data.c.dec,
                                             Gaia_DR2.ra,
                                             Gaia_DR2.dec).asc())
                .limit(1)
                .offset(1))

        decollided = (sample_data
                      .select(sample_data.c.catalogid,
                              sample_data.c.phot_g_mean_mag.alias('magnitude_g'))
                      .join(subq, peewee.JOIN.LEFT_LATERAL, on=peewee.SQL('true'))
                      .where(subq.c.source_id.is_null())
                      .with_cte(sample_data))

        return decollided
