#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-26
# @Filename: guide.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db import catalogdb

from .base import BaseCarton


class GuideCarton(BaseCarton):

    name = 'guide'
    category = 'guide'

    tile = True
    tile_region = None
    tile_num = 501

    def build_query(self):

        sample_data = (catalogdb.Gaia_DR2
                       .select(catalogdb.Gaia_DR2.source_id,
                               catalogdb.Gaia_DR2.ra,
                               catalogdb.Gaia_DR2.dec,
                               catalogdb.Gaia_DR2.phot_g_mean_mag,
                               catalogdb.Gaia_DR2.ref_epoch)
                       .join(catalogdb.GaiaDR2Clean)
                       .where((catalogdb.Gaia_DR2.phot_g_mean_mag > self.config['g_min']) &
                              (catalogdb.Gaia_DR2.phot_g_mean_mag < self.config['g_max']))
                       .cte('sample_data'))

        # We should use q3c_join_pm and q3c_dist_pm here.
        subq = (catalogdb.Gaia_DR2
                .select(catalogdb.Gaia_DR2.source_id)
                .where(catalogdb.Gaia_DR2.cone_search(sample_data.c.ra,
                                                      sample_data.c.dec,
                                                      self.config['min_separation']))
                .where(peewee.fn.q3c_join(sample_data.c.ra, sample_data.c.dec,
                                          catalogdb.Gaia_DR2.ra,
                                          catalogdb.Gaia_DR2.dec,
                                          self.config['min_separation']))
                .order_by(peewee.fn.q3c_dist(sample_data.c.ra,
                                             sample_data.c.dec,
                                             catalogdb.Gaia_DR2.ra,
                                             catalogdb.Gaia_DR2.dec).asc())
                .limit(1)
                .offset(1))

        decollided = (sample_data
                      .select(sample_data.c.source_id,
                              sample_data.c.phot_g_mean_mag.alias('magnitude_g'),
                              sample_data.c.ref_epoch.alias('epoch'))
                      .join(subq, peewee.JOIN.LEFT_LATERAL, on=peewee.SQL('true'))
                      .where(subq.c.source_id.is_null())
                      .with_cte(sample_data))

        return decollided
