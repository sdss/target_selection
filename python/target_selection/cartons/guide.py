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


class Guide(BaseCarton):

    name = 'guide'
    category = 'guide'

    def build_query(self):

        sample_data = (catalogdb.GaiaDR2Source
                       .select(catalogdb.GaiaDR2Source.source_id,
                               catalogdb.GaiaDR2Source.ra,
                               catalogdb.GaiaDR2Source.dec,
                               catalogdb.GaiaDR2Source.phot_g_mean_mag)
                       .join(catalogdb.GaiaDR2Clean)
                       .where((catalogdb.GaiaDR2Source.phot_g_mean_mag > 13) &
                              (catalogdb.GaiaDR2Source.phot_g_mean_mag < 17))
                       .cte('sample_data'))

        subq = (catalogdb.GaiaDR2Source
                .select(catalogdb.GaiaDR2Source.source_id)
                .where(catalogdb.GaiaDR2Source.cone_search(sample_data.c.ra,
                                                           sample_data.c.dec,
                                                           1 / 60.))
                .where(peewee.fn.q3c_join(sample_data.c.ra, sample_data.c.dec,
                                          catalogdb.GaiaDR2Source.ra,
                                          catalogdb.GaiaDR2Source.dec,
                                          1. / 60))
                .order_by(peewee.fn.q3c_dist(sample_data.c.ra,
                                             sample_data.c.dec,
                                             catalogdb.GaiaDR2Source.ra,
                                             catalogdb.GaiaDR2Source.dec).asc())
                .limit(1)
                .offset(1))

        decollided = (sample_data
                      .select(sample_data.c.source_id,
                              sample_data.c.phot_g_mean_mag.alias('magnitude_g'))
                      .join(subq, peewee.JOIN.LEFT_LATERAL, on=peewee.SQL('true'))
                      .where(subq.c.source_id.is_null())
                      .with_cte(sample_data))

        return decollided
