#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-02
# @Filename: guide.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db import catalogdb, SDSS5dbModel


catalogdb.database.set_profile('mako-tunnel')


class Guide(SDSS5dbModel):

    source_id = peewee.BigIntegerField()

    class Meta:
        db_table = 'guide'
        schema = 'catalogdb'


class TempSample(SDSS5dbModel):

    source_id = peewee.BigIntegerField()
    ra = peewee.DoubleField()
    dec = peewee.DoubleField()

    class Meta:
        db_table = 'temp_sample'


if catalogdb.database.table_exists('guide', 'catalogdb'):
    peewee.SchemaManager(Guide).truncate_table()
else:
    catalogdb.database.create_tables([Guide])

# catalogdb.database.create_tables([TempSample], temporary=True)


sample_data = (catalogdb.GaiaDR2Source
               .select(catalogdb.GaiaDR2Source.source_id,
                       catalogdb.GaiaDR2Source.ra,
                       catalogdb.GaiaDR2Source.dec)
               .join(catalogdb.GaiaDR2Clean)
               .where((catalogdb.GaiaDR2Source.phot_g_mean_mag > 13) &
                      (catalogdb.GaiaDR2Source.phot_g_mean_mag < 17))
            #    .where(catalogdb.GaiaDR2Source.cone_search(10, 10, 10))
               .cte('sample_data'))

# TempSample.insert_from(sample_data.select(), [TempSample.source_id,
#                                               TempSample.ra, TempSample.dec])

# catalogdb.database.execute_sql('CREATE INDEX temp_sample_q3c_idx '
#                                'ON temp_sample (q3c_ang2ipix(ra, dec));')
# catalogdb.database.execute_sql('CLUSTER temp_sample_q3c_idx on temp_sample;')

print('a')

# Get second nearest neighbour (because we are self-matching). In q3c_join the
# order must be Guide and then guide_alias. Not sure why but otherwise it takes forever.

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
print('b')
# Left lateral join with the subquery. The rows that don't have a nearest neighbour within the
# joining distance will have subq.pk == NULL. We select those.

decollided = (sample_data
              .select(sample_data.c.source_id)
              .join(subq, peewee.JOIN.LEFT_LATERAL, on=peewee.SQL('true'))
              .where(subq.c.source_id.is_null())
               .where(peewee.fn.q3c_poly_query(
                   sample_data.c.ra,
                   sample_data.c.dec,
                   peewee.SQL('\'((0, 0), (0, 20), (20, 20), (20, 0))\'::polygon')))

              .with_cte(sample_data))
print('c')
print(decollided)
Guide.insert_from(decollided, [Guide.source_id]).execute()

print(Guide.select().count())
