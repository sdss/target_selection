#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_halo.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (BestBrightest, Catalog,
                                             CatalogToAllWise, CatalogToTIC_v8,
                                             SkyMapperGaia, TIC_v8)

from . import BaseCarton


class MWM_Halo_Best_Brightest_Carton(BaseCarton):
    """MWM Halo Best & Brightest.

    Definition:
        Photometrically selected metal-poor stars to (hopefully) be
        observed with spare fibers. Selected using "v1" and "v2"
        from Schlaufman + Casey 2014.

        A sample of metal-poor giants selected based on optical and
        infrared photometry.

        Select all stars from B&B catalog in two rounds. First priority
        with version=2, second priority by version=1.

        Original selection criteria:

            0.45 < J-H < 0.6
            W3 > 8
            -0.04 < W1-W2 < 0.04
            J - W2 > 0.5
            e^z / (1 + e^z) > 0.13 where z is defined in Schlaufman+Casey
            (version=1): J-W2 > 0.5 (B-V -0.8) + 0.6
            (version=2): B-V < 1.2

            (250k sources for version=2, 396k sources for version=1)

    """

    name = 'mwm_halo_bb'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_filler'
    cadence = None
    priority = None

    def build_query(self, version_id, query_region=None):

        query = (CatalogToAllWise
                 .select(CatalogToAllWise.catalogid,
                         BestBrightest.designation,
                         BestBrightest.version)
                 .join(BestBrightest,
                       on=(BestBrightest.cntr == CatalogToAllWise.target_id))
                 .where(CatalogToAllWise.version_id == version_id,
                        CatalogToAllWise.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToAllWise, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """Define priority based on version."""

        for version in [1, 2]:
            priority = self.parameters[f'priority_version_{version}']
            (model
             .update({model.priority: priority})
             .where(model.version == version).execute())


class MWM_Halo_SkyMapper_Carton(BaseCarton):
    """MWM Halo SkyMapper.

    Definition: Photometrically selected metal-poor stars with [Fe/H] < -1.7
    (determined using SkyMapper Ca K photometry).

    """

    name = 'mwm_halo_sm'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_filler'
    cadence = None
    priority = 6050

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         SkyMapperGaia.gaia_source_id)
                 .join(TIC_v8)
                 .join(SkyMapperGaia,
                       on=(TIC_v8.gaia_int == SkyMapperGaia.gaia_source_id))
                 .where(SkyMapperGaia.feh < self.parameters['feh'],
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
