#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-02
# @Filename: mwm_planet.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (TESS_TOI, Catalog,
                                             CatalogToTIC_v8, TIC_v8,
                                             TwoMassPSC)

from . import BaseCarton


class MWM_TESS_Planets_Carton(BaseCarton):
    """TESS Planets.

    Definition:
        All planet-candidate TOIs + Community TOIs (CTIOs), TESS
        Candidate-Target List (CTL) selected by a TIC Priority cut.
        All cut between 7 < H < 12.

    All the candidates are compiled in catalogdb.tess_toi. Note that this list
    contains duplicates.

    """

    name = 'mwm_planet_tess'
    program = 'mwm_planet'
    category = 'science'
    mapper = 'MWM'
    cadence = None
    priority = 2600

    def build_query(self, version_id, query_region=None):

        query = (TESS_TOI
                 .select(CatalogToTIC_v8.catalogid,
                         TESS_TOI.ticid,
                         TESS_TOI.tess_disposition,
                         TwoMassPSC.h_m,
                         TESS_TOI.target_type)
                 .join(TIC_v8)
                 .join(TwoMassPSC)
                 .join_from(TIC_v8, CatalogToTIC_v8)
                 .where(TwoMassPSC.h_m > self.parameters['h_min'],
                        TwoMassPSC.h_m < self.parameters['h_max'],
                        CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .distinct([TESS_TOI.ticid]))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model, **kwargs):
        """Set priorities based on ``target_type``.

        exo_TOI: 2600
        exo_CTOI: 2605
        2min: 2010

        """

        priorities = {'exo_TOI': 2600,
                      'exo_CTOI': 2605,
                      '2min': 2610}

        for target_type, priority in priorities.items():
            (model
             .update({model.priority: priority})
             .where(model.target_type == target_type).execute())
