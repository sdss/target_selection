#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2023-04-25
# @Filename: mwm_bin_vis.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import JOIN

from sdssdb.peewee.sdss5db.catalogdb import (CatalogToGaia_DR3,
                                             CatalogToTwoMassPSC, Gaia_DR3,
                                             Galah_dr3, Lamost_dr6,
                                             RAVE_DR6_Gaia_DR3_XMatch,
                                             RAVE_DR6_Gauguin_Madera,
                                             SDSS_DR17_APOGEE_Allstarmerge,
                                             TwoMassPSC,
                                             Visual_binary_gaia_dr3)

from target_selection.cartons import BaseCarton


class MWM_Wide_Binaries_Base_Carton(BaseCarton):
    """MWM wide binaries carton."""

    mapper = 'MWM'
    category = 'science'
    program = 'mwm_filler'
    cadence = 'bright_1x1'
    priority = None

    def build_query(self, version_id, query_region=None):

        query = None

        for col in ['source_id1', 'source_id2']:
            source_col = getattr(Visual_binary_gaia_dr3, col)
            this_query = (CatalogToGaia_DR3
                          .select(CatalogToGaia_DR3.catalogid,
                                  Gaia_DR3.source_id,
                                  Visual_binary_gaia_dr3.source_id1,
                                  Visual_binary_gaia_dr3.source_id2,
                                  TwoMassPSC.h_m,
                                  SDSS_DR17_APOGEE_Allstarmerge.apogee_id,
                                  Galah_dr3.dr3_source_id.alias('galah_dr3_source_id'),
                                  Lamost_dr6.source_id.alias('lamost_dr6_source_id'),
                                  RAVE_DR6_Gauguin_Madera.rave_obs_id,
                                  Visual_binary_gaia_dr3.dr2_parallax1,
                                  Visual_binary_gaia_dr3.dr2_parallax2,
                                  Visual_binary_gaia_dr3.sep_au)
                          .join(Gaia_DR3)
                          .join(Visual_binary_gaia_dr3,
                                on=(Gaia_DR3.source_id == source_col))
                          .join_from(CatalogToGaia_DR3,
                                     CatalogToTwoMassPSC,
                                     on=(CatalogToGaia_DR3.catalogid ==
                                         CatalogToTwoMassPSC.catalogid))
                          .join(TwoMassPSC,
                                on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
                          .join(SDSS_DR17_APOGEE_Allstarmerge,
                                on=(TwoMassPSC.designation == peewee.fn.replace(
                                    SDSS_DR17_APOGEE_Allstarmerge.apogee_id, '2M', '')),
                                join_type=JOIN.LEFT_OUTER)
                          .join_from(Visual_binary_gaia_dr3, Galah_dr3,
                                     on=(source_col == Galah_dr3.dr3_source_id),
                                     join_type=JOIN.LEFT_OUTER)
                          .join_from(Visual_binary_gaia_dr3, Lamost_dr6,
                                     on=(source_col == Lamost_dr6.source_id),
                                     join_type=JOIN.LEFT_OUTER)
                          .join_from(Visual_binary_gaia_dr3, RAVE_DR6_Gaia_DR3_XMatch,
                                     on=(source_col == RAVE_DR6_Gaia_DR3_XMatch.gaiae3),
                                     join_type=JOIN.LEFT_OUTER)
                          .join(RAVE_DR6_Gauguin_Madera,
                                on=(RAVE_DR6_Gaia_DR3_XMatch.obsid ==
                                    RAVE_DR6_Gauguin_Madera.rave_obs_id),
                                join_type=JOIN.LEFT_OUTER)
                          .where(CatalogToGaia_DR3.version_id == version_id,
                                 CatalogToGaia_DR3.best >> True,
                                 CatalogToTwoMassPSC.version_id == version_id,
                                 CatalogToTwoMassPSC.best >> True)
                          .distinct(Gaia_DR3.source_id))

            if 'apogee' in self.name:
                this_query = this_query.where(TwoMassPSC.h_m < 11)

            if query is None:
                query = this_query
            else:
                query |= this_query

        return query


class MWM_Wide_Binaries_APOGEE_Carton(MWM_Wide_Binaries_Base_Carton):
    """Wide binaries APOGEE carton."""

    name = 'mwm_bin_vis_apogee'
    instrument = 'APOGEE'
