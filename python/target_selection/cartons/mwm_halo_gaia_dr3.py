#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-22
# @Filename: mwm_halo_gaia_dr3.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             Gaia_DR3,
                                             Gaia_dr3_vari_rrlyrae)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class MWM_mwm_halo_distant_rrl_Carton(BaseCarton):
    """ 5.1.30. mwm_halo_distant_rrl
Shorthand name: mwm_halo_distant_rrl
Existing carton code: replaces open fiber 2020 28c
Simplified Description of selection criteria:
all Gaia DR3 RRLs with BP < 18.8.
SELECT b.*, a.best_classification, a.pf, a.p1_o, a.zp_mag_g,
 a.average_rv, a.num_clean_epochs_g
FROM gaiadr3.vari_rrlyrae as a
INNER JOIN gaiadr3.gaia_source_lite as b ON a.source_id = b.source_id
WHERE b.phot_bp_mean_mag < 18.8
Return columns:
Metadata:
Priority: 3050 (may change after A/B tests)
Cadence: bright_1x1
Instrument: BOSS
can_offset = True
Lead contact: Alexander Ji
Assigned to: Pramod
Status: Ready for implementation
    """

    name = 'mwm_halo_distant_rrl'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'bright_1x1'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 3050
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_dr3_vari_rrlyrae.best_classification,
                         Gaia_dr3_vari_rrlyrae.pf,
                         Gaia_dr3_vari_rrlyrae.p1_o,
                         Gaia_dr3_vari_rrlyrae.zp_mag_g,
                         Gaia_dr3_vari_rrlyrae.average_rv,
                         Gaia_dr3_vari_rrlyrae.num_clean_epochs_g)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Gaia_dr3_vari_rrlyrae,
                       on=(Gaia_DR3.source_id == Gaia_dr3_vari_rrlyrae.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 18.8,))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_mwm_halo_distant_rrl_dark_Carton(BaseCarton):
    """  5.1.31. mwm_halo_distant_rrl_dark
Shorthand name: mwm_halo_distant_rrl_dark
Existing carton code:
Simplified Description of selection criteria:
all Gaia DR3 RRLs with BP < 18.8.
SELECT b.*, a.best_classification, a.pf, a.p1_o, a.zp_mag_g,
a.average_rv, a.num_clean_epochs_g
FROM gaiadr3.vari_rrlyrae as a
INNER JOIN gaiadr3.gaia_source_lite as b ON a.source_id = b.source_id
WHERE b.phot_bp_mean_mag < 18.8
Return columns:
Metadata:
Priority: 3049 (may change after A/B tests)
Cadence: dark_1x2
Instrument: BOSS
can_offset = True
Lead contact: Alexander Ji
    """

    name = 'mwm_halo_distant_rrl_dark'
    category = 'science'
    instrument = 'BOSS'
    cadence = 'dark_1x2'
    program = 'mwm_halo'
    mapper = 'MWM'
    priority = 3049
    can_offset = True

    def build_query(self, version_id, query_region=None):

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_dr3_vari_rrlyrae.best_classification,
                         Gaia_dr3_vari_rrlyrae.pf,
                         Gaia_dr3_vari_rrlyrae.p1_o,
                         Gaia_dr3_vari_rrlyrae.zp_mag_g,
                         Gaia_dr3_vari_rrlyrae.average_rv,
                         Gaia_dr3_vari_rrlyrae.num_clean_epochs_g)
                 .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .join(Gaia_dr3_vari_rrlyrae,
                       on=(Gaia_DR3.source_id == Gaia_dr3_vari_rrlyrae.source_id))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_bp_mean_mag < 18.8,))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
