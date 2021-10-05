#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2021-09-23
# @Filename: ops_bright_stars.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToTIC_v8,
                                             CatalogToTycho2, Gaia_DR2, TIC_v8,
                                             TwoMassPSC, TwoMassXSC, Tycho2)

from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class OPS_Gaia_Brightneighbors_Carton(BaseCarton):
    """6.1.1. Bright Gaia (G < 12) Stars
    Owner: Kevin Covey

    Shorthand name:
    ops_gaia_brightneighbors

    Simplified Description of selection criteria:
    "Select all objects from gaia DR2 with G < 13"

    Wiki page: NA

    Additional source catalogs needed: None

    Additional cross-matching needed:  None

    Return columns: Gaia id, Gaia RA, Gaia Dec, Gaia RA proper motion,
     Gaia Dec proper motion, G, BP, RP

    cadence options for these targets (list all options,
    even though no single target will receive more than one):
    Null (since this is a veto catalog,
    we want cadence, value, priority and instrument to all be Null).

    Notes:  Gaia magnitudes should be transformed onto the SDSS system
    for the targetdb.magnitude table in the standard way
    (including recording the appropriate value in the opt_prov field).

    """

    name = 'ops_gaia_brightneighbors'
    category = 'veto_location_boss'
    instrument = None
    cadence = None
    program = 'ops'
    mapper = None
    priority = None

    # target_selection propagates the following columns if they are generated
    # during the query (with exactly the below name)
    # g, r, i, z, h, j, k, bp, rp gaia_g, optical_prov
    #
    # If any g,r, i, z, optical_prov are missing, then the code will
    # try to find them in SDSS, PS1, and Gaia.
    #
    # Hence, we use the names gaia_g, bp, rp below.

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_dr2_ra'),
                         Gaia_DR2.dec.alias('gaia_dr2_dec'),
                         Gaia_DR2.pmra.alias('gaia_dr2_pmra'),
                         Gaia_DR2.pmdec.alias('gaia_dr2_pmdec'),
                         Gaia_DR2.phot_g_mean_mag.alias('gaia_g'),
                         Gaia_DR2.phot_bp_mean_mag.alias('bp'),
                         Gaia_DR2.phot_rp_mean_mag.alias('rp'))
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.phot_g_mean_mag < 13))

        # Gaia_DR2 peewee model class corresponds to
        # table catalogdb.gaia_dr2_source.

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class OPS_Tycho2_Brightneighbors_Carton(BaseCarton):
    """6.1.2. Bright Tycho2 (VT < 13) Stars
    Owner: Kevin Covey

    Shorthand name:
    ops_tycho2_brightneighbors

    Simplified Description of selection criteria:
    "Select all objects from Tycho2 with VT < 13"

    Wiki page: NA

    Additional source catalogs needed: None

    Additional cross-matching needed:  None

    Return columns: Tycho2 id, Tycho2 RA, Tycho2 Dec,
    Tycho2 RA proper motion, Tycho2 Dec proper motion, VT, BT

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    Null (since this is a veto catalog, we want cadence, value,
    priority and instrument to all be Null).

    Notes:  Tycho2 magnitudes VT, BT will be transformed to pseudo-gaia_g magnitudes
    calculated for the targetdb.magnitudes table
    using the transforms from Evans et al. (2018):
            G = VT - 0.02051 - 0.2706 * (BT - VT) +
            0.03394 * (BT - VT)^2 - 0.05937 * (BT - VT)^3
    all other magnitudes can be stored as 'null',
    and a new optical_prov entry should be used to indicate
    the source of these magnitudes (e.g., 'gaia_psfmag_tycho2')
    """

    name = 'ops_tycho2_brightneighbors'
    category = 'veto_location_boss'
    instrument = None
    cadence = None
    program = 'ops'
    mapper = None
    priority = None

    # The column tycho2.designation is the primary key of the table catalogdb.tycho2.
    # Hence below we use
    # on=(CatalogToTycho2.target_id == Tycho2.designation)
    #
    # optical_prov must be part of the query as shown below.
    # It cannot be added later as a column in post_process().

    def build_query(self, version_id, query_region=None):

        optical_prov = peewee.Value('gaia_psfmag_tycho2')
        query = (CatalogToTycho2
                 .select(CatalogToTycho2.catalogid,
                         Tycho2.tycid,
                         Tycho2.designation,
                         Tycho2.ramdeg.alias('tycho2_ra'),
                         Tycho2.demdeg.alias('tycho2_dec'),
                         Tycho2.pmra.alias('tycho2_pmra'),
                         Tycho2.pmde.alias('tycho2_pmde'),
                         Tycho2.vtmag,
                         Tycho2.btmag,
                         optical_prov.alias('optical_prov'))
                 .join(Tycho2, on=(CatalogToTycho2.target_id == Tycho2.designation))
                 .where(CatalogToTycho2.version_id == version_id,
                        CatalogToTycho2.best >> True,
                        Tycho2.vtmag < 13))

        if query_region:
            query = (query
                     .join_from(CatalogToTycho2, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):
        """
        Compute new column gaia_g from tycho2 vtmag and btmag.
        """

        self.database.execute_sql(
            "alter table sandbox.temp_ops_tycho2_brightneighbors " +
            " add column gaia_g double precision ;")

        # self.database.execute_sql(
        #    "alter table sandbox.temp_ops_tycho2_brightneighbors " +
        #    " add column optical_prov text ;")

        cursor = self.database.execute_sql(
            "select catalogid, vtmag, btmag from " +
            " sandbox.temp_ops_tycho2_brightneighbors ;")

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            vtmag = output[i][1]
            btmag = output[i][2]

            if (vtmag is not None):
                if(btmag is not None):
                    current_gaia_g = (vtmag - 0.02051 -
                                      0.2706 * (btmag - vtmag) +
                                      0.03394 * (btmag - vtmag)**2 -
                                      0.05937 * (btmag - vtmag)**3)
                    current_optical_prov = 'gaia_psfmag_tycho2a'
                else:
                    # Since btmag is None, we cannot use the above equation.
                    # Below equation sets gaia_g to a bright value
                    current_gaia_g = vtmag - 1
                    current_optical_prov = 'gaia_psfmag_tycho2b'
            else:
                raise TargetSelectionError(
                      'error: ' +
                      'ops_tycho2_brightneighbors post_process(): ' +
                      'vtmag is None')

            self.database.execute_sql(
                " update sandbox.temp_ops_tycho2_brightneighbors " +
                " set gaia_g = " + str(current_gaia_g) +
                " where catalogid = " + str(current_catalogid) + ";")

            self.database.execute_sql(
                " update sandbox.temp_ops_tycho2_brightneighbors " +
                " set optical_prov = '" + current_optical_prov + "'" +
                " where catalogid = " + str(current_catalogid) + ";")


class OPS_2MASS_PSC_Brightneighbors_Carton(BaseCarton):
    """6.2.  Bright 2MASS (H < 7) Point Sources
    Owner: Kevin Covey

    Shorthand name:
    ops_2mass_psc_brightneighbors

    Simplified Description of selection criteria:
    "Select all objects from the 2MASS Point Source Catalog with H < 7"

    Wiki page: NA

    Additional source catalogs needed: None

    Additional cross-matching needed:  None

    Return columns:  2MASS ID, 2MASS RA, 2MASS Dec, J, H, K

    cadence options for these targets
    (list all options, even though no single target will receive more than one):
    Null (since this is a veto catalog, we want cadence, value,
    priority and instrument to all be Null).
    """

    name = 'ops_2mass_psc_brightneighbors'
    category = 'veto_location_apogee'
    instrument = None
    cadence = None
    program = 'ops'
    mapper = None
    priority = None

    def build_query(self, version_id, query_region=None):

        # We do not select pmra and pmdec below because
        # twomass_psc table does not have pmra and pmdec.
        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         TwoMassPSC.designation.alias('twomass_psc_designation'),
                         TwoMassPSC.ra.alias('twomass_psc_ra'),
                         TwoMassPSC.decl.alias('twomass_psc_dec'),
                         TwoMassPSC.j_m.alias('twomass_psc_j_m'),
                         TwoMassPSC.h_m.alias('twomass_psc_h_m'),
                         TwoMassPSC.k_m.alias('twomass_psc_k_m'))
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 7))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class OPS_2MASS_XSC_Brightneighbors_Carton(BaseCarton):
    """6.3.  Bright 2MASS (H < 7) Extended Sources
Owner: Kevin Covey

Shorthand name:
ops_2mass_xsc_brightneighbors

Simplified Description of selection criteria:
"Select all objects from the 2MASS Extended Source Catalog with h_m_k20fe < 7"

Wiki page: NA

Additional source catalogs needed: None

Additional cross-matching needed:  None

Return columns: 2MASS ID, 2MASS RA, 2MASS Dec,
j_m_k20fe, h_m_k20fe, k_m_k20fe

cadence options for these targets
(list all options, even though no single target will receive more than one):
Null (since this is a veto catalog, we want cadence, value,
priority and instrument to all be Null).

    """

    name = 'ops_2mass_xsc_brightneighbors'
    category = 'veto_location_apogee'
    instrument = None
    cadence = None
    program = 'ops'
    mapper = None
    priority = None

    def build_query(self, version_id, query_region=None):

        # We do not select pmra and pmdec below because
        # twomass_xsc table does not have pmra and pmdec.
        # For table twomass_xsc, below use TIC_v8.twomass
        # instead of TIC_v8.twomass_psc.
        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid,
                         TwoMassXSC.designation.alias('twomass_xsc_designation'),
                         TwoMassXSC.ra.alias('twomass_xsc_ra'),
                         TwoMassXSC.decl.alias('twomass_xsc_dec'),
                         TwoMassXSC.j_m_k20fe.alias('twomass_xsc_j_m_k20fe'),
                         TwoMassXSC.h_m_k20fe.alias('twomass_xsc_h_m_k20fe'),
                         TwoMassXSC.k_m_k20fe.alias('twomass_xsc_k_m_k20fe'))
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassXSC, on=(TIC_v8.twomass == TwoMassXSC.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassXSC.h_m_k20fe < 7))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
