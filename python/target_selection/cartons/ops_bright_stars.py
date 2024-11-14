#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2021-09-23
# @Filename: ops_bright_stars.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    CatalogToTycho2,
    Gaia_DR3,
    TwoMassPSC,
    Tycho2,
)

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
    "Select all objects from gaia DR3 with G < 13"

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

    name = "ops_gaia_brightneighbors"
    category = "veto_location_boss"
    instrument = None
    cadence = None
    program = "ops"
    mapper = None
    priority = None
    can_offset = False

    # target_selection propagates the following columns if they are generated
    # during the query (with exactly the below name)
    # g, r, i, z, h, j, k, bp, rp gaia_g, optical_prov
    #
    # If any g,r, i, z, optical_prov are missing, then the code will
    # try to find them in SDSS, PS1, and Gaia.
    #
    # Hence, we use the names gaia_g, bp, rp below.

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.pmra.alias("gaia_dr3_pmra"),
                Gaia_DR3.pmdec.alias("gaia_dr3_pmdec"),
                Gaia_DR3.phot_g_mean_mag.alias("gaia_g"),
                Gaia_DR3.phot_bp_mean_mag.alias("bp"),
                Gaia_DR3.phot_rp_mean_mag.alias("rp"),
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                Gaia_DR3.phot_g_mean_mag < 13,
            )
        )

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

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
    the source of these magnitudes.
    For the rows with the above transform, optical_prov = 'gaia_psfmag_tycho2a'.

    If VT is not null and BT is null then the Tycho2 VT magnitude will be
    transformed to pseudo gaia_g magnitudes calculated for
    the targetdb.magnitudes table using the ad-hoc transform:
                  G = VT - 1
     For the rows with the ad-hoc transform, optical_prov = 'gaia_psfmag_tycho2b'.
    """

    name = "ops_tycho2_brightneighbors"
    category = "veto_location_boss"
    instrument = None
    cadence = None
    program = "ops"
    mapper = None
    priority = None
    can_offset = False

    # The column tycho2.designation is the primary key of the table catalogdb.tycho2.
    # Hence below we use
    # on=(CatalogToTycho2.target_id == Tycho2.designation)
    #
    # We will not fill the g, r, i, z columns.
    # However, we will create these columns in the query here.
    # This is because if these columns exist then add_optical_magnitudes()
    # will return right away.
    # So add_optical_magnitudes() will not create the g, r, i, z, and
    # optical_prov columns and
    # it will not run queries to get g, r, i, z.
    #
    # We will fill the gaia_g and optical_prov columns in post_process().
    #
    # These columns must be part of the query as shown below.
    # They cannot be added later in post_process().

    def build_query(self, version_id, query_region=None):
        # In the below code,
        # in cast('text'), 'text' refers to the PostgreSQL column type 'text'
        # and
        # in cast('real'), 'real' refers to the PostgreSQL column type 'real'
        optical_prov = peewee.Value("sdss_psfmag_tycho2").cast("text")
        g = peewee.Value(None).cast("real")
        r = peewee.Value(None).cast("real")
        i = peewee.Value(None).cast("real")
        z = peewee.Value(None).cast("real")
        gaia_g = peewee.Value(None).cast("real")

        query = (
            CatalogToTycho2.select(
                CatalogToTycho2.catalogid,
                Tycho2.tycid,
                Tycho2.designation,
                Tycho2.ramdeg.alias("tycho2_ra"),
                Tycho2.demdeg.alias("tycho2_dec"),
                Tycho2.pmra.alias("tycho2_pmra"),
                Tycho2.pmde.alias("tycho2_pmde"),
                Tycho2.vtmag,
                Tycho2.btmag,
                optical_prov.alias("optical_prov"),
                g.alias("g"),
                r.alias("r"),
                i.alias("i"),
                z.alias("z"),
                gaia_g.alias("gaia_g"),
            )
            .join(Tycho2, on=(CatalogToTycho2.target_id == Tycho2.designation))
            .where(
                CatalogToTycho2.version_id == version_id,
                CatalogToTycho2.best >> True,
                Tycho2.vtmag < 13,
            )
        )

        if query_region:
            query = query.join_from(CatalogToTycho2, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query

    def post_process(self, model):
        """
        Compute new column gaia_g from tycho2 vtmag and btmag.
        """

        cursor = self.database.execute_sql(
            "select catalogid, vtmag, btmag from " + " sandbox.temp_ops_tycho2_brightneighbors ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            vtmag = output[i][1]
            btmag = output[i][2]

            if vtmag is not None:
                if btmag is not None:
                    current_gaia_g = (
                        vtmag
                        - 0.02051
                        - 0.2706 * (btmag - vtmag)
                        + 0.03394 * (btmag - vtmag) ** 2
                        - 0.05937 * (btmag - vtmag) ** 3
                    )
                    current_optical_prov = "gaia_psfmag_tycho2a"
                else:
                    # Since btmag is None, we cannot use the above equation.
                    # Below equation sets gaia_g to a bright value
                    current_gaia_g = vtmag - 1
                    current_optical_prov = "gaia_psfmag_tycho2b"
            else:
                raise TargetSelectionError(
                    "error: " + "ops_tycho2_brightneighbors post_process(): " + "vtmag is None"
                )

            self.database.execute_sql(
                " update sandbox.temp_ops_tycho2_brightneighbors "
                + " set gaia_g = "
                + str(current_gaia_g)
                + " where catalogid = "
                + str(current_catalogid)
                + ";"
            )

            self.database.execute_sql(
                " update sandbox.temp_ops_tycho2_brightneighbors "
                + " set optical_prov = '"
                + current_optical_prov
                + "'"
                + " where catalogid = "
                + str(current_catalogid)
                + ";"
            )


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

    name = "ops_2mass_psc_brightneighbors"
    category = "veto_location_apogee"
    instrument = None
    cadence = None
    program = "ops"
    mapper = None
    priority = None
    can_offset = False

    def build_query(self, version_id, query_region=None):
        # We do not select pmra and pmdec below because
        # twomass_psc table does not have pmra and pmdec.
        query = (
            CatalogToTwoMassPSC.select(
                CatalogToTwoMassPSC.catalogid,
                TwoMassPSC.designation.alias("twomass_psc_designation"),
                TwoMassPSC.ra.alias("twomass_psc_ra"),
                TwoMassPSC.decl.alias("twomass_psc_dec"),
                TwoMassPSC.j_m.alias("twomass_psc_j_m"),
                TwoMassPSC.h_m.alias("twomass_psc_h_m"),
                TwoMassPSC.k_m.alias("twomass_psc_k_m"),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .where(
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                TwoMassPSC.h_m < 7,
            )
        )

        if query_region:
            query = query.join_from(CatalogToTwoMassPSC, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query
