#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2024-12-20
# @Filename: mwm_openfiber.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    Gaia_DR3,
)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class Openfibertargets_mwm_mdwarfs_plato_apogee_Carton(BaseCarton):
    """
    Shorthand name: openfibertargets_mwm_mdwarfs_plato_apogee
    Link: Stellar Parameters and Abundances for PLATO M-dwarf targets
    Simplified Description of selection criteria:
    ADQL query in Gaia archive for Gaia DR3:

    select sq.source_id, sq.ra, sq.dec
    from (
    select source_id, ra, dec,
    bp_rp as bp_rp_0,
    phot_g_mean_mag - (-0.17276 + 0.47885*(bp_rp) + -0.71953power(bp_rp,2) +
        0.24374power(bp_rp,3) + -0.04458power(bp_rp,4) +
        0.00317power(bp_rp,5)) as v0mag,
    phot_g_mean_mag + 5 * LOG10(parallax) - 10 as g0absmag
    from gaiadr3.gaia_source
    where
    l>=255.9-23 and l<=255.9+23 and b>=-24.6-23 and b<=-24.6+23
    and parallax > 0
    ) as sq
    where
    sq.v0mag <= 16
    and sq.g0absmag > -8.62sq.bp_rp_0 +24.96
    and sq.g0absmag >= 2.334sq.bp_rp_0 +2.259

    Return columns: Gaia DR3 galactic coordinates,
        Gaia BP-RP, parallax, G magnitude
    Metadata: can_offset=False
    Priority: 6085
    Cadence: bright_1x1
    Instrument: APOGEE
    Program: open_fiber
    Lead contact: Diogo Souto
    """

    name = "openfibertargets_mwm_mdwarfs_plato_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    program = "open_fiber"
    mapper = "MWM"
    priority = 6085
    can_offset = False

    # In the below query,
    # Gaia_DR3.l is for galactic longitude
    # Gaia_DR3.b is for galactic latitude
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
                Gaia_DR3.bp_rp,
                Gaia_DR3.parallax,
                Gaia_DR3.l,
                Gaia_DR3.b,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                Gaia_DR3.parallax > 0,
                Gaia_DR3.l >= 255.9 - 23,
                Gaia_DR3.l <= 255.9 + 23,
                Gaia_DR3.b >= -24.6 - 23,
                Gaia_DR3.b <= -24.6 + 23,
                Gaia_DR3.phot_g_mean_mag + 5 * peewee.fn.LOG10(Gaia_DR3.parallax) - 10
                > -8.62 * Gaia_DR3.bp_rp + 24.96,
                Gaia_DR3.phot_g_mean_mag + 5 * peewee.fn.LOG10(Gaia_DR3.parallax) - 10
                >= 2.334 * Gaia_DR3.bp_rp + 2.259,
                Gaia_DR3.phot_g_mean_mag
                - (
                    -0.17276
                    + 0.47885 * (Gaia_DR3.bp_rp)
                    + (-0.71953) * peewee.fn.power(Gaia_DR3.bp_rp, 2)
                    + 0.24374 * peewee.fn.power(Gaia_DR3.bp_rp, 3)
                    + (-0.04458) * peewee.fn.power(Gaia_DR3.bp_rp, 4)
                    + 0.00317 * peewee.fn.power(Gaia_DR3.bp_rp, 5)
                )
                < 16,
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
