#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2024-07-01
# @Filename: mwm_hypervelocity.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import Catalog, CatalogToGaia_DR3, Gaia_DR3

from target_selection.cartons import BaseCarton


# from target_selection.exceptions import TargetSelectionError


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class Openfibertargets_mwm_hypervelocity_stars_boss_Carton(BaseCarton):
    """
    Shorthand name: openfibertargets_mwm_hypervelocity_stars_boss
    Simplified Description of selection criteria:
    select * from gaiadr3.gaia_source
    where phot_g_mean_mag < 19.5
    and ruwe < 1.4
    and pm > 50
    and ((4.74pm/parallax > 900) or (4.74pm/(parallax + parallax_error) > 600))
    Return columns: From Gaia DR3: phot_g_mean_mag, ra, dec, pm, ruwe,
    parallax, parallax_error, pmra, pmdec
    Metadata: can_offset=False
    Priority: 6085
    Cadence: bright_1x1 for sources with G < 18.5, and dark_1x1 for fainter sources
    Instrument: BOSS
    Lead contact: Kareem El-Badry
    Wiki page:
    https://sdss-wiki.atlassian.net/wiki/spaces/OPS/pages/13707660/Cartons+for+post+v1
    Additional source catalogs needed: Gaia_DR3

    """

    name = "openfibertargets_mwm_hypervelocity_stars_boss"
    category = "science"
    instrument = "BOSS"
    cadence = None  # cadence is assigned in post_process()
    program = "open_fiber"
    mapper = "MWM"
    priority = 6085
    can_offset = False

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag.alias("gaia_dr3_phot_g_mean_mag"),
                Gaia_DR3.pm,
                Gaia_DR3.pmra,
                Gaia_DR3.pmdec,
                Gaia_DR3.ruwe,
                Gaia_DR3.parallax,
                Gaia_DR3.parallax_error,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                Gaia_DR3.phot_g_mean_mag < 19.5,
                Gaia_DR3.ruwe < 1.4,
                Gaia_DR3.pm > 50,
                (
                    (4.74 * Gaia_DR3.pm / Gaia_DR3.parallax > 900)
                    | (4.74 * Gaia_DR3.pm / (Gaia_DR3.parallax + Gaia_DR3.parallax_error) > 600)
                ),  # noqa E501
            )
        )

        # There can be cases in which the same catalogid has multiple entries
        # in a catalog_to_x table since the same physical object
        # may match with multiple catalogids.
        # Hence, we have the below code in the above query.
        #                CatalogToGaia_DR3.best >> True,

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

    def post_process(self, model):
        """
        cadence options for these targets:
        bright_1x1 for sources with G < 18.5, and dark_1x1 for fainter sources
        """

        cursor = self.database.execute_sql(
            "select catalogid, gaia_dr3_phot_g_mean_mag from "
            + " sandbox.temp_openfibertargets_mwm_hypervelocity_stars_boss ;"
        )

        output = cursor.fetchall()

        for i in range(len(output)):
            current_catalogid = output[i][0]
            current_g = output[i][1]

            if current_g < 18.5:
                current_cadence = "bright_1x1"
            else:
                current_cadence = "dark_1x1"

            if current_cadence is not None:
                self.database.execute_sql(
                    " update sandbox.temp_openfibertargets_mwm_hypervelocity_stars_boss "
                    + " set cadence = '"
                    + current_cadence
                    + "'"
                    " where catalogid = " + str(current_catalogid) + ";"
                )
