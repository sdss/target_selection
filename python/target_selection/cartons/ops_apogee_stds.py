#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-11-14
# @Filename: ops_apogee_stds.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import Catalog, CatalogToTIC_v8, Gaia_DR2, TIC_v8, TwoMassPSC

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee Model name ---> postgres table name
# TwoMassPSC(CatalogdbModel) --->'twomass_psc'
#
# The CatalogTo* peewee model names are not explicit in catalogdb.py.
# The names are constructed by prepending CatalogTo
# to the corresponding model name which is in catalogdb.py.
# For example:
# CatalogToTIC_v8--->'catalog_to_tic_v8'


class OPS_APOGEE_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_apogee_stds
    Selection criteria:
    Shorthand name: ops_apogee_stds

    Simplified Description of selection criteria:
    parent cataloga ew 2MASS point source catalog and GAIA-DR2,
    after restricting to:

    7 < H < 11,
    -0.25 < J - K < 0.3,
    j_msigcom < 0.1, h_msigcom < 0.1, k_msigcom < 0.1,
    rd_flg = 2 or 1 in all three (JHK) bands
    ph_qual = A or B in all three (JHK) bands
    gal_contam = 0 (no galaxy contamination)
    prox > 6.0
    cc_flg = 000 (zero in all three (JHK) bands)
    ext_key = null (no match in extended source catalog)
    GAIA PM<100, GAIA Palarrax<10
    Once it has been limited to sources meeting the above criteria,
    the sky is then divided into an nside = 128 HEALPIX skymap.
    The 5 bluest sources (i.e. 5 smallest J-K) in each healpix are then
    selected and saved for the output carton.

    (note – the HEALPIX binning and sub-querying might be useful
    to implement as a function that can be run for future cartons;
    I can imagine that this capability might prove useful
    for updated BOSS calibrator cartons, and
    maybe for searching for guide stars, or
    even maybe some kind of science cartons that
    depend on sampling things smoothly in a spatial sense?)

    All HEALPix use the nested ordering.
    number_healpix_pixels = 12 *Nside*Nside
    Hence number_healpix_pixels = 12*128*128 = 196608

    Return columns: All the filter columns plus healpix_128
    Lead contact: Kevin Covey
    """

    name = "ops_std_apogee"
    category = "standard_apogee"
    cadence = None
    program = "ops_std"
    priority = 5500
    mapper = None
    instrument = "APOGEE"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        query = (
            Catalog.select(
                CatalogToTIC_v8.catalogid,
                Catalog.ra,
                Catalog.dec,
                TwoMassPSC.h_m,
                TwoMassPSC.j_m,
                TwoMassPSC.k_m,
                TwoMassPSC.h_msigcom,
                TwoMassPSC.j_msigcom,
                TwoMassPSC.k_msigcom,
                TwoMassPSC.rd_flg,
                TwoMassPSC.ph_qual,
                TwoMassPSC.gal_contam,
                TwoMassPSC.prox,
                TwoMassPSC.cc_flg,
                TwoMassPSC.ext_key,
                peewee.fn.healpix_ang2ipix_nest(128, Catalog.ra, Catalog.dec).alias("healpix_128"),
                (TwoMassPSC.j_m - TwoMassPSC.k_m).alias("j_minus_k"),
            )
            .join(CatalogToTIC_v8, on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                TwoMassPSC.h_m > 7,
                TwoMassPSC.h_m < 11,
                (TwoMassPSC.j_m - TwoMassPSC.k_m) > -0.25,
                TwoMassPSC.j_msigcom < 0.1,
                TwoMassPSC.h_msigcom < 0.1,
                TwoMassPSC.k_msigcom < 0.1,
                (TwoMassPSC.rd_flg == "111")
                | (TwoMassPSC.rd_flg == "211")
                | (TwoMassPSC.rd_flg == "121")
                | (TwoMassPSC.rd_flg == "221")
                | (TwoMassPSC.rd_flg == "112")
                | (TwoMassPSC.rd_flg == "212")
                | (TwoMassPSC.rd_flg == "122")
                | (TwoMassPSC.rd_flg == "222"),
                (TwoMassPSC.ph_qual == "AAA")
                | (TwoMassPSC.ph_qual == "BAA")
                | (TwoMassPSC.ph_qual == "ABA")
                | (TwoMassPSC.ph_qual == "BBA")
                | (TwoMassPSC.ph_qual == "AAB")
                | (TwoMassPSC.ph_qual == "BAB")
                | (TwoMassPSC.ph_qual == "ABB")
                | (TwoMassPSC.ph_qual == "BBB"),
                TwoMassPSC.gal_contam == 0,
                TwoMassPSC.prox > 6.0,
                TwoMassPSC.cc_flg == "000",
                TwoMassPSC.ext_key.is_null(),
                peewee.fn.sqrt(peewee.fn.pow(Gaia_DR2.pmra, 2) + peewee.fn.pow(Gaia_DR2.pmdec, 2))
                < 100,
                Gaia_DR2.parallax < 10,
            )
        )

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = query.where(
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
        Select the 5 bluest sources (i.e. 5 smallest J-K) in each healpix pixel.
        """

        self.database.execute_sql("update sandbox.temp_ops_std_apogee " + "set selected = false")

        cursor = self.database.execute_sql(
            "select catalogid, healpix_128, j_minus_k from "
            + " sandbox.temp_ops_std_apogee "
            + " order by healpix_128 asc, j_minus_k asc;"
        )

        output = cursor.fetchall()

        list_of_catalog_id = [0] * len(output)
        nside = 128
        total_number_healpix_pixels = 12 * nside * nside
        count = [0] * total_number_healpix_pixels
        current_target = 0
        for i in range(len(output)):
            current_healpix = output[i][1]
            if count[current_healpix] < 5:
                count[current_healpix] = count[current_healpix] + 1
                list_of_catalog_id[current_target] = output[i][0]
                current_target = current_target + 1

        max_target = current_target
        for k in range(max_target + 1):
            self.database.execute_sql(
                " update sandbox.temp_ops_std_apogee set selected = true "
                + " where catalogid = "
                + str(list_of_catalog_id[k])
                + ";"
            )
