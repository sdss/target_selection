#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2022-12-11
# @Filename: mwm_legacy_ir2opt.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    CatalogToTwoMassPSC,
    Gaia_DR3,
    SDSS_DR17_APOGEE_Allstarmerge,
    TwoMassPSC,
)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee model : postgres table name
# Catalog : catalogdb.catalog
# CatalogToTIC_v8 : catalogdb.catalog_to_tic_v8
# Gaia_DR2 : catalogdb.gaia_dr2_source
# SDSS_DR17_APOGEE_Allstarmerge : catalogdb.sdss_dr17_apogee_allstarmerge
# TIC_v8 : catalogdb.tic_v8
# TwoMassPSC : catalogdb.twomass_psc
#
#
# We associate each APOGEE target with a 2MASS source since the
# vast majority of APOGEE targets were selected from the 2MASS PSC.
# APOGEE_ID is the same as
# the “designation” column of the 2MASS PSC;
#
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# For sdss_dr17_apogee_allstarmerge we
# have to strip off the “2M” in the APOGEE_ID.
# However, some apogee_id values do not have 2M in the front.
# (old: for sdss_apogeeallstarmerge_r13 we also had to strip off the
#  “2M” in the APOGEE_ID.)
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# sdss5db=# select apogee_id from catalogdb.SDSS_DR17_APOGEE_Allstarmerge
# where apogee_id like '2M%';
#      apogee_id
# ---------------------
# 2M00000002+7417074
# 2M00000019-1924498
# etc.
# sdss5db=# select apogee_id from catalogdb.SDSS_DR17_APOGEE_Allstarmerge
# where apogee_id not like '2M%';
#    apogee_id
# --------------------
# 19140272-1554055
# 19155129-1617591
#
#
# Historical Note:
# The v0.5 version of this carton used catalogdb.sdss_apogeeallstarmerge_r13.
#


class MWM_Legacy_ir2opt_boss_Carton(BaseCarton):
    """
    Shorthand name: mwm_legacy_ir2opt_boss

    APOGEE targets to be observed with BOSS fibres.
    To be used as a filler for MWM-led plates.

    Simplified Description of selection criteria: Select all Gaia targets that
    have an APOGEE counterpart (i.e., have an entry in sdss_dr17_apogee_allstarmerge)
    with 13 < G < 18, BP > 13, and RP > 13.

    Wiki page: NA

    Additional source catalogs needed: Gaia DR2, sdss_dr17_apogee_allstarmerge

    Additional cross-matching needed: NA

    Return columns: catalog_id, source_id, apogee_id, phot_g_mean_mag,
    phot_rp_mean_mag, phot_bp_mean_mag
    """

    name = "mwm_legacy_ir2opt_boss"
    category = "science"
    instrument = "BOSS"
    cadence = "bright_1x1"
    program = "mwm_legacy"
    mapper = "MWM"
    priority = 6100
    can_offset = True

    # In the below query, we use replace() instead of ltrim() since
    # ltrim('2M20', '2M') will also trim the second 2.

    def build_query(self, version_id, query_region=None):
        query = (
            Catalog.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_DR3.phot_bp_mean_mag,
                Gaia_DR3.phot_rp_mean_mag,
                TwoMassPSC.pts_key,
                TwoMassPSC.designation,
            )
            .join(CatalogToGaia_DR3, on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .switch(CatalogToGaia_DR3)
            .join(
                CatalogToTwoMassPSC,
                on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid),
            )
            .join(TwoMassPSC, on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
            .join(
                SDSS_DR17_APOGEE_Allstarmerge,
                on=(
                    TwoMassPSC.designation
                    == peewee.fn.replace(SDSS_DR17_APOGEE_Allstarmerge.apogee_id, "2M", "")
                ),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToTwoMassPSC.version_id == version_id,
                CatalogToTwoMassPSC.best >> True,
                Gaia_DR3.phot_g_mean_mag.between(8, 18),
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
