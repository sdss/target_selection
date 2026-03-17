#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-21
# @Filename: mwm_bin_gaia.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToGaia_DR3,
    Gaia_DR3,
    Gaia_dr3_nss_two_body_orbit,
)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class MWM_bin_gaia_astb_apogee_Carton(BaseCarton):
    """5.1.21. mwm_bin_gaia_astb_apogee
    Shorthand name:  mwm_bin_gaia_astb_apogee
    Existing carton code: adql_query_mwm_bin_gaia_astb_apogee.txt (this is
    ADQL code that runs on the Gaia archive) Simplified Description of
    selection criteria:  binaries from gaiadr3.nss_two_body_orbit satisfying
    (period < 1000) & (8 < phot_g_mean_mag < 13) & ((nss_solution_type =
    'Orbital') or (nss.nss_solution_type = 'AstroSpectroSB1')).
    This gives 62,738 sources.

    Return columns:
     query for full carton (returns 62,738)
     select gs.source_id, gs.ra, gs.dec, gs.phot_g_mean_mag
     from gaiadr3.nss_two_body_orbit as nss,
      gaiadr3.gaia_source as gs
      where (nss.nss_solution_type = 'Orbital' or
      nss.nss_solution_type = 'AstroSpectroSB1') and
      gs.source_id = nss.source_id and gs.phot_g_mean_mag between 8 and 13
      and period < 1000
    Metadata:
    Priority: 2550-2559
    Cadence: 2 visits separated by at least 100 days
     – new cadence    cadence assessed to be incompatible with existing set.
      Plan to test implementation with simplified cadence
      (bright_2x1_long_v2  - and allow partial completion / extra epoch eligibility)
    Instrument: APOGEE
    can_offset = True
    Lead contact: Kareem El-Badry
    """

    name = "mwm_bin_gaia_astb_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_2x1_long_v2"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2550
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_dr3_nss_two_body_orbit.period,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(
                Gaia_dr3_nss_two_body_orbit,
                on=(Gaia_DR3.source_id == Gaia_dr3_nss_two_body_orbit.source_id),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "Orbital")
                | (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "AstroSpectroSB1"),
                Gaia_DR3.phot_g_mean_mag.between(8, 13),
                Gaia_dr3_nss_two_body_orbit.period < 1000,
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


class MWM_bin_gaia_astb_boss_Carton(BaseCarton):
    """5.1.22. mwm_bin_gaia_astb_boss
    Shorthand name:  mwm_bin_gaia_astb_boss
    Existing carton code: adql_query_mwm_bin_gaia_astb_boss.txt
    Simplified Description of selection criteria:  binaries from
    gaiadr3.nss_two_body_orbit satisfying (period < 1000) & (13 <
    phot_g_mean_mag < 16) & ((nss_solution_type = 'Orbital') or
    (nss.nss_solution_type = 'AstroSpectroSB1')).
    This gives 87,470 sources.

    Return columns:
     query for full carton (returns 87470)
     select gs.source_id, gs.ra,
     gs.dec, gs.phot_g_mean_mag from gaiadr3.nss_two_body_orbit as nss,
     gaiadr3.gaia_source as gs where (nss.nss_solution_type = 'Orbital' or
     nss.nss_solution_type = 'AstroSpectroSB1')and gs.source_id =
     nss.source_id and gs.phot_g_mean_mag between 13 and 16 and period < 1000
    Metadata:
    Priority: 2560-2569
    Cadence:
    2 visits separated by at least 100 days – new cadence    cadence
    assessed to be incompatible with existing set. Plan to test
    implementation with simplified cadence (bright_2x1 / dark_2x1 - and
    allow partial completion / extra epoch eligibility)
    Instrument: BOSS
    can_offset = True
    Lead contact: Kareem El-Badry
    """

    name = "mwm_bin_gaia_astb_boss"
    category = "science"
    instrument = "BOSS"
    cadence = "bright_2x1_long_v2"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2560
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_dr3_nss_two_body_orbit.period,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(
                Gaia_dr3_nss_two_body_orbit,
                on=(Gaia_DR3.source_id == Gaia_dr3_nss_two_body_orbit.source_id),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "Orbital")
                | (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "AstroSpectroSB1"),
                Gaia_DR3.phot_g_mean_mag.between(13, 16),
                Gaia_dr3_nss_two_body_orbit.period < 1000,
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


class MWM_bin_gaia_sb_apogee_Carton(BaseCarton):
    """5.1.23. mwm_bin_gaia_sb_apogee
    Shorthand name:  mwm_bin_gaia_sb_apogee Existing carton code:
    adql_query_mwm_bin_gaia_sb_apogee.txt Simplified Description of
    selection criteria:  binaries from gaiadr3.nss_two_body_orbit satisfying
    (period < 1000) & (8 < phot_g_mean_mag < 11.5) & ((nss_solution_type =
    'SB1') or (nss_solution_type = 'SB1C')).
    This gives 51,810 sources.

    Return columns:
     query for full carton (returns 51810)
     select gs.source_id, gs.ra,
     gs.dec, gs.phot_g_mean_mag from gaiadr3.nss_two_body_orbit as nss,
     gaiadr3.gaia_source as gswhere ((nss.nss_solution_type = 'SB1') or
     (nss.nss_solution_type = 'SB1C'))and gs.source_id = nss.source_idand
     gs.phot_g_mean_mag between 8 and 11.5 and period < 1000
    Metadata:
    Priority: 2570-2579
    Cadence: bright_1x1
    Instrument: APOGEE
    can_offset = True
    Lead contact: Kareem El-Badry
    """

    name = "mwm_bin_gaia_sb_apogee"
    category = "science"
    instrument = "APOGEE"
    cadence = "bright_1x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2570
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_dr3_nss_two_body_orbit.period,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(
                Gaia_dr3_nss_two_body_orbit,
                on=(Gaia_DR3.source_id == Gaia_dr3_nss_two_body_orbit.source_id),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "SB1")
                | (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "SB1C"),
                Gaia_DR3.phot_g_mean_mag.between(8, 11.5),
                Gaia_dr3_nss_two_body_orbit.period < 1000,
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


class MWM_bin_gaia_sb_boss_Carton(BaseCarton):
    """5.1.24. mwm_bin_gaia_sb_boss
    Shorthand name:  mwm_bin_gaia_sb_boss Existing carton code:
    adql_query_mwm_bin_gaia_sb_boss.txt Simplified Description of selection
    criteria:  binaries from gaiadr3.nss_two_body_orbit satisfying (period <
    1000) & (11.5 < phot_g_mean_mag < 14) & ((nss_solution_type = 'SB1') or
    (nss_solution_type = 'SB1C')).
    This gives 111,852 sources.

    Return columns:
     query for full carton (returns 111,852)
     select gs.source_id, gs.ra,
     gs.dec, gs.phot_g_mean_mag from gaiadr3.nss_two_body_orbit as nss,
     gaiadr3.gaia_source as gswhere ((nss.nss_solution_type = 'SB1') or
     (nss.nss_solution_type = 'SB1C'))and gs.source_id = nss.source_idand
     gs.phot_g_mean_mag between 11.5 and 14 and period < 1000
    Metadata:
    Priority: 2580-2589
    Cadence: bright_1x1
    Instrument: BOSS
    can_offset = True
    Lead contact: Kareem El-Badry
    """

    name = "mwm_bin_gaia_sb_boss"
    category = "science"
    instrument = "BOSS"
    cadence = "bright_1x1"
    program = "mwm_bin"
    mapper = "MWM"
    priority = 2580
    can_offset = True

    def build_query(self, version_id, query_region=None):
        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                Gaia_DR3.source_id,
                Gaia_DR3.ra.alias("gaia_dr3_ra"),
                Gaia_DR3.dec.alias("gaia_dr3_dec"),
                Gaia_DR3.phot_g_mean_mag,
                Gaia_dr3_nss_two_body_orbit.period,
            )
            .join(Gaia_DR3, on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
            .join(
                Gaia_dr3_nss_two_body_orbit,
                on=(Gaia_DR3.source_id == Gaia_dr3_nss_two_body_orbit.source_id),
            )
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "SB1")
                | (Gaia_dr3_nss_two_body_orbit.nss_solution_type == "SB1C"),
                Gaia_DR3.phot_g_mean_mag.between(11.5, 14),
                Gaia_dr3_nss_two_body_orbit.period < 1000,
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
