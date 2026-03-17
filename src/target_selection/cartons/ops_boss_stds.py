#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-12-14
# @Filename: ops_boss_stds.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

# additional imports required by ops_boss_stds_ps1dr2
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToLegacy_Survey_DR8,
    CatalogToLegacy_Survey_DR10,
    CatalogToPanstarrs1,
    CatalogToSDSS_DR13_PhotoObj_Primary,
    CatalogToTIC_v8,
    Gaia_DR2,
    Legacy_Survey_DR8,
    Legacy_Survey_DR10,
    Panstarrs1,
    TIC_v8,
    TwoMassPSC,
    eBOSS_Target_v5,
)

from target_selection.cartons import BaseCarton
from target_selection.mag_flux import AB2Jy


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee Model name ---> postgres table name
# Gaia_DR2(CatalogdbModel)--->'gaia_dr2_source'
# TwoMassPSC(CatalogdbModel) --->'twomass_psc'
# eBOSS_Target_v5(CatalogdbModel)--->'ebosstarget_v5'
# Legacy_Survey_DR8(CatalogdbModel)--->'legacy_survey_dr8'
#
# The CatalogTo* peewee model names are not explicit in catalogdb.py.
# The names are constructed by prepending CatalogTo
# to the corresponding model name which is in catalogdb.py.
# For example:
# CatalogToSDSS_DR13_PhotoObj_Primary--->'catalog_to_sdss_dr13_photoobj_primary'
# CatalogToLegacy_Survey_DR8--->'catalogdb.catalog_to_legacy_survey_dr8'
#
# In the carton code, peewee.fn.log() is calling
# the PostgreSQL log() which is a base 10 logarithm.
# The below link has more details:
# https://www.postgresql.org/docs/12/functions-math.html


class OPS_BOSS_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds
    Selection criteria:
    --- criteria for ops_BOSS_stds ---
    #calculate distance modulus (could convert to use BailerJones distance)
    distMod = 5.*np.log10(1000./gaia_DR2.parallax)-5.

    #calculate absolute g and k magnitudes
    abs_gmag = gaia_DR2.phot_g_mean_mag - distMod
    abs_kmag = twomass_psc.k_m - distMod

    meet_std_criteria =
    np.where( ( ( (gaia_DR2.phot_bp_mean_mag -
                   gaia_DR2.phot_rp_mean_mag) >= 0.65) &
                  ( (gaia_DR2.phot_bp_mean_mag -
                     gaia_DR2.phot_rp_mean_mag) <= 0.8) &
                  ( (abs_gmag) >= 3.5) & ( (abs_gmag) <= 5.5) )
             )
    Also add the below condition:
             Gaia_DR2.phot_g_mean_mag > 13 AND
             Gaia_DR2.phot_g_mean_mag < 18.5

    Lead contact: Kevin Covey
    """

    name = "ops_std_boss"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5425
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        distMod = 5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0
        abs_gmag = Gaia_DR2.phot_g_mean_mag - distMod

        # We need a join with Catalog because
        # when using coordinates in the cartons we want to always
        # use Catalog.ra and Catalog.dec.
        # This is because those coordinates have
        # all been put in a common epoch 2015.5.
        # Below, we are using Catalog.ra and Catalog.dec
        # inside the "if query_region" block.

        query = (
            Catalog.select(CatalogToTIC_v8.catalogid)
            .join(CatalogToTIC_v8, on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                Gaia_DR2.parallax > 0,
                (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag) >= 0.65,
                (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag) <= 0.8,
                abs_gmag >= 3.5,
                abs_gmag <= 5.5,
                Gaia_DR2.phot_g_mean_mag > 13,
                Gaia_DR2.phot_g_mean_mag < 18.5,
            )
        )
        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            # Here we do not need a join with Catalog since query already contains
            # a join with Catalog above.
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


class OPS_BOSS_Red_Stds_Deredden_Carton(BaseCarton):
    """
        Shorthand name: ops_boss_red_stds_deredden

        Selection Criteria:
        This carton has selection criteria as defined in the below SQL code.

        SELECT c.catalogid, c.ra, c.dec, c.pmra, c.pmdec,
                twomass.j_m, twomass.h_m, twomass.k_m,
                gaia.phot_g_mean_mag, gaia.phot_bp_mean_mag, gaia.phot_rp_mean_mag,
                 gaia.parallax, gaia.pmra, gaia.pmdec

        FROM catalog c

        INNER JOIN catalog_to_tic_v8 ctic USING (catalogid)

        INNER JOIN tic_v8 tic ON tic.id = ctic.target_id

        INNER JOIN gaia_dr2_source gaia ON gaia.source_id = tic.gaia_int

        INNER JOIN twomass_psc twomass ON twomass.designation = tic.twomass_psc

        WHERE gaia.parallax > 0 AND gaia.phot_g_mean_mag < 18 AND

        /* enforce bp - rp edges of reddening strip in bp-rp vs. M_G space */

        gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag > 0.95 AND

        gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag < 3.625 AND

        /* enforce top and bottom of reddening strip in bp-rp vs. M_G space */

        gaia.phot_g_mean_mag - 5.*log(1000./gaia.parallax) + 5. >
         1.15+(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag)*1.89 AND

        gaia.phot_g_mean_mag - 5.*log(1000./gaia.parallax) + 5. <
         3.4+(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag)*1.89 AND

        /* enforce limits in empirical g-k color (make sure the source is at
         least a little reddened relative to a standard F star) */

        gaia.phot_g_mean_mag - twomass.k_m > 2.5 AND

        /* enforce limits in dereddened g-k color */

        (gaia.phot_g_mean_mag - 1.89*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) )
         - (twomass.k_m-0.186*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) ) > 1.2 AND

        (gaia.phot_g_mean_mag - 1.89*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) )
         - (twomass.k_m-0.186*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) ) < 1.8 AND

        /* enforce limits in dereddened j-k color */

        (twomass.j_m-0.582*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) )
         - (twomass.k_m-0.186*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) ) > 0.15 AND

        (twomass.j_m-0.582*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) )
         - (twomass.k_m-0.186*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) ) < 0.45 AND

        /* enforce limits in dereddened absolute k magnitude */

        twomass.k_m - 5.*log(1000./gaia.parallax) + 5. -
        0.186*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) > 2 AND

        twomass.k_m - 5.*log(1000./gaia.parallax) + 5. -
         0.186*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag - 0.7125) < 3.25 AND

        /* enforce proper motion limit */

        SQRT(gaia.pmra^2 + gaia.pmdec^2) > 3.5 AND

        (gaia.phot_rp_mean_mag + 5.*log(SQRT(gaia.pmra^2 + gaia.pmdec^2))-10.)
         > 10.11-1.43*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag) AND

        (gaia.phot_rp_mean_mag + 5.*log(SQRT(gaia.pmra^2 + gaia.pmdec^2))-10.)
         < 13.11-1.43*(gaia.phot_bp_mean_mag - gaia.phot_rp_mean_mag)
    ;
        Lead contact: Kevin Covey
    """

    name = "ops_std_boss_red"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5450
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        # calculate proper motion amplitude and
        # term to construct reduced proper motion diagram
        pm_amp = peewee.fn.sqrt(Gaia_DR2.pmra * Gaia_DR2.pmra + Gaia_DR2.pmdec * Gaia_DR2.pmdec)
        pm_Mod = 5.0 * peewee.fn.log(pm_amp) - 10.0

        # calculate a parallax-based distance modulus
        distMod = 5.0 * peewee.fn.log(1000.0 / Gaia_DR2.parallax) - 5.0

        # calculate the absolute g magnitude
        abs_gmag = Gaia_DR2.phot_g_mean_mag - distMod
        abs_kmag = TwoMassPSC.k_m - distMod

        # infer reddenings associated with assumed Fstar colors
        fstar_bp_rp = 0.725
        ag = 1.890 * (Gaia_DR2.bp_rp - fstar_bp_rp)
        aj = 0.582 * (Gaia_DR2.bp_rp - fstar_bp_rp)
        ak = 0.186 * (Gaia_DR2.bp_rp - fstar_bp_rp)

        query = (
            Catalog.select(
                CatalogToTIC_v8.catalogid,
                Catalog.ra,
                Catalog.dec,
                Gaia_DR2.source_id,
                Gaia_DR2.phot_g_mean_mag,
                Gaia_DR2.phot_bp_mean_mag,
                Gaia_DR2.phot_rp_mean_mag,
                TwoMassPSC.j_m,
                TwoMassPSC.k_m,
                ag.alias("ag"),
                ak.alias("ak"),
            )
            .join(CatalogToTIC_v8, on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .switch(TIC_v8)
            .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                # require a parallax so that absolute magnitudes are meaningful
                Gaia_DR2.parallax > 0,
                # require a slightly non-zero pm amplitude to
                # a) ensure reduced proper motion is meaningful, and
                # b) modestly filter quasars as candidates
                pm_amp > 3.5,
                # require Gaia G < 18 to provide reasonable SNR, and
                # G-K > 2 to ensure measured colors are at least
                # slightly reddened w.r.t regular f stars
                Gaia_DR2.phot_g_mean_mag < 18,
                Gaia_DR2.phot_g_mean_mag - TwoMassPSC.k_m > 2.5,
                # enforce bp - rp edges of reddening strip in
                # bp-rp vs. M_G space
                Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag > 0.95,
                Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag < 3.625,
                # enforce top and bottom of reddening strip
                # in bp-rp vs. M_G space
                abs_gmag > (1.15 + ag),
                abs_gmag < (3.4 + ag),
                # enforce limits in dereddened g-k color
                (Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak) > 1.2,
                (Gaia_DR2.phot_g_mean_mag - ag) - (TwoMassPSC.k_m - ak) < 1.8,
                # enforce limits in dereddened j-k color
                (TwoMassPSC.j_m - aj) - (TwoMassPSC.k_m - ak) > 0.15,
                (TwoMassPSC.j_m - aj) - (TwoMassPSC.k_m - ak) < 0.45,
                # enforce limits in dereddened absolute k magnitude
                abs_kmag - ak > 2,
                abs_kmag - ak < 3.25,
                # enforce halo-like motions in the reduced proper motion
                # diagram (selected in a strip consistent
                # w/ reddening the existing eBOSS standards)
                (Gaia_DR2.phot_rp_mean_mag + pm_Mod)
                > (10.11 - 1.43 * (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag)),
                (Gaia_DR2.phot_rp_mean_mag + pm_Mod)
                < (13.11 - 1.43 * (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag)),
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


class OPS_eBOSS_Stds_Carton(BaseCarton):
    """
    Shorthand name: ops_eboss_stds
    Selection Criteria:
    The code of this carton is based on the below SQL.
    SELECT DISTINCT e.objid_targeting
       FROM ebosstarget_v5 e JOIN catalog_to_sdss_dr13_photoobj_primary s
       ON s.target_id = e.objid_targeting
       WHERE (e.eboss_target1 & pow(2, 50)::bigint) > 0 OR
             (e.eboss_target1 & pow(2, 51)::bigint) > 0 OR
             (e.eboss_target1 & pow(2, 52)::bigint) > 0
             and s.best is true and s.version_id = 21;

    Lead contact: Kevin Covey
    """

    name = "ops_std_eboss"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5300
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        selection_condition = (
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 50).cast("bigint")) > 0)
            | (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 51).cast("bigint")) > 0)
            | (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 52).cast("bigint")) > 0)
        )

        # We have distinct(eBOSS_Target_v5.objid_targeting) at the end
        # since the table ebosstarget_v5 has duplicate values.
        query = (
            Catalog.select(Catalog.catalogid)
            .join(
                CatalogToSDSS_DR13_PhotoObj_Primary,
                on=(Catalog.catalogid == CatalogToSDSS_DR13_PhotoObj_Primary.catalogid),
            )
            .join(
                eBOSS_Target_v5,
                on=(
                    CatalogToSDSS_DR13_PhotoObj_Primary.target_id
                    == eBOSS_Target_v5.objid_targeting
                ),
            )
            .where(
                CatalogToSDSS_DR13_PhotoObj_Primary.version_id == version_id,
                CatalogToSDSS_DR13_PhotoObj_Primary.best >> True,
                selection_condition,
            )
            .distinct(eBOSS_Target_v5.objid_targeting)
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


class OPS_BOSS_Stds_TIC_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds_tic

    Simplified Description of selection criteria: parent catalog is TIC.

    13 < G < 17   (Gaia mag between 13 and 17)
    6000 < TIC teff < 8000     (temp between 6000-8000 K, so in
     late A - F - early G spectral range)
    3 < TIC log g < 5.5    (select to physically meaningful log g
     values to exclude nulls etc.)
    Once TIC has been limited to sources meeting the above criteria,
     the sky is then divided into an nside = 128 HEALPIX skymap.
    The 10 highest gravity sources in each healpix are then selected and
     saved for the output carton.

    All HEALPix use the nested ordering.

    Return columns: All the filter columns (Gaiamag, teff, logg),
    all other optical magnitudes in the
    TIC (bmag, vmag, umag, gmag, rmag, imag, zmag) plus healpix_128  and
    2MASS ID + H magnitude (if it exists, null if not;
    saved in case we need to use as fillers for high latitude BHM plates)

    Pseudo SQL (optional):

    Cadence options for these targets (list all options,
    even though no single target will receive more than one): N/A

    Lead contact:  Kevin Covey (but maybe also Marina Kounkel
    & Hector Javier Ibarra-Medel )
    """

    name = "ops_std_boss_tic"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5400
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        query = (
            Catalog.select(
                CatalogToTIC_v8.catalogid,
                Catalog.ra,
                Catalog.dec,
                TIC_v8.hmag,
                TIC_v8.twomass_psc,
                Gaia_DR2.phot_g_mean_mag,
                TIC_v8.teff,
                TIC_v8.logg.alias("logg"),
                TIC_v8.bmag,
                TIC_v8.vmag,
                TIC_v8.umag,
                TIC_v8.gmag,
                TIC_v8.rmag,
                TIC_v8.imag,
                TIC_v8.zmag,
                peewee.fn.healpix_ang2ipix_nest(128, Catalog.ra, Catalog.dec).alias("healpix_128"),
            )
            .join(CatalogToTIC_v8, on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
            .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
            .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
            .where(
                CatalogToTIC_v8.version_id == version_id,
                CatalogToTIC_v8.best >> True,
                Gaia_DR2.phot_g_mean_mag > 13,
                Gaia_DR2.phot_g_mean_mag < 17,
                TIC_v8.teff > 6000,
                TIC_v8.teff < 8000,
                TIC_v8.logg > 3,
                TIC_v8.logg < 5.5,
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
        The 10 highest gravity sources in each healpix are then
        selected and saved for the output carton.in each healpix pixel.
        """

        self.database.execute_sql("update sandbox.temp_ops_std_boss_tic " + "set selected = false")

        cursor = self.database.execute_sql(
            "select catalogid, healpix_128, logg from "
            + " sandbox.temp_ops_std_boss_tic "
            + " order by healpix_128 asc, logg desc;"
        )

        output = cursor.fetchall()

        list_of_catalog_id = [0] * len(output)
        nside = 128
        total_number_healpix_pixels = 12 * nside * nside
        count = [0] * total_number_healpix_pixels
        current_target = 0
        for i in range(len(output)):
            current_healpix = output[i][1]
            if count[current_healpix] < 10:
                count[current_healpix] = count[current_healpix] + 1
                list_of_catalog_id[current_target] = output[i][0]
                current_target = current_target + 1

        max_target = current_target
        for k in range(max_target + 1):
            self.database.execute_sql(
                " update sandbox.temp_ops_std_boss_tic set selected = true "
                + " where catalogid = "
                + str(list_of_catalog_id[k])
                + ";"
            )


class OPS_BOSS_Stds_LSDR8_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds_lsdr8
    – CRITERIA STILL IN VETTING STAGE -
    probably need to repeat all
    of the below in de-reddened colour space

    Comments: Spectrophotometric standards suitable for use by BOSS in dark time.

    Simplified Description of selection criteria:
    This simplified description is very high level.
    The carton is based on the pseudo SQL below.

    Parent catalog is legacy_survey_dr8 as ls, gaia_dr2_source as g2 with criteria:

    g2.G > 15.5
    15.95 < ls.r < 18.05
    -0.5 < g2.parallax < 1.0
    0.254 < ls.g - ls.r < 0.448
     0.024 < ls.r - ls.z < 0.190
    0.619 < g2.bp - g2.rp < 0.863
    0.0 < g2.G - ls.r < 0.10
    ls.type = 'PSF'
    ls.ref_cat = 'G2'
    ls.nobs_g > 2 && ls.nobs_r > = 2 && ls.nobs_z >= 2
    ls.maskbits = 0
    ls.gaia_duplicated_source = False
    Wiki page:
    All-sky BOSS standards#skyBOSSstandards-TransferringeBOSSselectiontolegacy_survey_dr8

    Return columns: TBD

    Pseudo SQL (optional):

    SELECT
        c.catalogid,l.ls_id,
        c.ra,c.dec,
        l.flux_g,l.flux_r,l.flux_z,l.flux_w1,
        l.flux_ivar_g,l.flux_ivar_r,l.flux_ivar_z,l.flux_ivar_w1,
        l.gaia_phot_g_mean_mag, l.gaia_phot_bp_mean_mag, l.gaia_phot_rp_mean_mag,
        l.parallax,l.parallax_ivar
    FROM catalog AS c
    JOIN catalog_to_legacy_survey_dr8 AS c2l
        ON c.catalogid = c2l.catalogid
    JOIN legacy_survey_dr8 AS l
        ON c2l.target_id = l.ls_id
    WHERE l.type = 'PSF'
        AND ref_cat = 'G2'
        AND gaia_phot_g_mean_mag > 15.5
        AND ((22.5 - 2.5*log10(greatest(1e-9,l.flux_r))) BETWEEN 15.95 AND 18.05)
        AND (l.parallax BETWEEN -0.5 AND 1.0)
        AND ((-2.5*log10(greatest(1e-9,l.flux_g)/
        greatest(1e-9,l.flux_r))) BETWEEN 0.254 AND 0.448)
        AND ((-2.5*log10(greatest(1e-9,l.flux_r)/
        greatest(1e-9,l.flux_z))) BETWEEN 0.024 AND 0.190)
        AND ((l.gaia_phot_bp_mean_mag-l.gaia_phot_rp_mean_mag)
         BETWEEN 0.619 AND 0.863)
        AND ((l.gaia_phot_g_mean_mag - (22.5-2.5*log10(greatest(1e-9,l.flux_r))))
         BETWEEN 0.0 AND 0.10)
        AND gaia_duplicated_source = false
        AND l.nobs_g >=2
        AND l.nobs_r >=2
        AND l.nobs_z >=2
        AND maskbits = 0
        AND c.version_id = ...
        AND c2l.best = true
    Cadence options for these targets (list all options,
     even though no single target will receive more than one): N/A

    Lead contact:  Tom Dwelly
    """

    name = "ops_std_boss_lsdr8"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5350
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()

        # an alias to simplify accessing the query parameters:
        pars = self.parameters

        # safety catch to avoid log-of-zero and divide by zero errors.
        # => use a flux in nano-maggies below which we give up
        nMgy_min = 1e-3  # equiv to AB=30

        # Below line is used to avoid divide by zero or log of zero,
        #     peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_g)
        # Below peewee.fn.log is log to the base 10.
        # peewee.fn.log(peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_r))

        # transform the legacysurvey grz into sdss psfmag griz
        # use transforms decribed here:
        # https://wiki.sdss.org/display/OPS/All-sky+BOSS+standards#All-skyBOSSstandards-TransformingphotometryofeBOSS-likestandardsintoSDSSsystem  # noqa
        # extract coeffs from fit logs via:
        # gawk 'BEGIN {print("coeffs = {")} /POLYFIT/{printf("\"%s%d\": %s,\n", substr($3,length($3)), $8, $10)} END {print("}")}'  ops_std_eboss/lsdr8_mag_to_sdss_psfmag_*.log  # noqa
        coeffs = {
            "g2": 0.193896,
            "g1": -0.051181,
            "g0": 0.032614,
            "i2": 0.044794,
            "i1": -0.513119,
            "i0": -0.021466,
            "r2": -0.034595,
            "r1": 0.132328,
            "r0": 0.011408,
            "z2": -0.381446,
            "z1": 0.135980,
            "z0": -0.020589,
        }

        g0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g))
        r0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r))
        z0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z))
        g_r = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g) / peewee.fn.greatest(nMgy_min, ls.flux_r)
        )
        r_z = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r) / peewee.fn.greatest(nMgy_min, ls.flux_z)
        )

        g = g0 + coeffs["g0"] + coeffs["g1"] * g_r + coeffs["g2"] * g_r * g_r
        r = r0 + coeffs["r0"] + coeffs["r1"] * g_r + coeffs["r2"] * g_r * g_r
        i = r0 + coeffs["i0"] + coeffs["i1"] * r_z + coeffs["i2"] * r_z * r_z
        z = z0 + coeffs["z0"] + coeffs["z1"] * r_z + coeffs["z2"] * r_z * r_z

        g0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g / ls.mw_transmission_g)
        )
        r0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
        )
        z0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_z / ls.mw_transmission_z)
        )

        g_r_dered = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g / ls.mw_transmission_g)
            / peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
        )
        r_z_dered = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
            / peewee.fn.greatest(nMgy_min, ls.flux_z / ls.mw_transmission_z)
        )

        # lsdr8 quotes ebv so we can go straight to E(G_BP-G_RP) and E(G_BP-G) directly
        # using the Stassun et al relation:
        # E(BP-RP) = 1.31 * E(B-V)
        # NO! # but we first need to recalibrate the SFD E(B-V) following Schlafly&Finkbeiner2011
        # NO! # such that SF11 E(B-V) = SFD EB(-V) * 0.884
        # NO! E_b_v_corr = 0.884
        # NO! E_bp_rp = 1.31 * E_b_v_corr
        E_bp_rp = 1.31
        R_gaia_g = 1.890 * E_bp_rp
        R_gaia_bp = 2.429 * E_bp_rp
        E_bp_g = R_gaia_bp - R_gaia_g  # = 0.7061

        bp_rp_dered = ls.gaia_phot_bp_mean_mag - ls.gaia_phot_rp_mean_mag - ls.ebv * E_bp_rp
        bp_g_dered = ls.gaia_phot_bp_mean_mag - ls.gaia_phot_g_mean_mag - ls.ebv * E_bp_g

        # bp_rp_dered = (
        #     (
        #         ls.gaia_phot_bp_mean_mag
        #         + 2.5 * peewee.fn.log(ls.mw_transmission_g)
        #     ) - (
        #         ls.gaia_phot_rp_mean_mag
        #         + 2.5 * peewee.fn.log(0.5 * (ls.mw_transmission_r + ls.mw_transmission_z))
        #     )
        # )
        # bp_g_dered = (
        #     (
        #         ls.gaia_phot_bp_mean_mag
        #         + 2.5 * peewee.fn.log(ls.mw_transmission_g)
        #     ) - (
        #         ls.gaia_phot_g_mean_mag
        #         + 2.5 * peewee.fn.log(ls.mw_transmission_r)
        #     )
        # )

        g_r_dered_nominal = pars["g_r_dered_nominal"]
        r_z_dered_nominal = pars["r_z_dered_nominal"]
        bp_rp_dered_nominal = pars["bp_rp_dered_nominal"]
        bp_g_dered_nominal = pars["bp_g_dered_nominal"]

        dered_dist2 = (
            (g_r_dered - g_r_dered_nominal) * (g_r_dered - g_r_dered_nominal)
            + (r_z_dered - r_z_dered_nominal) * (r_z_dered - r_z_dered_nominal)
            + (bp_rp_dered - bp_rp_dered_nominal) * (bp_rp_dered - bp_rp_dered_nominal)
            + (bp_g_dered - bp_g_dered_nominal) * (bp_g_dered - bp_g_dered_nominal)
        )
        dered_dist_max2 = pars["dered_dist_max"] * pars["dered_dist_max"]

        optical_prov = peewee.Value("sdss_psfmag_from_lsdr8")

        query = (
            Catalog.select(
                Catalog.catalogid,
                Catalog.ra,
                Catalog.dec,
                ls.ls_id,
                # ls.flux_g,
                # ls.flux_r,
                # ls.flux_z,
                # ls.flux_w1,
                # ls.flux_ivar_g,
                # ls.flux_ivar_r,
                # ls.flux_ivar_z,
                # ls.flux_ivar_w1,
                g0.alias("ls8_mag_g"),
                r0.alias("ls8_mag_r"),
                z0.alias("ls8_mag_z"),
                g_r.alias("ls8_mag_g_r"),
                r_z.alias("ls8_mag_r_z"),
                g0_dered.alias("ls8_mag_dered_g"),
                r0_dered.alias("ls8_mag_dered_r"),
                z0_dered.alias("ls8_mag_dered_z"),
                g_r_dered.alias("ls8_mag_dered_g_r"),
                r_z_dered.alias("ls8_mag_dered_r_z"),
                bp_rp_dered.alias("gdr2_mag_dered_bp_rp"),
                bp_g_dered.alias("gdr2_mag_dered_bp_g"),
                dered_dist2.alias("dered_dist2"),
                optical_prov.alias("optical_prov"),
                g.alias("g"),
                r.alias("r"),
                i.alias("i"),
                z.alias("z"),
                ls.gaia_phot_g_mean_mag.alias("gaia_g"),
                ls.gaia_phot_bp_mean_mag.alias("bp"),
                ls.gaia_phot_rp_mean_mag.alias("rp"),
                ls.ebv.alias("ls8_ebv"),
                ls.mw_transmission_g.alias("ls8_mw_transmission_g"),
                ls.mw_transmission_r.alias("ls8_mw_transmission_r"),
                ls.mw_transmission_z.alias("ls8_mw_transmission_z"),
                ls.parallax,
                ls.parallax_ivar,
                ls.nobs_g.alias("ls8_nobs_g"),
                ls.nobs_r.alias("ls8_nobs_r"),
                ls.nobs_z.alias("ls8_nobs_z"),
                # ls.maskbits,
            )
            .join(c2ls, on=(Catalog.catalogid == c2ls.catalogid))
            .join(ls, on=(c2ls.target_id == ls.ls_id))
            .where(
                c2ls.version_id == version_id,
                c2ls.best >> True,
                ls.type == "PSF",
                ls.ref_cat == "G2",
                ls.gaia_phot_g_mean_mag > pars["mag_gaia_g_min"],
                ls.parallax < pars["parallax_max"],
                ls.parallax
                > (
                    pars["parallax_min_at_g16"]
                    + (ls.gaia_phot_g_mean_mag - 16.0) * pars["parallax_min_slope"]
                ),
                ls.gaia_duplicated_source >> False,
                ls.nobs_g >= 1,  # TODO increase to >= 2 in lsdr9
                ls.nobs_r >= 1,  # TODO increase to >= 2 in lsdr9
                ls.nobs_z >= 1,  # TODO increase to >= 2 in lsdr9
                ls.flux_g > nMgy_min,
                ls.flux_r > nMgy_min,
                ls.flux_z > nMgy_min,
                ls.maskbits == 0,
                dered_dist2 < dered_dist_max2,
                r0.between(pars["mag_ls_r_min"], pars["mag_ls_r_max"]),
                #
                # g_r.between(pars['ls_g_r_min'], pars['ls_g_r_max']),
                # r_z.between(pars['ls_r_z_min'], pars['ls_r_z_max']),
                # (ls.gaia_phot_bp_mean_mag - ls.gaia_phot_rp_mean_mag)
                # .between(pars['gaia_bp_rp_min'], pars['gaia_bp_rp_max']),
                # (ls.gaia_phot_g_mean_mag - r0)
                # .between(pars['gaia_g_ls_r_min'], pars['gaia_g_ls_r_max'])
                # # (22.5 -
                # #  2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r)))  # noqa: E501
                # # .between(pars['mag_ls_r_min'], pars['mag_ls_r_max']),
                # # (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g) /  # noqa: E501
                # #                       peewee.fn.greatest(nMgy_min, ls.flux_r)))  # noqa: E501
                # # .between(pars['ls_g_r_min'], pars['ls_g_r_max']),
                # # (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r) /  # noqa: E501
                # #                       peewee.fn.greatest(nMgy_min, ls.flux_z)))  # noqa: E501
                # # .between(pars['ls_r_z_min'], pars['ls_r_z_max']),
                # # (ls.gaia_phot_g_mean_mag -
                # #  (22.5 -
                # #   2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r))))  # noqa: E501
                # # .between(pars['gaia_g_ls_r_min'], pars['gaia_g_ls_r_max'])
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
                ),
                peewee.fn.q3c_radial_query(
                    ls.ra, ls.dec, query_region[0], query_region[1], query_region[2]
                ),
            )
        return query


# ############################################################################


class OPS_BOSS_Stds_LSDR10_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds_lsdr10

    Comments: Spectrophotometric standards suitable for use by BOSS in dark time.

              Currently this is a copy of OPS_BOSS_Stds_LSDR8_Carton() with only name changes

    Lead contact:  Tom Dwelly
    """

    name = "ops_std_boss_lsdr10"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5350
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        ls = Legacy_Survey_DR10.alias()
        c2ls = CatalogToLegacy_Survey_DR10.alias()

        # an alias to simplify accessing the query parameters:
        pars = self.parameters

        # safety catch to avoid log-of-zero and divide by zero errors.
        # => use a flux in nano-maggies below which we give up
        nMgy_min = 1e-3  # equiv to AB=30

        # transform the legacysurvey grz into sdss psfmag griz
        # use transforms decribed here:
        # https://wiki.sdss.org/display/OPS/All-sky+BOSS+standards#All-skyBOSSstandards-TransformingphotometryofeBOSS-likestandardsintoSDSSsystem  # noqa
        # extract coeffs from fit logs via:
        # gawk 'BEGIN {print("coeffs = {")} /POLYFIT/{printf("\"%s%d\": %s,\n", substr($3,length($3)), $8, $10)} END {print("}")}'  ops_std_eboss/lsdr8_mag_to_sdss_psfmag_*.log  # noqa
        coeffs = {
            "g2": 0.193896,
            "g1": -0.051181,
            "g0": 0.032614,
            "i2": 0.044794,
            "i1": -0.513119,
            "i0": -0.021466,
            "r2": -0.034595,
            "r1": 0.132328,
            "r0": 0.011408,
            "z2": -0.381446,
            "z1": 0.135980,
            "z0": -0.020589,
        }

        g0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_g))
        r0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_r))
        i0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_i))
        z0 = 22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.flux_z))
        g_r = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g) / peewee.fn.greatest(nMgy_min, ls.flux_r)
        )
        r_i = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r) / peewee.fn.greatest(nMgy_min, ls.flux_i)
        )
        i_z = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_i) / peewee.fn.greatest(nMgy_min, ls.flux_z)
        )
        r_z = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r) / peewee.fn.greatest(nMgy_min, ls.flux_z)
        )

        g = g0 + coeffs["g0"] + coeffs["g1"] * g_r + coeffs["g2"] * g_r * g_r
        r = r0 + coeffs["r0"] + coeffs["r1"] * g_r + coeffs["r2"] * g_r * g_r
        i = r0 + coeffs["i0"] + coeffs["i1"] * r_z + coeffs["i2"] * r_z * r_z
        z = z0 + coeffs["z0"] + coeffs["z1"] * r_z + coeffs["z2"] * r_z * r_z

        g0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g / ls.mw_transmission_g)
        )
        r0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
        )
        i0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_i / ls.mw_transmission_i)
        )
        z0_dered = 22.5 - 2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_z / ls.mw_transmission_z)
        )

        g_r_dered = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_g / ls.mw_transmission_g)
            / peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
        )
        r_i_dered = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
            / peewee.fn.greatest(nMgy_min, ls.flux_i / ls.mw_transmission_i)
        )
        i_z_dered = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_i / ls.mw_transmission_i)
            / peewee.fn.greatest(nMgy_min, ls.flux_z / ls.mw_transmission_z)
        )
        r_z_dered = -2.5 * peewee.fn.log(
            peewee.fn.greatest(nMgy_min, ls.flux_r / ls.mw_transmission_r)
            / peewee.fn.greatest(nMgy_min, ls.flux_z / ls.mw_transmission_z)
        )

        # lsdr8 quotes ebv so we can go straight to E(G_BP-G_RP) and E(G_BP-G) directly
        # using the Stassun et al relation:
        # E(BP-RP) = 1.31 * E(B-V)
        # NO! # but we first need to recalibrate the SFD E(B-V) following Schlafly&Finkbeiner2011
        # NO! # such that SF11 E(B-V) = SFD EB(-V) * 0.884
        # NO! E_b_v_corr = 0.884
        # NO! E_bp_rp = 1.31 * E_b_v_corr
        E_bp_rp = 1.31
        R_gaia_g = 1.890 * E_bp_rp
        R_gaia_bp = 2.429 * E_bp_rp
        E_bp_g = R_gaia_bp - R_gaia_g  # = 0.7061

        bp_rp_dered = ls.gaia_phot_bp_mean_mag - ls.gaia_phot_rp_mean_mag - ls.ebv * E_bp_rp
        bp_g_dered = ls.gaia_phot_bp_mean_mag - ls.gaia_phot_g_mean_mag - ls.ebv * E_bp_g

        g_r_dered_nominal = pars["g_r_dered_nominal"]
        r_z_dered_nominal = pars["r_z_dered_nominal"]
        bp_rp_dered_nominal = pars["bp_rp_dered_nominal"]
        bp_g_dered_nominal = pars["bp_g_dered_nominal"]

        dered_dist2 = (
            (g_r_dered - g_r_dered_nominal) * (g_r_dered - g_r_dered_nominal)
            + (r_z_dered - r_z_dered_nominal) * (r_z_dered - r_z_dered_nominal)
            + (bp_rp_dered - bp_rp_dered_nominal) * (bp_rp_dered - bp_rp_dered_nominal)
            + (bp_g_dered - bp_g_dered_nominal) * (bp_g_dered - bp_g_dered_nominal)
        )
        dered_dist_max2 = pars["dered_dist_max"] * pars["dered_dist_max"]

        optical_prov = peewee.Value("sdss_psfmag_from_lsdr8")

        query = (
            Catalog.select(
                Catalog.catalogid,
                Catalog.ra,
                Catalog.dec,
                ls.ls_id,
                g0.alias("ls10_mag_g"),
                r0.alias("ls10_mag_r"),
                i0.alias("ls10_mag_i"),
                z0.alias("ls10_mag_z"),
                g_r.alias("ls10_mag_g_r"),
                r_i.alias("ls10_mag_r_i"),
                i_z.alias("ls10_mag_i_z"),
                r_z.alias("ls10_mag_r_z"),
                g0_dered.alias("ls10_mag_dered_g"),
                r0_dered.alias("ls10_mag_dered_r"),
                i0_dered.alias("ls10_mag_dered_i"),
                z0_dered.alias("ls10_mag_dered_z"),
                g_r_dered.alias("ls10_mag_dered_g_r"),
                r_i_dered.alias("ls10_mag_dered_r_i"),
                i_z_dered.alias("ls10_mag_dered_i_z"),
                r_z_dered.alias("ls10_mag_dered_r_z"),
                bp_rp_dered.alias("gaia_mag_dered_bp_rp"),
                bp_g_dered.alias("gaia_mag_dered_bp_g"),
                dered_dist2.alias("dered_dist2"),
                optical_prov.alias("optical_prov"),
                g.alias("g"),
                r.alias("r"),
                i.alias("i"),
                z.alias("z"),
                ls.gaia_phot_g_mean_mag.alias("gaia_g"),
                ls.gaia_phot_bp_mean_mag.alias("bp"),
                ls.gaia_phot_rp_mean_mag.alias("rp"),
                ls.ebv.alias("ls10_ebv"),
                ls.mw_transmission_g.alias("ls10_mw_transmission_g"),
                ls.mw_transmission_r.alias("ls10_mw_transmission_r"),
                ls.mw_transmission_i.alias("ls10_mw_transmission_i"),
                ls.mw_transmission_z.alias("ls10_mw_transmission_z"),
                ls.parallax,
                ls.parallax_ivar,
                ls.nobs_g.alias("ls10_nobs_g"),
                ls.nobs_r.alias("ls10_nobs_r"),
                ls.nobs_i.alias("ls10_nobs_i"),
                ls.nobs_z.alias("ls10_nobs_z"),
                ls.ref_cat.alias("ls10_ref_cat"),
                ls.ref_id.alias("ls10_ref_id"),
                # ls.maskbits,
            )
            .join(c2ls, on=(Catalog.catalogid == c2ls.catalogid))
            .join(ls, on=(c2ls.target_id == ls.ls_id))
            .where(
                c2ls.version_id == version_id,
                c2ls.best >> True,
                ls.type == "PSF",
                ((ls.ref_cat == "G2") | (ls.ref_cat == "GE")),
                ls.gaia_phot_g_mean_mag > pars["mag_gaia_g_min"],
                ls.parallax < pars["parallax_max"],
                ls.parallax
                > (
                    pars["parallax_min_at_g16"]
                    + (ls.gaia_phot_g_mean_mag - 16.0) * pars["parallax_min_slope"]
                ),
                ls.gaia_duplicated_source >> False,
                ls.nobs_g >= 1,  # TODO increase to >= 2 in lsdr9/10?
                ls.nobs_r >= 1,  # TODO increase to >= 2 in lsdr9/10?
                ls.nobs_z >= 1,  # TODO increase to >= 2 in lsdr9/10?
                ls.flux_g > nMgy_min,
                ls.flux_r > nMgy_min,
                ls.flux_z > nMgy_min,
                ls.maskbits == 0,
                dered_dist2 < dered_dist_max2,
                r0.between(pars["mag_ls_r_min"], pars["mag_ls_r_max"]),
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
                ),
                peewee.fn.q3c_radial_query(
                    ls.ra, ls.dec, query_region[0], query_region[1], query_region[2]
                ),
            )
        return query


# ############################################################################


# ############################################################################
class OPS_BOSS_Stds_PS1DR2_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds_ps1dr2

    Comments: Spectrophotometric standards suitable for use by BOSS in dark time,
          selected from panstarrs1-dr2 + gaia-dr2.

    Lead contact:  Tom Dwelly
    """

    name = "ops_std_boss_ps1dr2"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5351
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        ps = Panstarrs1.alias()
        c2ps = CatalogToPanstarrs1.alias()
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        # an alias to simplify accessing the query parameters:
        pars = self.parameters

        # transform the panstarrs1-dr2 griz into sdss psfmag griz
        # use transforms decribed here:
        # https://wiki.sdss.org/display/OPS/All-sky+BOSS+standards#All-skyBOSSstandards-TransformingphotometryofeBOSS-likestandardsintoSDSSsystem  # noqa
        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ pe=""; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  ops_std_eboss/ps1dr2_chp_psf_to_sdss_psfmag_?_results.log  # noqa

        coeffs = {
            "g2": 0.115563,
            "g1": 0.068765,
            "g0": 0.012047,
            "i2": -0.385214,
            "i1": 0.149677,
            "i0": -0.026127,
            "r2": -0.070151,
            "r1": 0.070129,
            "r0": -0.014197,
            "z2": -2.141255,
            "z1": 0.147746,
            "z0": -0.034845,
        }

        # start from ps1dr2 chp psf mags
        g_r = ps.g_chp_psf - ps.r_chp_psf
        r_i = ps.r_chp_psf - ps.i_chp_psf
        i_z = ps.i_chp_psf - ps.z_chp_psf

        # compute apparent sdss psfmags
        g = ps.g_chp_psf + coeffs["g0"] + coeffs["g1"] * g_r + coeffs["g2"] * g_r * g_r
        r = ps.r_chp_psf + coeffs["r0"] + coeffs["r1"] * g_r + coeffs["r2"] * g_r * g_r
        i = ps.i_chp_psf + coeffs["i0"] + coeffs["i1"] * r_i + coeffs["i2"] * r_i * r_i
        z = ps.z_chp_psf + coeffs["z0"] + coeffs["z1"] * i_z + coeffs["z2"] * i_z * i_z

        # dereddining steps

        # extinction terms for Panstarrs
        # Use R_b from Schlafly & Finkbeiner 2011, Table 6
        # https://ui.adsabs.harvard.edu/abs/2011ApJ...737..103S/abstract
        # assume R_V = 3.1
        # R_b = A_V / E(B-V)
        # for Panstarrs1 grizy bands we have:
        R_g = 3.172
        R_r = 2.271
        R_i = 1.682
        R_z = 1.322
        # R_y = 1.087
        E_g_r = R_g - R_r
        E_r_i = R_r - R_i
        E_i_z = R_i - R_z

        # extinction terms for Gaia
        # Use table 3 from Wang & Chen 2019
        # https://ui.adsabs.harvard.edu/abs/2019ApJ...877..116W/abstract
        # R_G = 1.890
        # R_BP = 2.429
        # R_RP = 1.429
        # These R_b are expressed with E(BP-RP) as the denominator - so need to convert to
        # using E(B-V) as denominator instead

        # an email from Keivan Stassun helps with this:
        #   Hi Jennifer, we had to work
        #   these out for the TIC paper because we used Gaia G mags and
        #   Bp-Rp colors for everything... you can see the details in
        #   Section 2.3.3 of Stassun et al (2019), but here's the final
        #   relations we adopted:
        #   E(${G}_{\mathrm{BP}}$ − ${G}_{\mathrm{RP}}$) = 1.31 E(B − V)
        #   AG = 2.72 E(B − V)

        E_bp_rp = 1.31
        R_gaia_g = 1.890 * E_bp_rp
        R_gaia_bp = 2.429 * E_bp_rp
        # R_gaia_rp = 1.429 * 1.31
        E_bp_g = R_gaia_bp - R_gaia_g  # = 0.7061

        # use ebv from tic_v8 match
        g_r_dered = g_r - tic.ebv * E_g_r
        r_i_dered = r_i - tic.ebv * E_r_i
        i_z_dered = i_z - tic.ebv * E_i_z
        bp_rp_dered = tic.gaiabp - tic.gaiarp - tic.ebv * E_bp_rp
        bp_g_dered = tic.gaiabp - tic.gaiamag - tic.ebv * E_bp_g

        g_r_dered_nominal = pars["g_r_dered_nominal"]
        r_i_dered_nominal = pars["r_i_dered_nominal"]
        i_z_dered_nominal = pars["i_z_dered_nominal"]
        bp_rp_dered_nominal = pars["bp_rp_dered_nominal"]
        bp_g_dered_nominal = pars["bp_g_dered_nominal"]

        dered_dist2 = (
            (g_r_dered - g_r_dered_nominal) * (g_r_dered - g_r_dered_nominal)
            + (r_i_dered - r_i_dered_nominal) * (r_i_dered - r_i_dered_nominal)
            + (i_z_dered - i_z_dered_nominal) * (i_z_dered - i_z_dered_nominal)
            + (bp_rp_dered - bp_rp_dered_nominal) * (bp_rp_dered - bp_rp_dered_nominal)
            + (bp_g_dered - bp_g_dered_nominal) * (bp_g_dered - bp_g_dered_nominal)
        )

        optical_prov = peewee.Value("sdss_psfmag_from_ps1dr2")

        ext_flags = 8388608 + 16777216
        dered_dist_max2 = pars["dered_dist_max"] * pars["dered_dist_max"]

        # the following are just to bracket
        # the result to make the query run faster
        r_stk_psf_flux_min = AB2Jy(pars["mag_ps_r_max"] + 0.2)
        r_stk_psf_flux_max = AB2Jy(pars["mag_ps_r_min"] - 0.2)

        query = (
            c2ps.select(
                c2ps.catalogid,
                tic.ra,
                tic.dec,
                ps.catid_objid.alias("ps1_catid_objid"),
                tic.gaia_int.alias("gaia_source"),
                ps.g_chp_psf.alias("ps1dr2_chp_psfmag_g"),
                ps.r_chp_psf.alias("ps1dr2_chp_psfmag_r"),
                ps.i_chp_psf.alias("ps1dr2_chp_psfmag_i"),
                ps.z_chp_psf.alias("ps1dr2_chp_psfmag_z"),
                g_r.alias("ps1dr2_chp_psfmag_g_r"),
                r_i.alias("ps1dr2_chp_psfmag_r_i"),
                i_z.alias("ps1dr2_chp_psfmag_i_z"),
                tic.ebv.alias("tic_ebv"),
                g_r_dered.alias("ps1dr2_chp_psfmag_g_r_dered"),
                r_i_dered.alias("ps1dr2_chp_psfmag_r_i_dered"),
                i_z_dered.alias("ps1dr2_chp_psfmag_i_z_dered"),
                bp_rp_dered.alias("gdr2_mag_dered_bp_rp"),
                bp_g_dered.alias("gdr2_mag_dered_bp_g"),
                dered_dist2.alias("dered_dist2"),
                optical_prov.alias("optical_prov"),
                g.alias("g"),
                r.alias("r"),
                i.alias("i"),
                z.alias("z"),
                tic.gaiamag.alias("gaia_g"),
                tic.gaiabp.alias("bp"),
                tic.gaiarp.alias("rp"),
                tic.jmag.alias("j"),
                tic.hmag.alias("h"),
                tic.kmag.alias("k"),
                tic.plx.alias("parallax"),
                tic.e_plx.alias("parallax_error"),
            )
            .join(ps, on=(c2ps.target_id == ps.catid_objid))
            .switch(c2ps)
            .join(c2tic, on=(c2ps.catalogid == c2tic.catalogid))
            .join(tic, on=(c2tic.target_id == tic.id))
            .where(
                c2ps.version_id == version_id,
                c2tic.version_id == version_id,
                c2ps.best >> True,
                c2tic.best >> True,
                ps.flags.bin_and(ext_flags) == 0,
                tic.plx < pars["parallax_max"],
                tic.plx
                > (
                    pars["parallax_min_at_g16"] + (tic.gaiamag - 16.0) * pars["parallax_min_slope"]
                ),
                dered_dist2 < dered_dist_max2,
                ps.r_chp_psf.between(pars["mag_ps_r_min"], pars["mag_ps_r_max"]),
                # the following are just attempts to bracket
                # the result to make the query run faster
                tic.gaiamag.between(pars["mag_gaia_g_min"], pars["mag_gaia_g_max"]),
                ps.r_stk_psf_flux.between(r_stk_psf_flux_min, r_stk_psf_flux_max),
                tic.ebv < pars["tic_ebv_max"],
                bp_rp_dered.between(pars["bp_rp_dered_min"], pars["bp_rp_dered_max"]),
            )
        )

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = query.join_from(CatalogToPanstarrs1, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


# ############################################################################
class OPS_BOSS_Stds_GaiaDR2_Carton(BaseCarton):
    """
    Shorthand name: ops_boss_stds_gdr2

    Comments: Spectrophotometric standards suitable for use by BOSS in dark time,
              selected from tic+gaia-dr2.

    Lead contact:  Tom Dwelly
    """

    name = "ops_std_boss_gdr2"
    category = "standard_boss"
    cadence = None
    program = "ops_std"
    priority = 5352
    mapper = None
    instrument = "BOSS"
    can_offset = False

    def build_query(self, version_id, query_region=None):
        tic = TIC_v8.alias()
        c2tic = CatalogToTIC_v8.alias()

        # an alias to simplify accessing the query parameters:
        pars = self.parameters

        # transform the gaia g,bp,rp into sdss psfmag griz
        # use transforms decribed here:
        # https://wiki.sdss.org/display/OPS/All-sky+BOSS+standards#All-skyBOSSstandards-TransformingphotometryofeBOSS-likestandardsintoSDSSsystem  # noqa
        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ pe=""; printf("\"%s%d%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  ops_std_eboss/gdr2_mag_to_sdss_psfmag_?_results.log  # noqa

        coeffs = {
            "g2": 0.226514,
            "g1": 0.373358,
            "g0": -0.073834,
            "i2": 0.038586,
            "i1": -0.505039,
            "i0": 0.216803,
            "r2": 0.212874,
            "r1": -0.381950,
            "r0": 0.156923,
            "z2": -0.246274,
            "z1": -0.372790,
            "z0": 0.235517,
        }

        bp_rp = tic.gaiabp - tic.gaiarp
        # compute apparent sdss psfmags
        g = tic.gaiamag + coeffs["g0"] + coeffs["g1"] * bp_rp + coeffs["g2"] * bp_rp * bp_rp
        r = tic.gaiamag + coeffs["r0"] + coeffs["r1"] * bp_rp + coeffs["r2"] * bp_rp * bp_rp
        i = tic.gaiamag + coeffs["i0"] + coeffs["i1"] * bp_rp + coeffs["i2"] * bp_rp * bp_rp
        z = tic.gaiamag + coeffs["z0"] + coeffs["z1"] * bp_rp + coeffs["z2"] * bp_rp * bp_rp

        # dereddining steps

        # extinction terms for Gaia
        # Use Stassun+19 (Ticv8) presciption
        E_bp_rp = 1.31
        R_gaia_g = 1.890 * E_bp_rp
        R_gaia_bp = 2.429 * E_bp_rp
        R_gaia_rp = 1.429 * E_bp_rp
        E_bp_g = R_gaia_bp - R_gaia_g  # = 0.7061
        E_g_rp = R_gaia_g - R_gaia_rp  # = 0.6039

        # use ebv from tic_v8 match
        bp_rp_dered = tic.gaiabp - tic.gaiarp - tic.ebv * E_bp_rp
        bp_g_dered = tic.gaiabp - tic.gaiamag - tic.ebv * E_bp_g
        g_rp_dered = tic.gaiamag - tic.gaiarp - tic.ebv * E_g_rp

        bp_rp_dered_nominal = pars["bp_rp_dered_nominal"]
        bp_g_dered_nominal = pars["bp_g_dered_nominal"]
        g_rp_dered_nominal = pars["g_rp_dered_nominal"]

        dered_dist2 = (
            (bp_rp_dered - bp_rp_dered_nominal) * (bp_rp_dered - bp_rp_dered_nominal)
            + (bp_g_dered - bp_g_dered_nominal) * (bp_g_dered - bp_g_dered_nominal)
            + (g_rp_dered - g_rp_dered_nominal) * (g_rp_dered - g_rp_dered_nominal)
        )

        optical_prov = peewee.Value("sdss_psfmag_from_gaia")

        dered_dist_max2 = pars["dered_dist_max"] * pars["dered_dist_max"]

        query = (
            Catalog.select(
                Catalog.catalogid,
                Catalog.ra,
                Catalog.dec,
                tic.id.alias("tic_id"),  # extra
                tic.gaia_int.alias("gaia_source"),  # extra
                tic.ebv.alias("tic_ebv"),  # extra
                bp_rp_dered.alias("gdr2_mag_dered_bp_rp"),  # extra
                bp_g_dered.alias("gdr2_mag_dered_bp_g"),  # extra
                g_rp_dered.alias("gdr2_mag_dered_g_rp"),  # extra
                dered_dist2.alias("dered_dist2"),  # extra
                optical_prov.alias("optical_prov"),
                g.alias("g"),
                r.alias("r"),
                i.alias("i"),
                z.alias("z"),
                tic.gaiamag.alias("gaia_g"),
                tic.gaiabp.alias("bp"),
                tic.gaiarp.alias("rp"),
                tic.jmag.alias("j"),
                tic.hmag.alias("h"),
                tic.kmag.alias("k"),
                tic.plx.alias("parallax"),  # extra
                tic.e_plx.alias("parallax_error"),  # extra
                tic.gallong.alias("tic_gal_l"),  # extra
                tic.gallat.alias("tic_gal_b"),  # extra
            )
            .join(c2tic, on=(Catalog.catalogid == c2tic.catalogid))
            .join(tic, on=(c2tic.target_id == tic.id))
            .where(
                c2tic.version_id == version_id,
                c2tic.best >> True,
                tic.plx < pars["parallax_max"],
                tic.plx
                > (
                    pars["parallax_min_at_g16"] + (tic.gaiamag - 16.0) * pars["parallax_min_slope"]
                ),
                dered_dist2 < dered_dist_max2,
                tic.gaiamag.between(pars["mag_gaia_g_min"], pars["mag_gaia_g_max"]),
                # the following are just to bracket the result to make the query run faster
                tic.gaiabp.between(pars["mag_gaia_bp_min"], pars["mag_gaia_bp_max"]),
                tic.gaiarp.between(pars["mag_gaia_rp_min"], pars["mag_gaia_rp_max"]),
                tic.ebv < pars["ebv_max"],
                ~(tic.gallat.between(-10.0, 10.0)),
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
                ),
            )

        return query
