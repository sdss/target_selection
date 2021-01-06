#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-12-14
# @Filename: ops_boss_stds.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog, CatalogToLegacy_Survey_DR8, CatalogToSDSS_DR13_PhotoObj_Primary,
    CatalogToTIC_v8, Gaia_DR2, Legacy_Survey_DR8, TIC_v8, TwoMassPSC,
    eBOSS_Target_v5)

from target_selection.cartons import BaseCarton


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

    name = 'ops_std_boss'
    category = 'standard'
    cadence = None
    program = 'ops_std'
    priority = 5425
    mapper = None

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

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid)
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.parallax > 0,
                        (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag) >= 0.65,
                        (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag) <= 0.8,
                        abs_gmag >= 3.5,
                        abs_gmag <= 5.5,
                        Gaia_DR2.phot_g_mean_mag > 13,
                        Gaia_DR2.phot_g_mean_mag < 18.5))
        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            # Here we do not need a join with Catalog since query already contains
            # a join with Catalog above.
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
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

    name = 'ops_std_boss_red'
    category = 'standard'
    cadence = None
    program = 'ops_std'
    priority = 5450
    mapper = None

    def build_query(self, version_id, query_region=None):

        # calculate proper motion amplitude and
        # term to construct reduced proper motion diagram
        pm_amp = peewee.fn.sqrt(Gaia_DR2.pmra * Gaia_DR2.pmra +
                                Gaia_DR2.pmdec * Gaia_DR2.pmdec)
        pm_Mod = 5.0 * peewee.fn.log(pm_amp) - 10.

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

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid, Catalog.ra, Catalog.dec,
                         Gaia_DR2.source_id, Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag, Gaia_DR2.phot_rp_mean_mag,
                         TwoMassPSC.j_m, TwoMassPSC.k_m,
                         ag.alias('ag'), ak.alias('ak'))
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
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
                        (Gaia_DR2.phot_rp_mean_mag + pm_Mod) >
                        (10.11 - 1.43 * (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag)),
                        (Gaia_DR2.phot_rp_mean_mag + pm_Mod) <
                        (13.11 - 1.43 * (Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag))))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
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

    name = 'ops_std_eboss'
    category = 'standard'
    cadence = None
    program = 'ops_std'
    priority = 5300
    mapper = None

    def build_query(self, version_id, query_region=None):

        selection_condition = (
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 50).
             cast("bigint")) > 0) |
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 51).
             cast("bigint")) > 0) |
            (eBOSS_Target_v5.eboss_target1.bin_and(peewee.fn.pow(2, 52).
             cast("bigint")) > 0))

        # We have distinct(eBOSS_Target_v5.objid_targeting) at the end
        # since the table ebosstarget_v5 has duplicate values.
        query = (Catalog
                 .select(Catalog.catalogid)
                 .join(CatalogToSDSS_DR13_PhotoObj_Primary,
                       on=(Catalog.catalogid ==
                           CatalogToSDSS_DR13_PhotoObj_Primary.catalogid))
                 .join(eBOSS_Target_v5,
                       on=(CatalogToSDSS_DR13_PhotoObj_Primary.target_id ==
                           eBOSS_Target_v5.objid_targeting))
                 .where(CatalogToSDSS_DR13_PhotoObj_Primary.version_id == version_id,
                        CatalogToSDSS_DR13_PhotoObj_Primary.best >> True,
                        selection_condition)
                 .distinct(eBOSS_Target_v5.objid_targeting))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
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

    name = 'ops_std_boss_tic'
    category = 'standard'
    cadence = None
    program = 'ops_std'
    priority = 5400
    mapper = None

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToTIC_v8.catalogid, Catalog.ra, Catalog.dec,
                         TwoMassPSC.h_m, TwoMassPSC.pts_key,
                         Gaia_DR2.phot_g_mean_mag,
                         TIC_v8.teff, TIC_v8.logg.alias('logg'),
                         TIC_v8.bmag, TIC_v8.vmag, TIC_v8.umag,
                         TIC_v8.gmag, TIC_v8.rmag, TIC_v8.imag,
                         TIC_v8.zmag,
                         peewee.fn.healpix_ang2ipix_nest(
                             128, Catalog.ra, Catalog.dec).alias('healpix_128'))
                 .join(CatalogToTIC_v8,
                       on=(Catalog.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(Gaia_DR2,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.phot_g_mean_mag > 13,
                        Gaia_DR2.phot_g_mean_mag < 17,
                        TIC_v8.teff > 6000,
                        TIC_v8.teff < 8000,
                        TIC_v8.logg > 3,
                        TIC_v8.logg < 5.5))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query

    def post_process(self, model):
        """
        The 10 highest gravity sources in each healpix are then
        selected and saved for the output carton.in each healpix pixel.
        """

        self.database.execute_sql("update sandbox.temp_ops_std_boss_tic " +
                                  "set selected = false")

        cursor = self.database.execute_sql(
            "select catalogid, healpix_128, logg from " +
            " sandbox.temp_ops_std_boss_tic " +
            " order by healpix_128 asc, logg desc;")

        output = cursor.fetchall()

        list_of_catalog_id = [0] * len(output)
        nside = 128
        total_number_healpix_pixels = 12 * nside * nside
        count = [0] * total_number_healpix_pixels
        current_target = 0
        for i in range(len(output)):
            current_healpix = output[i][1]
            if(count[current_healpix] < 10):
                count[current_healpix] = count[current_healpix] + 1
                list_of_catalog_id[current_target] = output[i][0]
                current_target = current_target + 1

        max_target = current_target
        for k in range(max_target + 1):
            self.database.execute_sql(
                " update sandbox.temp_ops_std_boss_tic set selected = true " +
                " where catalogid = " + str(list_of_catalog_id[k]) + ";")


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

    name = 'ops_std_boss_lsdr8'
    category = 'standard'
    cadence = None
    program = 'ops_std'
    priority = 5350
    mapper = None

    def build_query(self, version_id, query_region=None):

        # safety catch to avoid log-of-zero and divide by zero errors.
        # => use a flux in nano-maggies below which we give up
        nMgy_min = 1e-3  # equiv to AB=30

        # Below line is used to avoid divide by zero or log of zero,
        #     peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_g)
        # Below peewee.fn.log is log to the base 10.
        # peewee.fn.log(peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_r))

        # an alias to simplify accessing the query parameters:
        pars = self.parameters

        query = (Catalog
                 .select(Catalog.catalogid, Catalog.ra, Catalog.dec,
                         Legacy_Survey_DR8.ls_id,
                         Legacy_Survey_DR8.flux_g,
                         Legacy_Survey_DR8.flux_r,
                         Legacy_Survey_DR8.flux_z,
                         Legacy_Survey_DR8.flux_w1,
                         Legacy_Survey_DR8.flux_ivar_g,
                         Legacy_Survey_DR8.flux_ivar_r,
                         Legacy_Survey_DR8.flux_ivar_z,
                         Legacy_Survey_DR8.flux_ivar_w1,
                         Legacy_Survey_DR8.gaia_phot_g_mean_mag,
                         Legacy_Survey_DR8.gaia_phot_bp_mean_mag,
                         Legacy_Survey_DR8.gaia_phot_rp_mean_mag,
                         Legacy_Survey_DR8.parallax,
                         Legacy_Survey_DR8.parallax_ivar)
                 .join(CatalogToLegacy_Survey_DR8,
                       on=(Catalog.catalogid == CatalogToLegacy_Survey_DR8.catalogid))
                 .join(Legacy_Survey_DR8,
                       on=(CatalogToLegacy_Survey_DR8.target_id ==
                           Legacy_Survey_DR8.ls_id))
                 .where(CatalogToLegacy_Survey_DR8.version_id == version_id,
                        CatalogToLegacy_Survey_DR8.best >> True,
                        Legacy_Survey_DR8.type == 'PSF',
                        Legacy_Survey_DR8.ref_cat == 'G2',
                        Legacy_Survey_DR8.gaia_phot_g_mean_mag > pars['mag_gaia_g_min'],
                        Legacy_Survey_DR8.parallax
                        .between(pars['parallax_min'], pars['parallax_max']),
                        (Legacy_Survey_DR8.gaia_phot_bp_mean_mag -
                         Legacy_Survey_DR8.gaia_phot_rp_mean_mag)
                        .between(pars['gaia_bp_rp_min'], pars['gaia_bp_rp_max']),
                        Legacy_Survey_DR8.gaia_duplicated_source >> False,
                        Legacy_Survey_DR8.nobs_g >= 2,
                        Legacy_Survey_DR8.nobs_r >= 2,
                        Legacy_Survey_DR8.nobs_z >= 2,
                        Legacy_Survey_DR8.flux_g > nMgy_min,
                        Legacy_Survey_DR8.flux_r > nMgy_min,
                        Legacy_Survey_DR8.flux_z > nMgy_min,
                        Legacy_Survey_DR8.maskbits == 0,
                        (22.5 -
                         2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_r)))  # noqa: E501
                        .between(pars['mag_ls_r_min'], pars['mag_ls_r_max']),
                        (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_g) /  # noqa: E501
                                              peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_r)))  # noqa: E501
                        .between(pars['ls_g_r_min'], pars['ls_g_r_max']),
                        (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_r) /  # noqa: E501
                                              peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_z)))  # noqa: E501
                        .between(pars['ls_r_z_min'], pars['ls_r_z_max']),
                        (Legacy_Survey_DR8.gaia_phot_g_mean_mag -
                         (22.5 -
                          2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, Legacy_Survey_DR8.flux_r))))  # noqa: E501
                        .between(pars['gaia_g_ls_r_min'], pars['gaia_g_ls_r_max'])
                        ))

        # Below ra, dec and radius are in degrees
        # query_region[0] is ra of center of the region
        # query_region[1] is dec of center of the region
        # query_region[2] is radius of the region
        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query
