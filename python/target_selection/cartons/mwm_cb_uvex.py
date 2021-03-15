#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-20
# @Filename: mwm_cb_uvex.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (UVOT_SSC_1, XMM_OM_SUSS_4_1,
                                             Catalog, CatalogToGUVCat,
                                             CatalogToTIC_v8,
                                             CatalogToUVOT_SSC_1,
                                             CatalogToXMM_OM_SUSS_4_1,
                                             Gaia_DR2,
                                             GeometricDistances_Gaia_DR2,
                                             GUVCat, TIC_v8, TwoMassPSC)

from target_selection.cartons import BaseCarton


fuv_mag = GUVCat.fuv_mag
nuv_mag = GUVCat.nuv_mag
fuv_magerr = GUVCat.fuv_magerr
nuv_magerr = GUVCat.nuv_magerr
fuv_nuv = fuv_mag - nuv_mag

BJ = GeometricDistances_Gaia_DR2

absg = (Gaia_DR2.phot_g_mean_mag - 5 * fn.log(BJ.r_est) + 5)
bp_rp = Gaia_DR2.bp_rp

gmagab = Gaia_DR2.phot_g_mean_mag - 25.688 + 25.792
fuvg = fuv_mag - gmagab

B_AB = Gaia_DR2.phot_bp_mean_mag - 25.351 + 25.386
R_AB = Gaia_DR2.phot_rp_mean_mag - 24.762 + 25.116
B_AB_R_AB = B_AB - R_AB

colour_absg = ((absg > 4.09) | (absg > 4.5457 * bp_rp + 4.0457))

pmra = Gaia_DR2.pmra
pmdec = Gaia_DR2.pmdec
pmra_error = Gaia_DR2.pmra_error
pmdec_error = Gaia_DR2.pmdec_error

pm = fn.sqrt(pmra * pmra + pmdec * pmdec)
pm_error = fn.sqrt(fn.pow((pmra_error * pmra / pm), 2) +
                   fn.pow((pmdec_error * pmdec / pm), 2))
logpmdpm = fn.log(pm / pm_error)

epar = Gaia_DR2.parallax / Gaia_DR2.parallax_error

astrometric_cuts = (BJ.r_lo <= 1500,
                    ~(((logpmdpm < 0.301) &
                       (epar > -1.4995 * logpmdpm - 4.05) &
                       (epar < 1.4995 * logpmdpm + 4.05)) |
                      (fn.pow((logpmdpm - 0.301), 2) / (0.39794 * 0.39794) +
                       fn.pow(epar, 2) / (4.5 * 4.5) <= 1)))


class MWM_CB_UVEX_BaseCarton(BaseCarton):

    def post_process(self, model, **kwargs):

        self.database.execute_sql(f'ALTER TABLE {self.path} ADD COLUMN value REAL;')
        model._meta.add_field('value', peewee.FloatField())

        # G < 17 => cadence = bright_1x1
        model.update(cadence='bright_1x1').where(model.phot_g_mean_mag < 17).execute()
        model.update(value=1).where(model.phot_g_mean_mag < 17).execute()

        # 17 < G < 19 => cadence = dark_1x2
        model.update(cadence='dark_1x2').where((model.phot_g_mean_mag > 17) &
                                               (model.phot_g_mean_mag < 19)).execute()
        model.update(value=2).where((model.phot_g_mean_mag > 17) &
                                    (model.phot_g_mean_mag < 19)).execute()

        # G > 19 => cadence = dark_1x3
        model.update(cadence='dark_1x3').where(model.phot_g_mean_mag > 19).execute()
        model.update(value=3).where(model.phot_g_mean_mag > 19).execute()

        return model


class MWM_CB_UVEX1_Carton(MWM_CB_UVEX_BaseCarton):
    """MWM Compact Binaries UV excess 1.

    Definition:

       Match Gaia - GALEX, use both bands of GALEX (FUV and NUV),
       apply color cuts to select objects with UV excess between main
       sequence and WD sequence. Remove objects with low parallax and
       proper motion accuracy to not get swamped by AGN and QSO.

    Pseudo-SQL:

        - Match Gaia DR2 with GALEX within 5 arcsec, keep the nearest match
          only, add distance and lower limit distance (columns r_est, r_lo;
          estimated and lower limit to estimated distance) from catalog
          gdr2_contrib/geometric_distance; use columns fuv_mag, nuv_mag,
          fuv_magerr, nuv_magerr from GALEX catalog.

        - Definitions:
            - pm = sqrt(pmra*pmra+pmdec*pmdec)
            - pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm) +
                         (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            - logpmdpm = log10(pm/pm_error)
            - gmagab = phot_g_mean_mag-25.688+25.792
            - bmagab = phot_bp_mean_mag-25.351+25.386
            - rmagab = phot_rp_mean_mag-24.762+25.116
            - B_AB-R_AB = bmagab-rmagab
            - absg=phot_g_mean_mag-5*log10(r_est)+5
            - nuvg=nuv_mag-gmagab
            - fuvg=fuv_mag-gmagab
            - fuvnub=fuv_mag-nuv_mag
            - bp_rp=phot_bp_mean_mag-phot_rp_mean_mag

        - General cut:
            - Gaia: visibility_periods_used >5

        - Colour and magnitude cuts:
            - nuv_magerr < 0.2 && nuv_mag>-100 &&
                (fuv_magerr >=0.2 || fuv_mag< -100)
            - absg < (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
            - ((nuv_mag-gmagab) < 2.25 ||
                ((nuv_mag-gmagab) < 6.725*(B_AB-R_AB)-1.735  &&
                 (nuv_mag-gmagab)< -0.983*(B_AB-R_AB)+8.24))

        - Astrometric cuts:
            - !((logpmdpm <0.301 &&
                (parallax/parallax_error)>-1.4995*logpmdpm-4.05 &&
                (parallax/parallax_error)<1.4995*logpmdpm+4.05) ||
               (logpmdpm-0.301)*(logpmdpm-0.301)/(0.39794*0.39794) +
               (parallax/parallax_error)^2/(4.5*4.5)<=1)
            - r_lo <= 1500

    """

    name = 'mwm_cb_uvex1'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def setup_transaction(self):

        self.database.execute_sql('SET LOCAL join_collapse_limit = 1;')
        super().setup_transaction()

    def build_query(self, version_id, query_region=None):

        colour_cuts = (fuv_magerr < 0.2,
                       nuv_magerr < 0.2,
                       fuv_mag > -100,
                       nuv_mag > -100,
                       colour_absg,
                       (absg > -1.11749253e-03 * fn.pow(fuv_nuv, 3) +
                        1.53748615e-02 * fn.pow(fuv_nuv, 2) +
                        3.66419895e-01 * fuv_nuv + 2.20026639),
                       ((fuvg < 6.08) | ((fuvg < 11.82 * B_AB_R_AB + 2.58) &
                                         (fuvg < -0.79 * B_AB_R_AB + 9.21))))

        query = (GUVCat
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_ra'),
                         Gaia_DR2.dec.alias('gaia_dec'),
                         Gaia_DR2.pmra,
                         Gaia_DR2.pmdec,
                         Gaia_DR2.pmra_error,
                         Gaia_DR2.pmdec_error,
                         Gaia_DR2.parallax,
                         Gaia_DR2.parallax_error,
                         Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         BJ.r_est,
                         BJ.r_lo,
                         GUVCat.objid.alias('guvcat_objid'),
                         GUVCat.ra.alias('guvcat_ra'),
                         GUVCat.dec.alias('guvcat_dec'),
                         GUVCat.nuv_mag,
                         GUVCat.fuv_mag,
                         GUVCat.nuv_magerr,
                         GUVCat.fuv_magerr)
                 .join(CatalogToGUVCat)
                 .join(CatalogToTIC_v8, on=(CatalogToGUVCat.catalogid ==
                                            CatalogToTIC_v8.catalogid))
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(BJ)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(CatalogToGUVCat.version_id == version_id,
                        CatalogToGUVCat.best >> True)
                 .where(Gaia_DR2.visibility_periods_used > 5)
                 .where(*colour_cuts)
                 .where(*astrometric_cuts))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_UVEX2_Carton(MWM_CB_UVEX_BaseCarton):
    """MWM Compact Binaries UV excess 2.

    Definition:

        Match Gaia - GALEX, use NUV channel of GALEX only, apply color cuts
        to select objects with UV excess between main sequence and WD sequence.
        Remove objects with low parallax and proper motion accuracy to not
        get swamped by AGN and QSO.

    Pseudo-SQL:

        - Match Gaia DR2 with GALEX within 5 arcsec, keep the nearest match
          only, add distance and lower limit distance (columns r_est, r_lo;
          estimated and lower limit to estimated distance) from catalog
          gdr2_contrib/geometric_distance; use columns fuv_mag, nuv_mag,
          fuv_magerr, nuv_magerr from GALEX catalog.

        - Definitions:
            - pm = sqrt(pmra*pmra+pmdec*pmdec)
            - pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm) +
                         (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            - logpmdpm = log10(pm/pm_error)
            - gmagab = phot_g_mean_mag-25.688+25.792
            - bmagab = phot_bp_mean_mag-25.351+25.386
            - rmagab = phot_rp_mean_mag-24.762+25.116
            - B_AB-R_AB = bmagab-rmagab
            - absg=phot_g_mean_mag-5*log10(r_est)+5
            - nuvg=nuv_mag-gmagab
            - fuvg=fuv_mag-gmagab
            - fuvnub=fuv_mag-nuv_mag
            - bp_rp=phot_bp_mean_mag-phot_rp_mean_mag

        - General cut:
            - Gaia: visibility_periods_used >5

        - Colour and magnitude cuts:
            - fuv_magerr < 0.2 && nuv_magerr < 0.2 &&
              fuv_mag > -100 && nuv_mag > -100
            - absg < (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
            - absg > -1.11749253e-03*fuv_nuv^3 +
                1.53748615e-02*fuv_nuv^2 +
                3.66419895e-01*fuv_nuv+2.20026639e+00
            - (fuvg < 6.08 || (fuvg < 11.82*(B_AB-R_AB) +
                               2.58 && fuvg < -0.79*(B_AB-R_AB)+9.21))

        - Astrometric cuts:
            - !((logpmdpm <0.301 &&
                (parallax/parallax_error)>-1.4995*logpmdpm-4.05 &&
                (parallax/parallax_error)<1.4995*logpmdpm+4.05) ||
               (logpmdpm-0.301)*(logpmdpm-0.301)/(0.39794*0.39794) +
               (parallax/parallax_error)^2/(4.5*4.5)<=1)
            - r_lo <= 1500

    """

    name = 'mwm_cb_uvex2'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def setup_transaction(self):

        self.database.execute_sql('SET LOCAL join_collapse_limit = 1;')
        super().setup_transaction()

    def build_query(self, version_id, query_region=None):

        colour_cuts = (nuv_magerr < 0.2,
                       nuv_mag > -100,
                       colour_absg,
                       (((nuv_mag - gmagab) < 2.25) |
                        (((nuv_mag - gmagab) < 6.725 * B_AB_R_AB - 1.735) &
                         ((nuv_mag - gmagab) < -0.983 * B_AB_R_AB + 8.24))))

        query = (GUVCat
                 .select(CatalogToGUVCat.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_ra'),
                         Gaia_DR2.dec.alias('gaia_dec'),
                         Gaia_DR2.pmra,
                         Gaia_DR2.pmdec,
                         Gaia_DR2.pmra_error,
                         Gaia_DR2.pmdec_error,
                         Gaia_DR2.parallax,
                         Gaia_DR2.parallax_error,
                         Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         BJ.r_est,
                         BJ.r_lo,
                         GUVCat.objid.alias('guvcat_objid'),
                         GUVCat.ra.alias('guvcat_ra'),
                         GUVCat.dec.alias('guvcat_dec'),
                         GUVCat.nuv_mag,
                         GUVCat.nuv_magerr)
                 .join(CatalogToGUVCat)
                 .join(CatalogToTIC_v8, on=(CatalogToGUVCat.catalogid ==
                                            CatalogToTIC_v8.catalogid))
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(BJ)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(CatalogToGUVCat.version_id == version_id,
                        CatalogToGUVCat.best >> True)
                 .where(Gaia_DR2.visibility_periods_used > 5)
                 .where(*colour_cuts)
                 .where(*astrometric_cuts))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_UVEX3_Carton(MWM_CB_UVEX_BaseCarton):
    """MWM Compact Binaries UV excess 3.

    Definition:

       Match Gaia - XMMOMSUOB (or OM SUSS v4.1), apply color cuts to select
       objects with UV excess between main sequence and WD sequence.
       Remove objects with low parallax and proper motion accuracy
       to not get swamped by AGN and QSO.

    Pseudo-SQL:

        - Apply some general UV quality flags and work on subset of entries
          in OM-SUSS:
            - (UVW2_AB_MAG >-100) || (UVW1_AB_MAG >-100) || (UVM2_AB_MAG>-100)
            - !(equals(substring(UVM2_QUALITY_FLAG_ST,0,1),"T") ||
                (equals(substring(UVM2_QUALITY_FLAG_ST,1,2),"T") &&
                 UVM2_SIGNIF <10) ||
                (equals(substring(UVM2_QUALITY_FLAG_ST,2,3),"T") &&
                 UVM2_SIGNIF <10) ||
                equals(substring(UVM2_QUALITY_FLAG_ST,6,7),"T") ||
                equals(substring(UVM2_QUALITY_FLAG_ST,7,8),"T") ||
                equals(substring(UVM2_QUALITY_FLAG_ST,8,9),"T"))

        - Match UV catalog with Gaia DR2 within 3 arcsec, keep the nearest
          match only, add distance and lower limit distance (columns r_est,
          r_lo; estimated and lower limit to estimated distance)
          from catalog gdr2_contrib/geometric_distance.

        - Definitions:
            - pm = sqrt(pmra*pmra+pmdec*pmdec)
            - pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm) +
                         (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            - logpmdpm = log10(pm/pm_error)
            - gmagab = phot_g_mean_mag-25.688+25.792
            - bmagab = phot_bp_mean_mag-25.351+25.386
            - rmagab = phot_rp_mean_mag-24.762+25.116
            - B_AB-R_AB = bmagab-rmagab
            - absg=phot_g_mean_mag-5*log10(r_est)+5
            - nuvg=nuv_mag-gmagab
            - fuvg=fuv_mag-gmagab
            - fuvnub=fuv_mag-nuv_mag
            - bp_rp=phot_bp_mean_mag-phot_rp_mean_mag

        - General cut:
            - Gaia: visibility_periods_used >5

        - Colour and magnitude cuts:
            - absg < (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
            - ((UVM2_AB_MAG - gmagab) < 2.25 ||
               ((UVM2_AB_MAG - gmagab) < 6 &&
                (UVM2_AB_MAG - gmagab) < 5.57377*(B_AB-R_AB)+0.2049))

        - Astrometric cuts:
            - !((logpmdpm <0.301 &&
                (parallax/parallax_error)>-1.4995*logpmdpm-4.05 &&
                (parallax/parallax_error)<1.4995*logpmdpm+4.05) ||
               (logpmdpm-0.301)*(logpmdpm-0.301)/(0.39794*0.39794) +
               (parallax/parallax_error)^2/(4.5*4.5)<=1)
            - r_lo <= 1500

    """

    name = 'mwm_cb_uvex3'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def build_query(self, version_id, query_region=None):

        uvm2_ab_mag = XMM_OM_SUSS_4_1.uvm2_ab_mag
        uvm2_quality_flag_st = XMM_OM_SUSS_4_1.uvm2_quality_flag_st
        uvm2_signif = XMM_OM_SUSS_4_1.uvm2_signif

        quality_flags = (
            (uvm2_ab_mag > -100),
            ~((fn.substr(uvm2_quality_flag_st, 1, 1) == 'T') |
              ((fn.substr(uvm2_quality_flag_st, 2, 1) == 'T') & (uvm2_signif < 10)) |
              ((fn.substr(uvm2_quality_flag_st, 3, 1) == 'T') & (uvm2_signif < 10)) |
              (fn.substr(uvm2_quality_flag_st, 7, 1) == 'T') |
              (fn.substr(uvm2_quality_flag_st, 8, 1) == 'T') |
              (fn.substr(uvm2_quality_flag_st, 9, 1) == 'T'))
        )

        quality_cte = (XMM_OM_SUSS_4_1
                       .select(XMM_OM_SUSS_4_1,
                               CatalogToXMM_OM_SUSS_4_1.catalogid)
                       .join(CatalogToXMM_OM_SUSS_4_1)
                       .where(CatalogToXMM_OM_SUSS_4_1.version_id == version_id,
                              CatalogToXMM_OM_SUSS_4_1.best >> True)
                       .where(*quality_flags)
                       .cte('quality_cte', materialized=True))

        colour_cuts = (
            colour_absg,
            ((quality_cte.c.uvm2_ab_mag - gmagab) < 2.25) |
            (((quality_cte.c.uvm2_ab_mag - gmagab) < 6) &
             ((quality_cte.c.uvm2_ab_mag - gmagab) < 5.57377 * (B_AB_R_AB) + 0.2049)))

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_ra'),
                         Gaia_DR2.dec.alias('gaia_dec'),
                         Gaia_DR2.pmra,
                         Gaia_DR2.pmdec,
                         Gaia_DR2.pmra_error,
                         Gaia_DR2.pmdec_error,
                         Gaia_DR2.parallax,
                         Gaia_DR2.parallax_error,
                         Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         gmagab.alias('gmagab'),
                         BJ.r_est,
                         BJ.r_lo,
                         quality_cte.c.srcnum,
                         quality_cte.c.obsid,
                         quality_cte.c.ra.alias('xmm_ra'),
                         quality_cte.c.dec.alias('xmm_dec'),
                         quality_cte.c.uvm2_signif,
                         quality_cte.c.uvm2_ab_mag.alias('uvm2_ab_mag'),
                         quality_cte.c.uvm2_quality_flag_st.alias('uvm2_quality_flag_st'))
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(BJ)
                 .join_from(Catalog, quality_cte,
                            on=(quality_cte.c.catalogid == Catalog.catalogid))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(Gaia_DR2.visibility_periods_used > 5)
                 .where(*colour_cuts)
                 .where(*astrometric_cuts)
                 .with_cte(quality_cte))

        if query_region:
            query = (query.where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                            Catalog.dec,
                                                            query_region[0],
                                                            query_region[1],
                                                            query_region[2])))

        return query


class MWM_CB_UVEX4_Carton(MWM_CB_UVEX_BaseCarton):
    """MWM Compact Binaries UV excess 4.

    Definition:

       Match Gaia - SWUVOTSSC (Swift UVOT catalog), apply color cuts to
       select objects with UV excess between main sequence and WD sequence.
       Remove objects with low parallax and proper motion accuracy to not
       get swamped by AGN and QSO.

    Pseudo-SQL:

        - Apply some general UV quality flags and work on subset of entries
        in UVOT catalog::
            - (UVW2_AB_MAG >-100) || (UVW1_AB_MAG >-100) || (UVM2_AB_MAG>-100)
            - !(equals(substring(UVM2_QUALITY_FLAG_BIN,8,9),"1") ||
                (equals(substring(UVM2_QUALITY_FLAG_BIN,7,8),"1") && UVM2_SIGNIF <10) ||
                (equals(substring(UVM2_QUALITY_FLAG_BIN,6,7),"1") && UVM2_SIGNIF <10) ||
                (equals(substring(UVM2_QUALITY_FLAG_BIN,3,4),"1") && UVM2_SIGNIF <10) ||
                equals(substring(UVM2_QUALITY_FLAG_BIN,2,3),"1"))

        - Match UV catalog with Gaia DR2 within 3 arcsec, keep the nearest
          match only, add distance and lower limit distance (columns r_est,
          r_lo; estimated and lower limit to estimated distance)
          from catalog gdr2_contrib/geometric_distance.

        - Definitions:
            - pm = sqrt(pmra*pmra+pmdec*pmdec)
            - pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm) +
                         (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            - logpmdpm = log10(pm/pm_error)
            - gmagab = phot_g_mean_mag-25.688+25.792
            - bmagab = phot_bp_mean_mag-25.351+25.386
            - rmagab = phot_rp_mean_mag-24.762+25.116
            - B_AB-R_AB = bmagab-rmagab
            - absg=phot_g_mean_mag-5*log10(r_est)+5
            - nuvg=nuv_mag-gmagab
            - fuvg=fuv_mag-gmagab
            - fuvnub=fuv_mag-nuv_mag
            - bp_rp=phot_bp_mean_mag-phot_rp_mean_mag

        - General cut:
            - Gaia: visibility_periods_used >5

        - Colour and magnitude cuts:
            - absg < (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
            - ((UVM2_AB_MAG - gmagab) < 2.25 ||
               ((UVM2_AB_MAG - gmagab) < 6 &&
                (UVM2_AB_MAG - gmagab) < 5.57377*(B_AB-R_AB)+0.2049))

        - Astrometric cuts:
            - !((logpmdpm <0.301 &&
                (parallax/parallax_error)>-1.4995*logpmdpm-4.05 &&
                (parallax/parallax_error)<1.4995*logpmdpm+4.05) ||
               (logpmdpm-0.301)*(logpmdpm-0.301)/(0.39794*0.39794) +
               (parallax/parallax_error)^2/(4.5*4.5)<=1)
            - r_lo <= 1500

    """

    name = 'mwm_cb_uvex4'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def build_query(self, version_id, query_region=None):

        uvw2_ab_mag = UVOT_SSC_1.uvw2_ab
        uvm2_ab_mag = UVOT_SSC_1.uvm2_ab
        uvw1_ab_mag = UVOT_SSC_1.uvw1_ab

        uvm2_quality_flag = UVOT_SSC_1.fuvm2
        uvm2_signif = UVOT_SSC_1.suvm2

        # uvXY_quality_flag is a 9-bit integer so character 8 (0-indexed) in
        # the binary string is the first bit, 7 is the second, etc.

        quality_flags = (
            (uvw2_ab_mag > -100) | (uvw1_ab_mag > -100) | (uvm2_ab_mag > -100),
            ~((uvm2_quality_flag.bin_and(2**0) > 0) |
              (uvm2_quality_flag.bin_and(2**1) > 0) & (uvm2_signif < 10) |
              (uvm2_quality_flag.bin_and(2**2) > 0) & (uvm2_signif < 10) |
              (uvm2_quality_flag.bin_and(2**5) > 0) & (uvm2_signif < 10) |
              (uvm2_quality_flag.bin_and(2**6) > 0))
        )

        quality_cte = (UVOT_SSC_1
                       .select(UVOT_SSC_1,
                               CatalogToUVOT_SSC_1.catalogid)
                       .join(CatalogToUVOT_SSC_1)
                       .where(CatalogToUVOT_SSC_1.version_id == version_id,
                              CatalogToUVOT_SSC_1.best >> True)
                       .where(*quality_flags)
                       .cte('quality_cte', materialized=True))

        colour_cuts = (
            colour_absg,
            ((quality_cte.c.uvm2_ab - gmagab) < 2.25) |
            (((quality_cte.c.uvm2_ab - gmagab) < 6) &
             ((quality_cte.c.uvm2_ab - gmagab) < 5.57377 * (B_AB_R_AB) + 0.2049)))

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_ra'),
                         Gaia_DR2.dec.alias('gaia_dec'),
                         Gaia_DR2.pmra,
                         Gaia_DR2.pmdec,
                         Gaia_DR2.pmra_error,
                         Gaia_DR2.pmdec_error,
                         Gaia_DR2.parallax,
                         Gaia_DR2.parallax_error,
                         Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         gmagab.alias('gmagab'),
                         BJ.r_est,
                         BJ.r_lo,
                         quality_cte.c.srcid,
                         quality_cte.c.obsid,
                         quality_cte.c.radeg.alias('uvot_ra'),
                         quality_cte.c.dedeg.alias('uvot_dec'),
                         quality_cte.c.suvm2,
                         quality_cte.c.uvm2_ab.alias('uvm2_ab'),
                         quality_cte.c.fuvm2.alias('fuvm2'))
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(BJ)
                 .join_from(Catalog, quality_cte,
                            on=(quality_cte.c.catalogid == Catalog.catalogid))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(Gaia_DR2.visibility_periods_used > 5)
                 .where(*colour_cuts)
                 .where(*astrometric_cuts)
                 .with_cte(quality_cte))

        if query_region:
            query = (query
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_UVEX5_Carton(MWM_CB_UVEX_BaseCarton):
    """MWM Compact Binaries UV excess 5.

    Definition:

       Match Gaia - GALEX, use both bands of GALEX (FUV and NUV), apply colour
       cuts to select objects with UV excess between main sequence and WD
       sequence. Remove objects with low parallax and proper motion accuracy to
       not get swamped by AGN and QSO.

    Pseudo-SQL:

        - Match Gaia DR2 with GALEX within 5 arcsec, keep the nearest match
          only, add distance and lower limit distance (columns r_est, r_lo;
          estimated and lower limit to estimated distance) from catalog
          gdr2_contrib/geometric_distance; use columns fuv_mag, nuv_mag,
          fuv_magerr, nuv_magerr from GALEX catalog.

        - Match with 2MASS.

        - Definitions:
            - pm = sqrt(pmra*pmra+pmdec*pmdec)
            - pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm) +
                         (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            - logpmdpm = log10(pm/pm_error)
            - gmr = phot_g_mean_mag - phot_rp_mean_mag
            - absg = phot_g_mean_mag-5*log10(r_est)+5
            - Main-sequence defined as polynomial:
                GMS = 0.00206868742*gmr**6 + 0.0401594518*gmr**5 -
                0.842512410*gmr**4 + 4.89384979*gmr**3 -
                12.3826637*gmr**2 + 17.0197205*gmr - 3.19987835

        - General cut:
            - Gaia: visibility_periods_used > 5

        - Colour and magnitude cuts:
            - Hmag < 15
            - fuv_magerr < 0.2 && nuv_magerr < 0.2 &&
              fuv_mag > -100 && nuv_mag > -100
            - abs(GMS-absg) <=0.5 && absg >=4.0866
            - r_est < 0.51 * pow(10, 0.2291*phot_g_mean_mag)

        - Astrometric cuts:
            - !((logpmdpm <0.301 &&
                (parallax/parallax_error)>-1.4995*logpmdpm-4.05 &&
                (parallax/parallax_error)<1.4995*logpmdpm+4.05) ||
               (logpmdpm-0.301)*(logpmdpm-0.301)/(0.39794*0.39794) +
               (parallax/parallax_error)^2/(4.5*4.5)<=1)
            - r_lo <= 1500

        This sequence yields 10,766 objects whose UV emission is though
        to arise from an unseen compact companion.

    """

    name = 'mwm_cb_uvex5'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def setup_transaction(self):

        self.database.execute_sql('SET LOCAL join_collapse_limit = 1;')
        super().setup_transaction()

    def build_query(self, version_id, query_region=None):

        bmr = Gaia_DR2.phot_bp_mean_mag - Gaia_DR2.phot_rp_mean_mag
        GMS = (0.00206868742 * fn.pow(bmr, 6) + 0.0401594518 * fn.pow(bmr, 5) -
               0.842512410 * fn.pow(bmr, 4) + 4.89384979 * fn.pow(bmr, 3) -
               12.3826637 * fn.pow(bmr, 2) + 17.0197205 * bmr - 3.19987835)

        colour_cuts = (TwoMassPSC.h_m < 15,
                       fn.abs(GMS - absg) <= 0.5,
                       absg >= 4.0866,
                       BJ.r_est < 0.51 * fn.pow(10, 0.2291 * Gaia_DR2.phot_g_mean_mag))

        query = (GUVCat
                 .select(CatalogToTIC_v8.catalogid,
                         Gaia_DR2.source_id,
                         Gaia_DR2.ra.alias('gaia_ra'),
                         Gaia_DR2.dec.alias('gaia_dec'),
                         Gaia_DR2.pmra,
                         Gaia_DR2.pmdec,
                         Gaia_DR2.pmra_error,
                         Gaia_DR2.pmdec_error,
                         Gaia_DR2.parallax,
                         Gaia_DR2.parallax_error,
                         Gaia_DR2.phot_g_mean_mag,
                         Gaia_DR2.phot_bp_mean_mag,
                         Gaia_DR2.phot_rp_mean_mag,
                         BJ.r_est,
                         BJ.r_lo,
                         GUVCat.objid.alias('guvcat_objid'),
                         GUVCat.ra.alias('guvcat_ra'),
                         GUVCat.dec.alias('guvcat_dec'),
                         GUVCat.nuv_mag,
                         GUVCat.fuv_mag,
                         GUVCat.nuv_magerr,
                         GUVCat.fuv_magerr,
                         GMS.alias('GMS'),
                         TwoMassPSC.ra,
                         TwoMassPSC.decl,
                         TwoMassPSC.h_m,
                         TwoMassPSC.pts_key)
                 .join(CatalogToGUVCat)
                 .join(CatalogToTIC_v8,
                       on=(CatalogToGUVCat.catalogid == CatalogToTIC_v8.catalogid))
                 .join(TIC_v8)
                 .join(TwoMassPSC)
                 .join_from(TIC_v8, Gaia_DR2)
                 .join(BJ)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True)
                 .where(CatalogToGUVCat.version_id == version_id,
                        CatalogToGUVCat.best >> True)
                 .where(Gaia_DR2.visibility_periods_used > 5)
                 .where(fuv_magerr < 0.2,
                        nuv_magerr < 0.2,
                        fuv_mag > -100,
                        nuv_mag > -100)
                 .where(*colour_cuts)
                 .where(*astrometric_cuts))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
