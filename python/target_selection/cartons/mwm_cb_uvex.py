#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-20
# @Filename: mwm_cb_uvex.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (XMM_OM_SUSS_4_1, Catalog,
                                             CatalogToGUVCat, CatalogToTIC_v8,
                                             CatalogToXMM_OM_SUSS_4_1,
                                             Gaia_DR2,
                                             GeometricDistances_Gaia_DR2,
                                             GUVCat, TIC_v8)

from . import BaseCarton


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

colour_absg = ((absg < 4.5457 * bp_rp + 9.9) &
               ((absg > 4.09) | (absg > 4.5457 * bp_rp + 4.0457)))

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


class MWM_CB_UVEX1_Carton(BaseCarton):
    """MWM Compact Binaries UV excess 1.

    Definition:

       Match Gaia - GALEX, use NUV channel of GALEX only, apply color cuts to
       select objects with UV excess between main sequence and WD sequence.
       Remove objects with low parallax and proper motion accuracy to not get
       swamped by AGN and QSO.

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
            - absg < 4.5457*bp_rp+9.9 &&
                (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
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
    program = 'CB'

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

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR2.source_id,
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
                         GUVCat.nuv_mag,
                         GUVCat.fuv_mag,
                         GUVCat.nuv_magerr,
                         GUVCat.fuv_magerr)
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(BJ)
                 .join_from(Catalog, CatalogToGUVCat)
                 .join(GUVCat)
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


class MWM_CB_UVEX2_Carton(BaseCarton):
    """MWM Compact Binaries UV excess 2.

    Definition:

       Match Gaia DR2 with GALEX within 5 arcsec, keep the nearest match only,
       add distance and lower limit distance (columns r_est, r_lo; estimated
       and lower limit to estimated distance) from catalog
       gdr2_contrib/geometric_distance

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
            - absg < 4.5457*bp_rp+9.9 &&
                (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
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
    program = 'CB'

    def build_query(self, version_id, query_region=None):

        colour_cuts = (nuv_magerr < 0.2,
                       nuv_mag > -100,
                       ((fuv_magerr >= 0.2) | (fuv_mag < -100)),
                       colour_absg,
                       (((nuv_mag - gmagab) < 2.25) |
                        (((nuv_mag - gmagab) < 6.725 * B_AB_R_AB - 1.735) &
                         ((nuv_mag - gmagab) < -0.983 * B_AB_R_AB + 8.24))))

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR2.source_id,
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
                         GUVCat.nuv_mag,
                         GUVCat.fuv_mag,
                         GUVCat.nuv_magerr,
                         GUVCat.fuv_magerr)
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(Gaia_DR2)
                 .join(BJ)
                 .join_from(Catalog, CatalogToGUVCat)
                 .join(GUVCat)
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


class MWM_CB_UVEX3_Carton(BaseCarton):
    """MWM Compact Binaries UV excess 3.

    Definition:

       Match Gaia DR2 with GALEX within 5 arcsec, keep the nearest match only,
       add distance and lower limit distance (columns r_est, r_lo; estimated
       and lower limit to estimated distance) from catalog
       gdr2_contrib/geometric_distance

    Pseudo-SQL:

        - Apply some general UV quality flags and work on subset of entries
          in OM-SUSS:
            - (UVW2_AB_MAG >-100) || (UVW1_AB_MAG >-100) || (UVM2_AB_MAG>-100)
            - !(equals(substring(UVW2_QUALITY_FLAG_ST,0,1),"T") ||
                equals(substring(UVM2_QUALITY_FLAG_ST,0,1),"T") ||
                equals(substring(UVW1_QUALITY_FLAG_ST,0,1),"T") ||
                (equals(substring(UVW2_QUALITY_FLAG_ST,1,2),"T") &&
                 UVW2_SIGNIF <10) ||
                (equals(substring(UVM2_QUALITY_FLAG_ST,1,2),"T") &&
                 UVM2_SIGNIF <10) ||
                (equals(substring(UVW1_QUALITY_FLAG_ST,1,2),"T") &&
                 UVW1_SIGNIF <10) ||
                (equals(substring(UVW2_QUALITY_FLAG_ST,2,3) ,"T") &&
                 UVW2_SIGNIF <10) ||
                (equals(substring(UVM2_QUALITY_FLAG_ST,2,3),"T") &&
                 UVM2_SIGNIF <10) ||
                (equals(substring(UVW1_QUALITY_FLAG_ST,2,3),"T") &&
                 UVW1_SIGNIF <10) ||
                (equals(substring(UVW2_QUALITY_FLAG_ST,6,7) ,"T")) ||
                equals(substring(UVM2_QUALITY_FLAG_ST,6,7),"T") ||
                equals(substring(UVW1_QUALITY_FLAG_ST,6,7),"T") ||
                equals(substring(UVW2_QUALITY_FLAG_ST,7,8),"T") ||
                equals(substring(UVM2_QUALITY_FLAG_ST,7,8),"T") ||
                equals(substring(UVW1_QUALITY_FLAG_ST,7,8),"T") ||
                equals(substring(UVW2_QUALITY_FLAG_ST,8,9),"T") ||
                equals(substring(UVM2_QUALITY_FLAG_ST,8,9),"T") ||
                equals(substring(UVW1_QUALITY_FLAG_ST,8,9),"T"))

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
            - absg < 4.5457*bp_rp+9.9 &&
                (absg > 4.09 || absg > 4.5457*bp_rp + 4.0457)  [colour_abs]
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
    program = 'CB'

    def build_query(self, version_id, query_region=None):

        uvw2_ab_mag = XMM_OM_SUSS_4_1.uvw2_ab_mag
        uvw1_ab_mag = XMM_OM_SUSS_4_1.uvw1_ab_mag
        uvm2_ab_mag = XMM_OM_SUSS_4_1.uvm2_ab_mag

        uvw2_quality_flag_st = XMM_OM_SUSS_4_1.uvw2_quality_flag_st
        uvm2_quality_flag_st = XMM_OM_SUSS_4_1.uvm2_quality_flag_st
        uvw1_quality_flag_st = XMM_OM_SUSS_4_1.uvw1_quality_flag_st

        uvw2_signif = XMM_OM_SUSS_4_1.uvw2_signif
        uvm2_signif = XMM_OM_SUSS_4_1.uvm2_signif
        uvw1_signif = XMM_OM_SUSS_4_1.uvw1_signif

        quality_flags = (
            (uvw2_ab_mag > -100) | (uvw1_ab_mag > -100) | (uvm2_ab_mag > -100),
            ~((fn.substr(uvw2_quality_flag_st, 1, 1) == 'T') |
              (fn.substr(uvm2_quality_flag_st, 1, 1) == 'T') |
              (fn.substr(uvw1_quality_flag_st, 1, 1) == 'T') |
              ((fn.substr(uvw2_quality_flag_st, 2, 1) == 'T') & (uvw2_signif < 10)) |
              ((fn.substr(uvm2_quality_flag_st, 2, 1) == 'T') & (uvm2_signif < 10)) |
              ((fn.substr(uvw1_quality_flag_st, 2, 1) == 'T') & (uvw1_signif < 10)) |
              ((fn.substr(uvw2_quality_flag_st, 3, 1) == 'T') & (uvw2_signif < 10)) |
              ((fn.substr(uvm2_quality_flag_st, 3, 1) == 'T') & (uvm2_signif < 10)) |
              ((fn.substr(uvw1_quality_flag_st, 3, 1) == 'T') & (uvw1_signif < 10)) |
              ((fn.substr(uvw2_quality_flag_st, 7, 1) == 'T') |
               (fn.substr(uvm2_quality_flag_st, 7, 1) == 'T') |
               (fn.substr(uvw1_quality_flag_st, 7, 1) == 'T') |
               (fn.substr(uvw2_quality_flag_st, 8, 1) == 'T') |
               (fn.substr(uvm2_quality_flag_st, 8, 1) == 'T') |
               (fn.substr(uvw1_quality_flag_st, 8, 1) == 'T') |
               (fn.substr(uvw2_quality_flag_st, 9, 1) == 'T') |
               (fn.substr(uvm2_quality_flag_st, 9, 1) == 'T') |
               (fn.substr(uvw1_quality_flag_st, 9, 1) == 'T')))
        )

        quality_cte = (XMM_OM_SUSS_4_1
                       .select(XMM_OM_SUSS_4_1,
                               CatalogToXMM_OM_SUSS_4_1.catalogid)
                       .join(CatalogToXMM_OM_SUSS_4_1)
                       .where(CatalogToXMM_OM_SUSS_4_1.version_id == version_id,
                              CatalogToXMM_OM_SUSS_4_1.best >> True)
                       .where(*quality_flags)
                       .cte('quality_cte', materialized=True))

        colour_cuts = (colour_absg,
                       ((quality_cte.c.uvm2_ab_mag - gmagab) < 2.25) |
                       (((quality_cte.c.uvm2_ab_mag - gmagab) < 6) &
                        ((quality_cte.c.uvm2_ab_mag - gmagab) < 5.57377 * (B_AB_R_AB) + 0.2049)))

        query = (Catalog
                 .select(Catalog.catalogid,
                         Gaia_DR2.source_id,
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
                         quality_cte.c.uvm2_ab_mag.alias('uvm2_ab_mag'),
                         quality_cte.c.uvw2_ab_mag.alias('uvw2_ab_mag'),
                         quality_cte.c.uvw1_ab_mag.alias('uvw1_ab_mag'),
                         quality_cte.c.uvm2_quality_flag_st.alias('uvm2_quality_flag_st'),
                         quality_cte.c.uvw2_quality_flag_st.alias('uvw2_quality_flag_st'),
                         quality_cte.c.uvw1_quality_flag_st.alias('uvw1_quality_flag_st'))
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
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
