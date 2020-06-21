#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-20
# @Filename: mwm_cb_uvex.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGUVCat,
                                             CatalogToTIC_v8, Gaia_DR2,
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
    """MWM Halo Best & Brightest.

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
                       ((fuvg < 6.08) | ((fuvg < 11.82 * (B_AB - R_AB) + 2.58) &
                                         (fuvg < -0.79 * (B_AB - R_AB) + 9.21))))

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
