#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2023-04-28
# @Filename: mwm_magcloud_rgb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import numpy
import pandas
from astropy.coordinates import SkyCoord
from astropy.units import mas, yr
from gala.coordinates import MagellanicStreamNidever08
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (CatalogToGaia_DR3, Gaia_DR3,
                                             Gaia_dr3_astrophysical_parameters,
                                             Gaia_Stellar_Parameters)

from target_selection.cartons import BaseCarton
from target_selection.cartons.mwm_magcloud_agb import roi_cut


class MWM_MagCloud_RGB_BOSS(BaseCarton):
    """MWM Magellanic clouds RGBs.

    Definition:

    Select AGB targets in the Magellanic Clouds using Gaia DR3 photometry and astrometry.
    Selection based on color cuts, proper motion cuts and parallax cuts. See wiki
    for David Nidever's code.

    """

    name = 'mwm_magcloud_rgb_boss'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_magcloud_rgb'
    instrument = 'BOSS'
    priority = 2819
    can_offset = True

    def build_query(self, version_id, query_region=None):

        # Parallax cut
        parallax_cut = ~((Gaia_DR3.parallax > 0) &
                         ((Gaia_DR3.parallax + 0.025) / Gaia_DR3.parallax_error > 5))

        # Rough cuts. Just to reduce the number of rows returned for post-process.
        bp_rp = Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag
        colour_cut = ((Gaia_DR3.phot_g_mean_mag <= self.parameters['gaiag_rough']) &
                      (bp_rp >= self.parameters['bp_rp_rough']))

        pm_cut = ((Gaia_DR3.pmra >= self.parameters['pmra_rough'][0]) &
                  (Gaia_DR3.pmra <= self.parameters['pmra_rough'][1]) &
                  (Gaia_DR3.pmdec >= self.parameters['pmdec_rough'][0]) &
                  (Gaia_DR3.pmdec <= self.parameters['pmdec_rough'][1]))

        ra, dec, radius = self.parameters['astro_rough']
        astro_cut = fn.q3c_radial_query(Gaia_DR3.ra, Gaia_DR3.dec, ra, dec, radius)

        query = (CatalogToGaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.ra,
                         Gaia_DR3.dec,
                         Gaia_DR3.l,
                         Gaia_DR3.b,
                         Gaia_DR3.pmra,
                         Gaia_DR3.pmdec,
                         Gaia_DR3.parallax,
                         Gaia_DR3.parallax_error,
                         Gaia_DR3.phot_g_mean_mag,
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_rp_mean_mag,
                         Gaia_Stellar_Parameters.chi2_opt,
                         Gaia_Stellar_Parameters.feh_confidence,
                         Gaia_Stellar_Parameters.gdr3_source_id,
                         Gaia_Stellar_Parameters.ln_prior,
                         Gaia_Stellar_Parameters.logg_confidence,
                         Gaia_Stellar_Parameters.quality_flags,
                         Gaia_Stellar_Parameters.teff_confidence,
                         Gaia_Stellar_Parameters.stellar_params_est_teff,
                         Gaia_Stellar_Parameters.stellar_params_est_fe_h,
                         Gaia_Stellar_Parameters.stellar_params_est_logg,
                         Gaia_Stellar_Parameters.stellar_params_est_e,
                         Gaia_Stellar_Parameters.stellar_params_est_parallax,
                         Gaia_Stellar_Parameters.stellar_params_err_teff,
                         Gaia_Stellar_Parameters.stellar_params_err_fe_h,
                         Gaia_Stellar_Parameters.stellar_params_err_logg,
                         Gaia_Stellar_Parameters.stellar_params_err_e,
                         Gaia_Stellar_Parameters.stellar_params_err_parallax,
                         Gaia_dr3_astrophysical_parameters.ag_gspphot,
                         Gaia_dr3_astrophysical_parameters.abp_gspphot,
                         Gaia_dr3_astrophysical_parameters.arp_gspphot)
                 .join(Gaia_DR3)
                 .join_from(Gaia_DR3, Gaia_dr3_astrophysical_parameters,
                            on=(Gaia_DR3.source_id == Gaia_dr3_astrophysical_parameters.source_id))
                 .join_from(Gaia_DR3, Gaia_Stellar_Parameters)
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_g_mean_mag < self.parameters['g_lim'],
                        parallax_cut,
                        colour_cut,
                        pm_cut,
                        astro_cut))

        return query

    def post_process(self, model, **kwargs):
        """Runs post-process."""

        data = pandas.read_sql(f'SELECT * from {self.path}', self.database)

        # Calculate Magellanic Stream coordinates
        coords = SkyCoord(ra=data.ra,
                          dec=data.dec,
                          unit='deg',
                          frame='icrs',
                          pm_ra_cosdec=data.pmra * mas / yr,
                          pm_dec=data.pmdec * mas / yr)
        mcoo = coords.transform_to(MagellanicStreamNidever08)

        data['pml_ms'] = mcoo.pm_L_cosB.value
        data['pmb_ms'] = mcoo.pm_B.value
        data['mlat'] = mcoo.B.value
        data['mlon'] = mcoo.L.value

        # Proper motion cut.
        pmdist = numpy.sqrt((data['pml_ms'] - 1.8)**2 + (data['pmb_ms'] - 0.40)**2)
        gdpm = pmdist < self.parameters['pmdist']
        data = data.loc[gdpm]

        # Make sure none of the needed quantities are NaN
        good = (numpy.isfinite(data['pmra']) &
                numpy.isfinite(data['pmdec']) &
                numpy.isfinite(data['parallax']) &
                numpy.isfinite(data['phot_g_mean_mag']) &
                numpy.isfinite(data['phot_bp_mean_mag']) &
                numpy.isfinite(data['phot_rp_mean_mag']) &
                numpy.isfinite(data['ag_gspphot']) &
                numpy.isfinite(data['abp_gspphot']) &
                numpy.isfinite(data['arp_gspphot']))
        data = data[good]

        # Get dereddened magnitudes
        gmag0 = data['phot_g_mean_mag'] - data['ag_gspphot']
        bp0 = data['phot_bp_mean_mag'] - data['abp_gspphot']
        rp0 = data['phot_rp_mean_mag'] - data['arp_gspphot']

        # Apply CMD cut
        bprpcut = self.parameters['bprpcut']
        gcut = self.parameters['gcut']
        _, cutind = roi_cut(bprpcut, gcut, bp0 - rp0, gmag0)
        data = data.iloc[cutind]

        valid_cids = data.catalogid.values

        (model
         .update({model.selected: False})
         .where(model.catalogid.not_in(valid_cids.tolist()))
         .execute())

        return super().post_process(model, **kwargs)
