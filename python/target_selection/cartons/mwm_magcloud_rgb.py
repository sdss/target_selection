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

from sdssdb.peewee.sdss5db.catalogdb import (CatalogToGaia_DR3,
                                             CatalogToTwoMassPSC, Gaia_DR3,
                                             TwoMassPSC)

from target_selection.cartons import BaseCarton


def roi_cut(xcut, ycut, x, y):
    """Use cuts in a 2D plane to select points from arrays.

    Parameters
    ----------
    xcut : numpy array
        Array of x-values for the cut.
    ycut : numpy array
        Array of y-values for the cut.
    x : numpy array or list
        Array of x-values that should be cut.
    y : numpy array or list
        Array of y-values that should be cut.

    Returns
    -------
    ind : numpy array
       The indices of values OUTSIDE the cut.
    cutind :
       The indices of values INSIDE the cut.

    """

    from matplotlib.path import Path

    tupVerts = list(zip(xcut, ycut))

    points = numpy.vstack((x, y)).T

    p = Path(tupVerts)  # make a polygon
    inside = p.contains_points(points)

    ind, = numpy.where(~inside)
    cutind, = numpy.where(inside)

    return ind, cutind


class MWM_MagCloud_RGB_Base(BaseCarton):
    """MWM Magellanic clouds RGBs.

    Definition:

    Select AGB targets in the Magellanic Clouds using Gaia DR3 photometry and astrometry.
    Selection based on color cuts, proper motion cuts and parallax cuts. See wiki
    for David Nidever's code.

    """

    mapper = 'MWM'
    category = 'science'
    program = 'mwm_magcloud_agb'
    can_offset = True

    def build_query(self, version_id, query_region=None):

        # Parallax cut
        parallax_cut = ~((Gaia_DR3.parallax > 0) &
                         (Gaia_DR3.parallax / Gaia_DR3.parallax_error > 5))

        # Rough cuts. Just to reduce the number of rows returned for post-process.
        bp_rp = Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag
        colour_cut = ((Gaia_DR3.phot_g_mean_mag <= self.parameters['gaiag_rough']) &
                      (bp_rp >= self.parameters['bp_rp_rough']))

        pm_cut = ((Gaia_DR3.pmra >= self.parameters['pmra_rough'][0]) &
                  (Gaia_DR3.pmra <= self.parameters['pmra_rough'][1]) &
                  (Gaia_DR3.pmdec >= self.parameters['pmdec_rough'][0]) &
                  (Gaia_DR3.pmdec <= self.parameters['pmdec_rough'][1]))

        ra, dec, radius = self.parameters['astro_rough']
        # astro_cut = fn.q3c_radial_query(Gaia_DR3.ra, Gaia_DR3.dec, ra, dec, radius)
        astro_cut = Gaia_DR3.dec < -40

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
                         TwoMassPSC.h_m,
                         TwoMassPSC.j_m,
                         TwoMassPSC.k_m)
                 .join(Gaia_DR3)
                 .join_from(CatalogToGaia_DR3, CatalogToTwoMassPSC,
                            on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid))
                 .join(TwoMassPSC)
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        CatalogToTwoMassPSC.best >> True,
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

        mlon_a = self.parameters['mlon_a']
        mlon_b = self.parameters['mlon_b']

        # CMD cut
        # AGB-O and AGB-C
        # Apply separate cuts to LMC and SMC

        # -- LMC --
        jkcut_lmc = numpy.array(self.parameters['jkcut_lmc'])
        hcut_lmc = numpy.array(self.parameters['hcut_lmc'])
        glmc = data['mlat'] > (mlon_a * data['mlon'] - mlon_b)
        data_lmc = data.loc[glmc]
        _, cutind = roi_cut(jkcut_lmc,
                            hcut_lmc,
                            data_lmc['j_m'] - data_lmc['k_m'],
                            data_lmc['h_m'])
        data_lmc = data_lmc.iloc[cutind]

        # -- SMC --
        jkcut_smc = jkcut_lmc + self.parameters['jkcut_smc_off']
        hcut_smc = hcut_lmc + self.parameters['hcut_smc_off']
        gsmc = data['mlat'] <= (mlon_a * data['mlon'] - mlon_b)
        data_smc = data.loc[gsmc]
        _, cutind = roi_cut(jkcut_smc,
                            hcut_smc,
                            data_smc['j_m'] - data_smc['k_m'],
                            data_smc['h_m'])

        data_smc = data_smc.iloc[cutind]

        valid_cids = numpy.hstack([data_lmc.catalogid.values,
                                   data_smc.catalogid.values])

        (model
         .delete()
         .where(model.catalogid.not_in(valid_cids.tolist()))
         .execute())

        return super().post_process(model, **kwargs)


class MWM_MagCloud_RGB_APOGEE(MWM_MagCloud_RGB_Base):
    """MWM Magellanic clouds RGBs. APOGEE carton."""

    name = 'mwm_magcloud_agb_apogee'
    instrument = 'APOGEE'
    priority = 2818

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(TwoMassPSC.h_m < self.parameters['h_lim'])

        return query


class MWM_MagCloud_RGB_BOSS(MWM_MagCloud_RGB_Base):
    """MWM Magellanic clouds RGBs. BOSS carton."""

    name = 'mwm_magcloud_agb_boss'
    instrument = 'BOSS'
    priority = 2819

    def build_query(self, version_id, query_region=None):

        query = super().build_query(version_id, query_region)
        query = query.where(Gaia_DR3 < self.parameters['g_lim'])

        return query
