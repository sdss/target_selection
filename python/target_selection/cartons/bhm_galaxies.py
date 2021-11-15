#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2021-11-02
# @Filename: bhm_galaxies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import JOIN
from peewee import fn
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    CatalogToLegacy_Survey_DR8,
    Legacy_Survey_DR8,
)
from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import AB2nMgy

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms
#
# This module provides the following BHM cartons:
# bhm_colr_galaxies_lsdr8
#

'''
'''


class BhmColrGalaxiesLsdr8Carton(BaseCarton):
    '''
    '''

    name = 'bhm_colr_galaxies_lsdr8'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_filler'
    tile = False
    cadence = 'dark_1x1'
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()
        ls = Legacy_Survey_DR8.alias()

        # set the Carton priority+values here - read from yaml
        priority = peewee.Value(int(self.parameters.get('priority', 10000)))
        value = peewee.Value(self.parameters.get('value', 0.0)).cast('float')
        inertial = peewee.Value(True)
        instrument = peewee.Value(self.instrument)
        #cadence = peewee.Value(self.parameters.get('cadence', self.cadence))

        dered_flux_z_min = AB2nMgy(self.parameters['dered_mag_z_max'])
        dered_fiberflux_z_min = AB2nMgy(self.parameters['dered_fibermag_z_max'])
        fiberflux_z_min = AB2nMgy(self.parameters['fibermag_z_max'])
        fiberflux_z_max = AB2nMgy(self.parameters['fibermag_z_min'])
        fiberflux_r_max = AB2nMgy(self.parameters['fibermag_r_min'])

        fiberflux_z_min_for_cadence1 = AB2nMgy(self.parameters['fibermag_z_for_cadence1'])
        fiberflux_z_min_for_cadence2 = AB2nMgy(self.parameters['fibermag_z_for_cadence2'])

        # compute transformed SDSS mags uniformly
        # transform the legacysurvey grz into sdss fiber2mag griz

        # extract coeffs from fit logs via:
        # awk 'BEGIN {print("coeffs = {")} /POLYFIT/{ if($3~/sdss_psfmag/){pe="p"} else if ($3~/sdss_fiber2mag/){pe="e"} else{pe="error"}; printf("\"%s%d_%s\": %s,\n", substr($3,length($3)), $8, pe, $10)} END {print("}")}'  bhm_spiders_clusters_lsdr8/lsdr8_fibermag_to_sdss_fiber2mag_?_results.log   # noqa
        coeffs = {
            "g2_e": -0.897719,
            "g1_e": 2.298300,
            "g0_e": -1.019299,
            "i2_e": -0.950114,
            "i1_e": 0.981972,
            "i0_e": -0.261645,
            "r2_e": -0.201741,
            "r1_e": 0.697128,
            "r0_e": -0.120926,
            "z2_e": -1.424312,
            "z1_e": 2.415301,
            "z0_e": -0.677163,
        }

        nMgy_min = 1e-3  # equiv to AB=30
        # extended - start from ls8 fiberfluxes
        g0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g)))
        r0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r)))
        z0_e = (22.5 - 2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_z)))
        g_r_e = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_g) /
                                      peewee.fn.greatest(nMgy_min, ls.fiberflux_r)))
        r_z_e = (-2.5 * peewee.fn.log(peewee.fn.greatest(nMgy_min, ls.fiberflux_r) /
                                      peewee.fn.greatest(nMgy_min, ls.fiberflux_z)))

        g_e = (g0_e + coeffs['g0_e'] + coeffs['g1_e'] * g_r_e + coeffs['g2_e'] * g_r_e * g_r_e)
        r_e = (r0_e + coeffs['r0_e'] + coeffs['r1_e'] * g_r_e + coeffs['r2_e'] * g_r_e * g_r_e)
        i_e = (r0_e + coeffs['i0_e'] + coeffs['i1_e'] * r_z_e + coeffs['i2_e'] * r_z_e * r_z_e)
        z_e = (z0_e + coeffs['z0_e'] + coeffs['z1_e'] * r_z_e + coeffs['z2_e'] * r_z_e * r_z_e)

        # validity checks
        valid = (g0_e.between(0.1, 29.9) &
                 r0_e.between(0.1, 29.9) &
                 z0_e.between(0.1, 29.9))

        opt_prov = peewee.Case(None, ((valid, 'sdss_fiber2mag_from_lsdr8'),), 'undefined')
        magnitude_g = peewee.Case(None, ((valid, g_e),), 'NaN')
        magnitude_r = peewee.Case(None, ((valid, r_e),), 'NaN')
        magnitude_i = peewee.Case(None, ((valid, i_e),), 'NaN')
        magnitude_z = peewee.Case(None, ((valid, z_e),), 'NaN')
        magnitude_gaia_g = peewee.Case(
            None,
            ((ls.gaia_phot_g_mean_mag.between(0.1, 29.9), ls.gaia_phot_g_mean_mag),),
            'NaN')
        magnitude_gaia_bp = peewee.Case(
            None,
            ((ls.gaia_phot_bp_mean_mag.between(0.1, 29.9), ls.gaia_phot_bp_mean_mag),),
            'NaN')
        magnitude_gaia_rp = peewee.Case(
            None,
            ((ls.gaia_phot_rp_mean_mag.between(0.1, 29.9), ls.gaia_phot_rp_mean_mag),),
            'NaN')

        opt_prov = peewee.Case(None, ((valid, 'sdss_fiber2mag_from_lsdr8'),), 'undefined')

        cadence = peewee.Case(
            None,
            (
                (ls.fiberflux_z > fiberflux_z_min_for_cadence1, self.parameters['cadence1']),
                (ls.fiberflux_z > fiberflux_z_min_for_cadence2, self.parameters['cadence2']),
            ),
            self.parameters['cadence3'])

        query = (
            c.select(
                c.catalogid.alias('catalogid'),
                ls.ls_id.alias('ls_id'),  # extra
                c.ra.alias('ra'),  # extra
                c.dec.alias('dec'),  # extra
                priority.alias('priority'),
                value.alias('value'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                opt_prov.alias('optical_prov'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                magnitude_gaia_g.alias('gaia_g'),
                magnitude_gaia_bp.alias('bp'),
                magnitude_gaia_rp.alias('rp'),
                inertial.alias('inertial'),
                g0_e.alias('ls8_fibermag_g'),  # extra
                r0_e.alias('ls8_fibermag_r'),  # extra
                z0_e.alias('ls8_fibermag_z'),  # extra
                ls.flux_g.alias('ls8_flux_g'),  # extra
                ls.flux_r.alias('ls8_flux_r'),  # extra
                ls.flux_z.alias('ls8_flux_z'),  # extra
                ls.ebv.alias('ls8_ebv'),  # extra
                ls.mw_transmission_z.alias('ls8_mw_transmission_z'),  # extra
            )
            .join(c2ls)
            .join(ls)
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                ls.type != 'PSF',
                ls.parallax == 0.0,
                ls.flux_z > dered_flux_z_min * ls.mw_transmission_z,
                ls.fiberflux_z > fiberflux_z_min,
                ls.fiberflux_z > dered_fiberflux_z_min * ls.mw_transmission_z,
                ls.fiberflux_r < fiberflux_r_max,
                ls.fiberflux_z < fiberflux_z_max,
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters['gaia_g_mag_limit'])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters['gaia_rp_mag_limit'])),
            )
            .distinct(c.catalogid)
        )

        # query_region[0] is ra of center of the region, degrees
        # query_region[1] is dec of center of the region, degrees
        # query_region[2] is radius of the region, degrees
        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query
