#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2023-04-28
# @Filename: mwm_magcloud_rgb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import dustmaps.sfd
import numpy
import pandas
from astropy.coordinates import SkyCoord
from astropy.units import mas, yr
from dustmaps.sfd import SFDQuery
from gala.coordinates import MagellanicStreamNidever08
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import CatalogToGaia_DR3, Gaia_DR3

from target_selection import log
from target_selection.cartons import BaseCarton
from target_selection.cartons.mwm_magcloud_agb import roi_cut
from target_selection.cartons.tools import create_table_as


class MWM_MagCloud_RGB_BOSS(BaseCarton):
    """MWM Magellanic clouds RGBs.

    Definition:

    Select AGB targets in the Magellanic Clouds using Gaia DR3 photometry and astrometry.
    Selection based on color cuts, proper motion cuts and parallax cuts. See wiki
    for David Nidever's code.

    """

    name = "mwm_magcloud_rgb_boss"
    mapper = "MWM"
    category = "science"
    program = "mwm_magcloud"
    instrument = "BOSS"
    priority = 2819
    cadence = "bright_1x1"
    can_offset = True

    def build_query(self, version_id, query_region=None):
        # Parallax cut
        parallax_cut = ~(
            (Gaia_DR3.parallax > 0) & ((Gaia_DR3.parallax + 0.025) / Gaia_DR3.parallax_error > 5)
        )

        # Rough cuts. Just to reduce the number of rows returned for post-process.
        bp_rp = Gaia_DR3.phot_bp_mean_mag - Gaia_DR3.phot_rp_mean_mag
        colour_cut = (Gaia_DR3.phot_g_mean_mag <= self.parameters["gaiag_rough"]) & (
            bp_rp >= self.parameters["bp_rp_rough"]
        )

        pm_cut = (
            (Gaia_DR3.pmra >= self.parameters["pmra_rough"][0])
            & (Gaia_DR3.pmra <= self.parameters["pmra_rough"][1])
            & (Gaia_DR3.pmdec >= self.parameters["pmdec_rough"][0])
            & (Gaia_DR3.pmdec <= self.parameters["pmdec_rough"][1])
        )

        ra, dec, radius = self.parameters["astro_rough"]
        astro_cut = fn.q3c_radial_query(Gaia_DR3.ra, Gaia_DR3.dec, ra, dec, radius)

        log.info("Running Q3C query ...")
        tmp_q3c, _ = create_table_as(
            Gaia_DR3.select(Gaia_DR3.source_id).where(astro_cut),
            "tmp_q3c",
            temporary=True,
            database=self.database,
        )

        query = (
            CatalogToGaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
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
            )
            .join(Gaia_DR3)
            .join(tmp_q3c, on=(tmp_q3c.c.source_id == Gaia_DR3.source_id))
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                Gaia_DR3.phot_g_mean_mag < self.parameters["g_lim"],
                parallax_cut,
                colour_cut,
                pm_cut,
            )
        )

        return query

    def post_process(self, model, **kwargs):
        """Runs post-process."""

        data = pandas.read_sql(f"SELECT * from {self.path}", self.database)

        # Calculate Magellanic Stream coordinates
        coords = SkyCoord(
            ra=data.ra.values,
            dec=data.dec.values,
            unit="deg",
            frame="icrs",
            pm_ra_cosdec=data.pmra.values * mas / yr,
            pm_dec=data.pmdec.values * mas / yr,
        )
        mcoo = coords.transform_to(MagellanicStreamNidever08)

        data["pml_ms"] = mcoo.pm_L_cosB.value
        data["pmb_ms"] = mcoo.pm_B.value
        data["mlat"] = mcoo.B.value
        data["mlon"] = mcoo.L.value

        # Calculate LMC and SMC radius
        lmcrad = numpy.sqrt((data["mlon"] + 0.6) ** 2 + (data["mlat"] - 3.6) ** 2)
        smcrad = numpy.sqrt((data["mlon"] + 15.53) ** 2 + (data["mlat"] + 11.58) ** 2)
        data["lmcrad"] = lmcrad
        data["smcrad"] = smcrad

        # Proper motion cut
        pmdist = numpy.sqrt((data["pml_ms"] - 1.8) ** 2 + (data["pmb_ms"] - 0.40) ** 2)
        gdpm = pmdist < self.parameters["pmdist"]
        data = data.loc[gdpm]

        # Make sure none of the needed quantities are NaN
        good = (
            numpy.isfinite(data["pmra"])
            & numpy.isfinite(data["pmdec"])
            & numpy.isfinite(data["parallax"])
            & numpy.isfinite(data["phot_g_mean_mag"])
            & numpy.isfinite(data["phot_bp_mean_mag"])
            & numpy.isfinite(data["phot_rp_mean_mag"])
        )
        data = data[good]

        # Get extinction values from SFD
        dustmaps.sfd.fetch()
        coo = SkyCoord(ra=data["ra"], dec=data["dec"], unit="degree", frame="icrs")
        sfd = SFDQuery()
        ebv = sfd(coo)

        # For high-extinction regions in the inner LMC/SMC SFD is not that reliable
        # but we can "rescale" ebv to get something reasonable
        # Selection for inner LMC/SMC with high extinction
        (bd,) = numpy.where(((data["lmcrad"] <= 4.5) | (data["smcrad"] <= 2.5)) & (ebv >= 0.1))
        ebv[bd] = ebv[bd] * 0.07  # rescale

        # Get dereddened magnitudes
        Av = 0.86 * 3.1 * ebv
        ag = 0.8509139481578785 * Av
        abp = 1.0901361069987485 * Av
        arp = 0.5894794093231802 * Av
        gmag0 = data["phot_g_mean_mag"] - ag
        bp0 = data["phot_bp_mean_mag"] - abp
        rp0 = data["phot_rp_mean_mag"] - arp

        # Apply CMD cut to dereddened magnitudes
        bprpcut = self.parameters["bprpcut"]
        gcut = self.parameters["gcut"]
        _, cutind = roi_cut(bprpcut, gcut, bp0 - rp0, gmag0)
        data = data.iloc[cutind]

        # Final magnitude cut
        (gcut,) = numpy.where(data["phot_g_mean_mag"] <= 17.5)
        data = data.iloc[gcut]

        valid_cids = data.catalogid.values

        # This way seems faster than updating from a list of values.
        values_cids = ",".join(f"({vc})" for vc in valid_cids)
        self.database.execute_sql(
            "CREATE TEMP TABLE valid_cids AS SELECT * "
            f"FROM (VALUES {values_cids}) "
            "AS t (catalogid)"
        )
        self.database.execute_sql("CREATE INDEX ON valid_cids (catalogid);")
        self.database.execute_sql(f"UPDATE {self.path} SET selected = false")
        self.database.execute_sql(
            f"UPDATE {self.path} SET selected = true "
            "FROM valid_cids vc "
            f"WHERE {self.path}.catalogid = vc.catalogid"
        )

        return super().post_process(model, **kwargs)
