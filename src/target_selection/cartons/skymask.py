#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-24
# @Filename: masking.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# flake8: noqa
# isort: skip_file

import os
import numpy as np
import healpy as hp
from pymoc import MOC
import pymangle as mangle
from astropy.io import fits
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u


known_masktypes = [
    "mangle",
    "hpxmoc",
    "circles",
]
known_senses = [
    "include",
    "exclude",
]


class SkyMask(object):
    """
    This is a generic masking class that carries out selections of targets based on their RA,DEC

    Parameters
    ==========

    filename     - filename (no path) of the mask description
                 - will look only inside the target_selection/python/target_selection/masks/ directory

    masktype     - variety of mask description
                 - must be in ["mangle", "hpxmoc", "circles", ]

    sense        - "include" or "exclude" objects inside this mask?

    name         - Nickname for information only

    col_lon     - used by "circles" masks - gives the column name containing the longitude information
    col_lat      - used by "circles" masks - gives the column name containing the latitude information
    radius       - used by "circles" masks (assumed to be in degrees)
                   - if a string then gives the column name containing the radius information
                   - if a float then assume a uniform radius for all circles

    """

    def __init__(
        self,
        filename=None,
        masktype="mangle",
        sense="include",
        name=None,
        col_lon="racen",
        col_lat="deccen",
        radius="radius",
    ):
        assert filename is not None
        assert len(filename) > 0
        self.filename = os.path.expanduser(filename)

        self.masktype = masktype.lower()
        assert self.masktype in known_masktypes

        self.sense = sense.lower()
        assert self.sense in known_senses

        self.col_lon = col_lon
        self.col_lat = col_lat
        self.radius = radius

        if name is not None:
            self.name = name
        else:
            self.name = self.filename

    def __str__(self):
        return f"SkyMask({self.masktype}) {self.name}"

    def apply_mangle(self, lon=None, lat=None):
        try:
            m = mangle.Mangle(self.filename)
        except:
            raise Exception(f"Unable to open mangle mask file: {self.filename}")

        m_mask = m.contains(lon, lat)

        return m_mask

    def apply_hpxmoc(self, lon=None, lat=None):
        try:
            moc = MOC()
            moc.read(self.filename, filetype="fits")
        except:
            raise Exception(f"Unable to find/open Healpix MOC file: {self.filename}")

        # get the moc nside at the max resolution of the MOC
        nside = hp.order2nside(moc.order)
        # get the healpix pixels indices of the targets at the max resolution
        idx = hp.ang2pix(nside, lon, lat, lonlat=True, nest=True)

        m_mask = moc.contains(idx)
        return m_mask

    def apply_circles(self, lon=None, lat=None):
        # read the list of circles from the file

        try:
            hdul = fits.open(self.filename)
        except:
            raise Exception(f"Failed to find/open circles file: {self.filename}")

        data = hdul[1].data

        try:
            circle_lon = data[self.col_lon]
            circle_lat = data[self.col_lat]
        except:
            raise Exception(
                f"Failed to find circle columns: {self.filename}[1][{self.col_lon},{self.col_lat}]"
            )

        if isinstance(self.radius, float):
            circle_radius = np.full(len(data), self.radius)
        elif isinstance(self.radius, str):
            try:
                circle_radius = data[self.radius]
            except:
                raise Exception(f"Failed to find circle column: {self.filename}[1][{self.radius}]")
        else:
            raise Exception(f"Cannot interpret radius: {self.radius}")

        hdul.close()

        ### need to use a single radius for all objects
        seplimit = u.Quantity(circle_radius.max(), unit="deg")

        coords_c = SkyCoord(circle_lon, circle_lat, frame="icrs", unit="deg")
        coords_t = SkyCoord(lon, lat, frame="icrs", unit="deg")
        idx_t, idx_c, d2d, d3d = search_around_sky(coords_t, coords_c, seplimit)

        ## filter the list of matches to use the per-circle radius
        i_t_ok = [i_t for i_t, i_c, d in zip(idx_t, idx_c, d2d) if d.deg < circle_radius[i_c]]

        m_mask = np.zeros(len(lon), np.bool)
        m_mask[i_t_ok] = True

        return m_mask

    def apply(self, lon=None, lat=None, flags=None):
        nlon = len(lon) if hasattr(lon, "__len__") else 1
        nlat = len(lat) if hasattr(lat, "__len__") else 1

        assert nlon == nlat
        if nlon == 0:
            return None

        if flags is not None:
            nflags = len(flags) if hasattr(flags, "__len__") else 1
            assert nflags == nlat

        if self.masktype == "mangle":
            m_mask = self.apply_mangle(lon, lat)
        elif self.masktype == "hpxmoc":
            m_mask = self.apply_hpxmoc(lon, lat)
        elif self.masktype == "circles":
            m_mask = self.apply_circles(lon, lat)
        else:
            raise Exception(f"Unknown masktype: {self.masktype}")

        # if necessary, adjust the sense of the sky mask
        if self.sense == "exclude":
            m_mask = ~m_mask

        # if input flags were supplied then do a logical AND with them
        if flags is None:
            return m_mask
        else:
            return m_mask & flags
