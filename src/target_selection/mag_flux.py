#!/usr/bin/env python
# encoding: utf-8
#
# @Author: Tom Dwelly
# @Date: Oct-2019
# @Filename: mag_flux.py
# @License: BSD 3-Clause
# @Copyright: Tom Dwelly

from __future__ import absolute_import, division, print_function, unicode_literals

from math import log10

import numpy as np


__invalid_mag_val = 99.9999
__invalid_mag_thresh = 90.0
""" see http://www.sdss3.org/dr8/algorithms/magnitudes.php """


# https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/pro/spec2d/readplugmap.pro
__psfflux_to_fiber2flux = {
    "u": (1.0 / 2.085),
    "g": (1.0 / 2.085),
    "r": (1.0 / 2.116),
    "i": (1.0 / 2.134),
    "z": (1.0 / 2.135),
}


def fluxmag(flux, zp=None):
    assert zp is not None
    if np.ndim(flux) > 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.nan_to_num(
                zp - 2.5 * np.log10(flux),
                copy=False,
                nan=__invalid_mag_val,
                posinf=__invalid_mag_val,
                neginf=__invalid_mag_val,
            )
    else:
        if flux <= 0.0:
            return __invalid_mag_val
        else:
            return zp - 2.5 * log10(flux)


def mag2flux(mag, zp=None):
    assert zp is not None
    if np.ndim(mag) > 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            flux = np.where(mag >= __invalid_mag_thresh, 0.0, 10.0 ** (-0.4 * (mag - zp)))
    else:
        if mag >= __invalid_mag_thresh:
            flux = 0.0
        else:
            flux = 10.0 ** (-0.4 * (mag - zp))
    return flux


def psfflux2fiber2flux(psfflux, filt=None):
    try:
        c = __psfflux_to_fiber2flux[filt]
    except:
        raise Exception(f"Unknown filter: {filt}")
    return psfflux * c


def fiber2flux2psfflux(fiber2flux, filt=None):
    try:
        c = __psfflux_to_fiber2flux[filt]
    except:
        raise Exception(f"Unknown filter: {filt}")
    return fiber2flux / c


def psfmag2fiber2mag(psfmag, filt=None):
    try:
        c = __psfflux_to_fiber2flux[filt]
    except:
        raise Exception(f"Unknown filter: {filt}")
    return psfmag - (2.5 * log10(c))


def fiber2mag2psfmag(fiber2mag, filt=None):
    try:
        c = __psfflux_to_fiber2flux[filt]
    except:
        raise Exception(f"Unknown filter: {filt}")
    return fiber2mag + (2.5 * log10(c))


def psfmag_minus_fiber2mag(filt=None):
    try:
        c = __psfflux_to_fiber2flux[filt]
    except:
        raise Exception(f"Unknown filter: {filt}")
    return 2.5 * log10(c)


def fiber2mag_minus_psfmag(filt=None):
    try:
        c = __psfflux_to_fiber2flux[filt]
    except:
        raise Exception(f"Unknown filter: {filt}")
    return -2.5 * log10(c)


# def nMgy2AB(flux):
#     return flux2mag(mag, 22.50)


# def Jy2AB(flux):
#     return flux2mag(mag, 8.90)


# def mJy2AB(flux):
#     return flux2mag(mag, 16.40)


# def uJy2AB(flux):
#     return flux2mag(mag, 23.90)


def AB2nMgy(mag):
    return mag2flux(mag, 22.50)


def AB2Jy(mag):
    return mag2flux(mag, 8.90)


def AB2mJy(mag):
    return mag2flux(mag, 16.40)


def AB2uJy(mag):
    return mag2flux(mag, 23.90)
