#!/usr/bin/env python
# encoding: utf-8
#
# @Author: Tom Dwelly
# @Date: Oct-2019
# @Filename: mag_flux.py
# @License: BSD 3-Clause
# @Copyright: Tom Dwelly


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
from math import log10

__invalid_mag_val = 99.9999
__invalid_mag_thresh = 90.0


def nMgy2AB(flux):
    ''' see http://www.sdss3.org/dr8/algorithms/magnitudes.php '''
    if np.ndim(flux) > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.nan_to_num(22.5 - 2.5*np.log10(flux), copy=False,
                                 nan=__invalid_mag_val,
                                 posinf=__invalid_mag_val,
                                 neginf=__invalid_mag_val)
    else:
        if flux <= 0.0:
            return __invalid_mag_val
        else:
            return 22.5 - 2.5*log10(flux)


def AB2nMgy(mag):
    ''' see http://www.sdss3.org/dr8/algorithms/magnitudes.php '''
    if np.ndim(mag) > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            flux = np.where(mag >= __invalid_mag_thresh, 0.0, 10.0**(-0.4*(mag-22.5)) )
    else:
        if mag >= __invalid_mag_thresh:
            flux = 0.0
        else:
            flux = 10.0**(-0.4*(mag-22.5))
    return flux
