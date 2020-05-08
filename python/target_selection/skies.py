#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-06
# @Filename: skies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import multiprocessing
from functools import partial

import healpy
import numpy
import pandas
from astropy.coordinates import SkyCoord, match_coordinates_sky
from peewee import PostgresqlDatabase

from target_selection import manager


def nested_regrade(pixels, nside_in, nside_out):
    r"""Returns the parent/children pixels from a given HealPix nested pixel.

    The HealPix nested mode follows a quadrilateral tree pixel scheme (see
    Figure 1 in Górski et al. 2005) When the resolution (nside) increases,
    each pixel is divided in four pixels. The numbering of such pixels follows
    a simple binary representation achieved by appending binary digits to the
    new pixels. So, for example, pixel 22 (``b10110``) becomes four new pixels
    with binary numbers ``b1011000, b1011001, b1011010, and b1011011``,
    corresponding to decimal ``91, 92, 93, and 94``.

    This function returns the child pixels from a given pixel when going from
    resolution ``nside_in`` to ``nside_out``, if ``nside_out > nside_in``, or
    the parent pixel if ``nside_out < nside_in``.

    Note that this function works only for pixels using the nested indexing and
    should not be used with the ring indexing.

    Parameters
    ----------
    pixels : int or ~numpy.ndarray
        The pixels for which we want to get the parents/children. Can be a
        single integer or an array of indices.
    nside_in : int
        The ``nside`` of the input pixels. Must be one of :math:`2^k` where
        :math:`k\in[0,1,2,\dots]`.
    nside_out : int
        The destination ``nside``.

    Returns
    -------
    output : `int` or `~numpy.ndarray`
        If ``nside_out < nside_in`` and ``pixels`` is a single value, this
        will be an integer with the parent of the input pixel. If ``pixels``
        is an array the output will be an array of the same size in which each
        element is the parent of the corresponding pixel in the input array.
        If ``nside_out > nside_in`` and ``pixels`` is an integer, the output
        will be an array with the child pixels. The size of the array will be
        :math:`2^{2l}` where ``l=k_out-k_in`` (:math:`n_{side}=2^k`).
        If the input in an array of pixels, the output will be a 2D array in
        which each row contains the child pixels of the corresponding pixel
        in the input array.

    """

    assert nside_in != nside_out, 'nside_in cannot be equal to nside_out.'

    k_in = numpy.log2(nside_in)
    k_out = numpy.log2(nside_out)

    assert k_in.is_integer() and k_out.is_integer(), \
        'nside_in or nside_out are not power of 2.'

    pixels = numpy.atleast_1d(pixels).astype(int)
    assert pixels.ndim == 1, 'dimension of input pixels is invalid.'

    # Npix = 12 * nside**2
    assert numpy.all(pixels <= 12 * nside_in**2), \
        'some pixel indices are greater than the maximum allowed for nside_in.'

    if k_in > k_out:

        degraded = numpy.right_shift(pixels, 2 * int(k_in - k_out))

        return degraded[0] if len(degraded) == 1 else degraded

    else:

        prograded = numpy.zeros((len(pixels), int(2**(2 * (k_out - k_in)))), dtype=int)
        prograded += pixels[numpy.newaxis].T
        for ii in range(int(k_out - k_in))[::-1]:

            prograded = numpy.left_shift(prograded, 2)

            repeats = numpy.repeat([0, 1, 2, 3], (4 ** ii) or 1)
            tiles = numpy.tile(repeats, prograded.shape[1] // len(repeats))

            prograded += tiles

        return prograded


def _process_tile(pix, database_params=None, query=None,
                  candidate_nside=None, tile_nside=None,
                  min_separation=None, downsample=False,
                  downsample_nside=None):
    """Processes a tile from a catalogue in parallel."""

    k = int(numpy.log2(candidate_nside))

    # We cannot pass the database itself. We need to create the connection
    # inside the function that gets parallelised.
    dbname = database_params.pop('dbname')
    database = PostgresqlDatabase(dbname, **database_params)

    targets = pandas.read_sql(query.format(pix=pix), database)
    targets['pix'] = healpy.ang2pix(candidate_nside, targets.ra, targets.dec,
                                    nest=True, lonlat=True)

    all_pix = nested_regrade(pix, tile_nside, candidate_nside)[0]
    mask_not_target = ~numpy.in1d(all_pix, targets.pix)

    pix = f'pix_{k}'
    candidates = pandas.DataFrame(all_pix[mask_not_target], columns=[pix])
    candidates['ra'], candidates['dec'] = healpy.pix2ang(candidate_nside,
                                                         candidates[pix],
                                                         nest=True,
                                                         lonlat=True)

    c1 = SkyCoord(ra=candidates.ra.to_numpy(), dec=candidates.dec.to_numpy(), unit='deg')
    c2 = SkyCoord(ra=targets.ra.to_numpy(), dec=targets.dec.to_numpy(), unit='deg')
    sep_arcsec = match_coordinates_sky(c1, c2, nthneighbor=1)[1].value * 3600.

    candidates['sep_neighbour'] = sep_arcsec

    valid_skies = candidates[candidates['sep_neighbour'] > min_separation].loc[:]

    if not downsample:
        return valid_skies

    valid_skies['down_pix'] = nested_regrade(valid_skies.loc[:, pix],
                                             candidate_nside,
                                             downsample_nside)

    # Each next order of a pixel has four nested pixels.
    k_tile = int(numpy.log2(tile_nside))
    k_downsample = int(numpy.log2(downsample_nside))
    n_downsample_pix = 4**(k_downsample - k_tile)

    n_skies_per_pix = (downsample // n_downsample_pix) + 1

    valid_skies_downsampled = (valid_skies.groupby('down_pix', axis=0)
                               .apply(lambda x: x.sample(n=(n_skies_per_pix
                                                            if len(x) > n_skies_per_pix
                                                            else len(x))))
                               .reset_index(drop=True))

    # valid_skies_downsampled.drop(columns='down_pix', inplace=True)

    return valid_skies_downsampled


def create_sky_catalogue(database, table, output, append=False, tile_nside=32,
                         candidate_nside=32768, min_separation=10,
                         ra_column='ra', dec_column='dec', n_cpus=None,
                         downsample=False, downsample_nside=256):
    """Identifies skies in a database catalogue.

    Skies are selected using the following procedure:

    - The sky is divided in HEALPix "tiles" of nside ``tile_nside``. For each
      tile the catalogue is queried to retrieve all the targets that lie in
      the tile footprint. This assumes that the ``healpix_ang2ipix_nest``
      Postgresql function from `pg_healpix
      <https://github.com/segasai/pg_healpix>`__ is available. It's recommended
      that an index is created for the tiling norder to speed the query.

    - Each tile is subsequently divided in pixels of nside ``candidate_nside``.
      Pixels that contain a catalogue target are removed as possible sky
      candidates. For the remaining pixels a nearest-neighbour search is done
      using a KDTree algorithm. Those pixels with neighbours closer than
      ``min_separation`` are rejected. The remaining pixels are considered
      valid skies.

    - If ``downsample`` is an integer, only that number of skies are returned
      for each tile. If so, the tile is divided in pixels of nside
      ``downsample_nside`` (which should be smaller than ``candidate_nside``
      but larger than ``tile_nside``) and for each downsample pixel
      ``int(downsample / downsample_npix) + 1`` skies are selected (or all the
      skies, if fewer are available), where ``downsample_npix`` is the number
      of pixels corresponding to the ``downsample_nside`` nside.

    - The process is parallelised for each tile and the results are compiled
      in a single Pandas data frame that is saved to an HDF5 file.

    All pixel values use the nested ordering.

    Parameters
    ----------
    database : .PeeweeDatabaseConnection
        A valid database connection.
    table : str
        Name of the table to query, optionally schema-qualified.
    output : str
        Path to the HDF5 file to write.
    append : bool
        If `True`, appends to the output file if it exists. Otherwise
        overwrites it.
    tile_nside : int
        The HEALPix nside to use to tile the all-sky catalogue.
    candidate_nside : int
        The HEALPix nside used to identify candidate pixels in each tile.
        Candidates are then checked to confirm that their closest neighbour
        is at least ``min_separation`` arcsec away.
    min_separation : int
        The minimum separation, in arcsec, between skies and their closest
        neighbour in the catalogue.
    n_cpus : int
        Number of CPUs to use in parallel. If `None`, defaults to the number
        of cores.
    downsample : bool or int
        The total number of skies to retrieve for each tile. If `False`,
        returns all candidate skies.
    downsample_nside : int
        The HEALPix nside used for downsampling. If ``downsample=True``, the
        resulting valid skies will be grouped by HEALPix pixels of this
        resolution. For each pixel a random sample will be drawn so that the
        total number of skies selected matches ``downsample``.

    """

    n_cpus = n_cpus or multiprocessing.cpu_count()

    key = 'data'
    hdf = pandas.HDFStore(output, complevel=9)

    if append is False and key in hdf:
        hdf.remove(key)

    assert database.connected, 'database is not connected.'

    query = (f'SELECT {ra_column}, {dec_column} FROM {table} '
             f'WHERE healpix_ang2ipix_nest({tile_nside}, {ra_column}, {dec_column}) = {{pix}};')

    n_tiles = healpy.nside2npix(tile_nside)
    pbar = manager.counter(total=n_tiles, desc='Tiles', unit='tiles')

    process_tile = partial(_process_tile,
                           database_params=database.connection_params.copy(),
                           query=query,
                           candidate_nside=candidate_nside,
                           tile_nside=tile_nside,
                           min_separation=min_separation,
                           downsample=downsample,
                           downsample_nside=downsample_nside)

    with multiprocessing.Pool(processes=n_cpus) as pool:
        for valid_skies in pool.imap_unordered(process_tile, range(n_tiles)):
            hdf.append(key, valid_skies)
            pbar.update()

    hdf.close()
