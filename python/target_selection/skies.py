#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-06
# @Filename: skies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import multiprocessing
import os
import warnings
from functools import partial

import healpy
import numpy
import pandas
from astropy.coordinates import SkyCoord, match_coordinates_sky
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from target_selection import log, manager
from target_selection.exceptions import TargetSelectionUserWarning


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

        prograded = numpy.zeros((len(pixels), int(2**(2 * (k_out - k_in)))),
                                dtype=int)
        prograded += pixels[numpy.newaxis].T
        for ii in range(int(k_out - k_in))[::-1]:

            prograded = numpy.left_shift(prograded, 2)

            repeats = numpy.repeat([0, 1, 2, 3], (4 ** ii) or 1)
            tiles = numpy.tile(repeats, prograded.shape[1] // len(repeats))

            prograded += tiles

        return prograded


def _process_tile(inputs, candidate_nside=None, tile_nside=None,
                  min_separation=None, mag_column=None, is_flux=False,
                  mag_threshold=None, downsample=False,
                  downsample_nside=None, seed=None):
    """Processes a tile from catalogue data."""

    tile, targets = inputs

    if len(targets) == 0:
        return False

    # For each target we calculate the corresponding pixel in candidate_nside
    # resolution.
    cpcol = f'pix_{candidate_nside}'
    targets[cpcol] = healpy.ang2pix(candidate_nside,
                                    targets.ra, targets.dec,
                                    nest=True, lonlat=True)

    # Get all the pixels in candidate_nside resolution that correspond to this
    # tile. Create a mask with all of them that are not occupied by a target.
    all_pix = nested_regrade(tile, tile_nside, candidate_nside)[0]
    mask_not_target = ~numpy.in1d(all_pix, targets[cpcol])

    # Create a DF with the candidate pixels and calculate their RAs and Decs.
    candidates = pandas.DataFrame(all_pix[mask_not_target], columns=[cpcol])
    candidates['ra'], candidates['dec'] = healpy.pix2ang(candidate_nside,
                                                         candidates[cpcol],
                                                         nest=True,
                                                         lonlat=True)

    # If we are using the magnitudes to correct the minimum separation, we
    # first convert the column to magnitudes if it's flux. Then select all the
    # targets with magnitude < mag_threshold and for each one mask out all
    # the candidate pixels that lie within the corrected separation.
    # TODO: We use a rectangular search here because doing a
    # SkyCoord.separation for each target would be very costly but there may
    # be a better way.
    if mag_column:

        if is_flux:
            mask_flux = targets[mag_column] > 0
            targets.loc[mask_flux, mag_column] = (
                22.5 - 2.5 * numpy.log10(targets[mask_flux][mag_column]))
            targets.loc[targets[mag_column] <= 0, mag_column] = numpy.nan

        btargets = targets[targets[mag_column] < mag_threshold]
        mask = numpy.ones(len(candidates), dtype=numpy.bool)

        for _, btarget in btargets.iterrows():
            sep_corr = min_separation + numpy.sqrt((mag_threshold -
                                                    btarget[mag_column]) / 0.2)
            sep_corr /= 3600.
            sep_corr_ra = sep_corr / numpy.cos(numpy.radians(btarget.dec))
            masked_out = ((candidates.ra > (btarget.ra - sep_corr_ra)) &
                          (candidates.ra < (btarget.ra + sep_corr_ra)) &
                          (candidates.dec > (btarget.dec - sep_corr)) &
                          (candidates.dec < (btarget.dec + sep_corr)))
            mask[masked_out] = False

        candidates = candidates.loc[mask, :].copy()

    # Determine the closest neighbout for each candidate pixel.
    c1 = SkyCoord(ra=candidates.ra.to_numpy(),
                  dec=candidates.dec.to_numpy(), unit='deg')
    c2 = SkyCoord(ra=targets.ra.to_numpy(),
                  dec=targets.dec.to_numpy(), unit='deg')

    target_idx, sep, _ = match_coordinates_sky(c1, c2, nthneighbor=1)

    matched_target = targets.iloc[target_idx, :]
    sep_arcsec = sep.value * 3600.

    # Add column with the separation.
    candidates['sep_neighbour'] = sep_arcsec

    if mag_column:
        candidates['mag_neighbour'] = matched_target.loc[:, mag_column].to_numpy()

    valid_skies = candidates.loc[candidates.sep_neighbour > min_separation, :].copy()
    valid_skies[f'tile_{tile_nside}'] = tile

    # If we are not downsampling, set the index and return.
    if downsample is False or downsample is None:
        if len(valid_skies) > 0:
            valid_skies.set_index(cpcol, inplace=True)
        else:
            return False
        return valid_skies

    # If downsampling, regrade each candidate pixel value to the nside we'll
    # use for downsampling.
    valid_skies['down_pix'] = nested_regrade(valid_skies[cpcol],
                                             candidate_nside,
                                             downsample_nside)

    if isinstance(downsample, int):

        assert downsample > 1, 'downsample must be > 1'

        # Calculate how many skies per downsample pixel we should get.
        # Each next order of a pixel has four nested pixels.
        k_tile = int(numpy.log2(tile_nside))
        k_downsample = int(numpy.log2(downsample_nside))
        n_downsample_pix = 4**(k_downsample - k_tile)

        n_skies_per_pix = (downsample // n_downsample_pix) + 1

        # Get the sky positions per downsampled tile.
        valid_skies_downsampled = (valid_skies.groupby('down_pix', axis=0)
                                   .apply(lambda x: x.sample(n=(n_skies_per_pix
                                                                if len(x) > n_skies_per_pix
                                                                else len(x)),
                                                             random_state=seed))
                                   .reset_index(drop=True))

    elif isinstance(downsample, (list, tuple, numpy.ndarray)):

        valid_skies_downsampled = valid_skies.loc[valid_skies[cpcol].isin(downsample), :]

    else:

        raise ValueError('invalid format for downsample.')

    if len(valid_skies_downsampled) > 0:
        valid_skies_downsampled.set_index(cpcol, inplace=True)
    else:
        return False

    return valid_skies_downsampled


def get_sky_table(database, table, output, tiles=None, append=False,
                  tile_nside=32, candidate_nside=32768, min_separation=10,
                  chunk_tiles=5, ra_column='ra', dec_column='dec',
                  mag_column=None, is_flux=False, mag_threshold=None,
                  downsample=False, downsample_nside=256, n_cpus=1, seed=None):
    r"""Identifies skies from a database table.

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
      valid skies. Alternatively, if ``mag_column`` and ``mag_threshold`` are
      defined, the separations to neighbours brighter than the threshold
      magnitude is corrected using the expression
      :math:`s^* = s + \sqrt{\dfrac{m_{thr}-m}{a}}` where
      :math:`s` is the minimum separation, :math:`m_{thr}` is the magnitude
      threshold, :math:`a=0.2` is a factor that takes into account the average
      seeing for APO and LCO.

    - If ``downsample`` is an integer, only that number of skies are returned
      for each tile. If so, the tile is divided in pixels of nside
      ``downsample_nside`` (which should be smaller than ``candidate_nside``
      but larger than ``tile_nside``) and for each downsample pixel
      ``int(downsample / downsample_npix) + 1`` skies are selected (or all the
      skies, if fewer are available), where ``downsample_npix`` is the number
      of pixels corresponding to the ``downsample_nside`` nside.

    - The process can be parallelised for each tile and the results are
      compiled in a single Pandas data frame that is saved to an HDF5 file.

    All pixel values use the nested ordering.

    Parameters
    ----------
    database : ~sdssdb.connection.PeeweeDatabaseConnection
        A valid database connection.
    table : str
        Name of the table to query, optionally schema-qualified.
    output : str
        Path to the HDF5 file to write.
    tiles : list
        A list of HEALPix pixels of nside ``tile_nside`` for which the sky
        selection will be done. If `None`, runs for all the pixels of nside
        ``tile_nside``.
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
    chunk_tiles : int
        Request the targets from this number of tiles at once, and then process
        them one by one. Reduces the number of queries at the expense of making
        each one return more rows and use more memory. This option is ignored
        if ``tiles`` is set.
    ra_column : str
        The name of the column in ``table`` that contains the Right Ascension
        coordinates, in degrees.
    dec_column : str
        The name of the column in ``table`` that contains the Declination
        coordinates, in degrees.
    mag_column : str
        The name of the column in ``table`` with the magnitude to be used to
        scale ``min_separation``.
    is_flux : bool
        If `True`, assumes the ``mag_column`` values are given in nanomaggies.
    mag_threshold : float
        The value below which the separation to neighbouring sources will be
        scaled.
    downsample : bool, int, or list
        The total number of skies to retrieve for each tile. If `False`,
        returns all candidate skies. ``downsample`` can also be a list of
        pixels with nside ``candidate_nside`` in which case only pixels in
        that list will be returned.
    downsample_nside : int
        The HEALPix nside used for downsampling. If ``downsample=True``, the
        resulting valid skies will be grouped by HEALPix pixels of this
        resolution. For each pixel a random sample will be drawn so that the
        total number of skies selected matches ``downsample``.
    n_cpus : int
        Number of CPUs to use for multiprocessing.
    seed : int
        The random state seed.

    """

    assert database.connected, 'database is not connected.'

    columns = (f'healpix_ang2ipix_nest('
               f'{tile_nside}, {ra_column}, {dec_column}) AS tile_{tile_nside}, '
               f'{ra_column}, {dec_column}')
    if mag_column and mag_threshold:
        columns += f', {mag_column}'

    if chunk_tiles is None:
        chunk_tiles = 1

    if tiles is None:
        tiles = numpy.arange(healpy.nside2npix(tile_nside))

    query = (f'SELECT {columns} FROM {table} '
             f'WHERE healpix_ang2ipix_nest('
             f'{tile_nside}, {ra_column}, {dec_column}) IN ({{values}}) '
             f'ORDER BY tile_{tile_nside};')

    pbar = manager.counter(total=len(tiles), desc='Tiles', unit='tiles')
    pbar.refresh()

    process_tile = partial(_process_tile,
                           candidate_nside=candidate_nside,
                           tile_nside=tile_nside,
                           min_separation=min_separation,
                           mag_column=mag_column, is_flux=is_flux,
                           mag_threshold=mag_threshold,
                           downsample=downsample,
                           downsample_nside=downsample_nside,
                           seed=seed)

    all_skies = None

    for tile_split in numpy.array_split(tiles, len(tiles) / chunk_tiles):

        values = ','.join(map(str, tile_split))

        targets = pandas.read_sql(query.format(values=values), database)
        targets.rename(columns={ra_column: 'ra', dec_column: 'dec'}, inplace=True)

        if n_cpus > 1:

            groups = targets.groupby(f'tile_{tile_nside}')

            # Increment the counter for each tile that didn't get any target.
            n_groups = len(groups)
            non_hit = len(tile_split) - n_groups
            pbar.update(non_hit)

            with multiprocessing.Pool(n_cpus) as pool:

                for valid_skies in pool.imap_unordered(process_tile, groups):
                    if valid_skies is not False and len(valid_skies) > 0:
                        if all_skies is None:
                            all_skies = valid_skies
                        else:
                            all_skies = all_skies.append(valid_skies)

                    pbar.update()

        else:

            for tile in tile_split:
                valid_skies = process_tile((tile, targets))
                if valid_skies is not False:
                    if all_skies is None and len(valid_skies) > 0:
                        all_skies = valid_skies
                    else:
                        all_skies = all_skies.append(valid_skies)

                pbar.update()

    key = 'data'
    hdf = pandas.HDFStore(output, mode='w')

    if append is False and key in hdf:
        hdf.remove(key)

    hdf.append(key, all_skies)

    hdf.close()


def plot_sky_density(file_or_data, nside, nside_plot=32, **kwargs):
    """Plots the number of skies as a HEALPix map.

    Parameters
    ----------
    file_or_data : str or ~pandas.DataFrame
        A HDF5 file with the sky catalogue or a Pandas data frame with the
        data. It is assumed that the file or data frame contains at least a
        column ``pix_X``, where ``X`` is the order corresponding to ``nside``.
    nside : int
        The HEALPix nside for the skies.
    nside_plot : int
        The HEALPix nside in which the data will be plotted.
    kwargs : dict
        Other keyword parameters to pass to healpy's
        `~healpy.visufunc.mollview`.

    Returns
    -------
    figure : `~matplotlib.figure.Figure`
        The Matplotlib figure with the plot.

    """

    if isinstance(file_or_data, str):
        data = pandas.read_hdf(file_or_data)
    else:
        data = file_or_data

    k = int(numpy.log2(nside))
    k_plot = int(numpy.log2(nside_plot))

    pix_plot = f'pix_{k_plot}'
    pix = f'pix_{k}'

    data[pix_plot] = nested_regrade(data[pix], nside, nside_plot)
    count = data.groupby(pix_plot).count()

    hmap = numpy.arange(healpy.nside2npix(nside_plot), dtype=numpy.float32)
    hmap[:] = healpy.UNSEEN
    hmap[count.index] = count[pix]

    figure = plt.figure()
    healpy.mollview(hmap, fig=figure.number, nest=True, **kwargs)

    return figure


def plot_skies(file_or_data, ra, dec, radius=1.5, targets=None,
               show_sky_buffer=False, buffer_radius=10.):
    """Plots the skies (and optionally targets) in a regions.

    Parameters
    ----------
    file_or_data : str or ~pandas.DataFrame
        A HDF5 file with the sky catalogue or a Pandas data frame with the
        data, including columns labelled ``ra`` and ``dec``.
    ra : float
        The right ascension of the centre of the field.
    dec : float
        The declination of the centre of the field.
    radius : float
        The FOV radius, in degrees.
    targets : list
        A target list or array as ``(ra, dec)``.
    show_sky_buffer : bool
        Plots the buffer regions around each sky.
    buffer_radius : float
        The sky buffer radius, in arcsec.

    Returns
    -------
    figure : `~matplotlib.figure.Figure`
        The Matplotlib figure with the plot.

    """

    if isinstance(file_or_data, str):
        data = pandas.read_hdf(file_or_data)
    else:
        data = file_or_data

    fig, ax = plt.subplots()

    cos_factor = numpy.cos(numpy.radians(dec))

    FOV = Ellipse((ra, dec), radius * 2 / cos_factor,
                  radius * 2, linewidth=2, color='None',
                  ec='k')
    ax.add_patch(FOV)

    # Plot skies
    centre = SkyCoord(ra=ra, dec=dec, unit='deg')
    sky_coords = SkyCoord(ra=data.ra.to_numpy(),
                          dec=data.dec.to_numpy(),
                          unit='deg')
    skies = data[sky_coords.separation(centre).value < radius]

    ax.scatter(skies.ra, skies.dec,
               marker='o', s=0.2, color='b', zorder=10, label='Skies')

    if show_sky_buffer:
        for _, sky in skies.iterrows():
            buffer = Ellipse((sky.ra, sky.dec),
                             buffer_radius / 3600. * 2,
                             buffer_radius / 3600. * 2,
                             color='y', ec='None',
                             alpha=0.2, zorder=0)
            ax.add_patch(buffer)

    # Plot targets
    if targets is not None:
        targets = numpy.array(targets)
        target_coords = SkyCoord(ra=targets[:, 0],
                                 dec=targets[:, 1],
                                 unit='deg')
        targets = targets[target_coords.separation(centre).value < radius]

        ax.scatter(targets[:, 0], targets[:, 1], marker='x', s=0.2,
                   color='r', zorder=20, label='Targets')

    ax.legend(loc='upper right')

    ax.set_xlim(ra - (radius - 0.1) / cos_factor,
                ra + (radius + 0.1) / cos_factor)
    ax.set_ylim(dec - radius - 0.1, dec + radius + 0.1)

    ax.set_xlabel('Right Ascension [deg]')
    ax.set_ylabel('Declination [deg]')

    ax.set_aspect(1 / cos_factor)

    fig.tight_layout()

    return fig


def create_sky_catalogue(database, tiles=None, **kwargs):
    """A script to generate a combined sky catalogue from multiple sources."""

    if not os.path.exists('gaia_skies.h5'):
        log.info('Procesing gaia_dr2_source.')
        get_sky_table(database, 'catalogdb.gaia_dr2_source', 'gaia_skies.h5',
                      mag_column='phot_g_mean_mag', mag_threshold=12,
                      downsample=2000, tiles=tiles, **kwargs)
    else:
        warnings.warn('Found file gaia_skies.h5', TargetSelectionUserWarning)

    # We use Gaia as the source for the downsampled candidates.
    gaia = pandas.read_hdf('gaia_skies.h5')
    downsample_list = gaia.index.drop_duplicates().values

    if not os.path.exists('ls8_skies.h5'):
        log.info('Procesing legacy_survey_dr8.')
        get_sky_table(database, 'catalogdb.legacy_survey_dr8', 'ls8_skies.h5',
                      mag_column='flux_g', is_flux=True, mag_threshold=12,
                      downsample=downsample_list, tiles=tiles, **kwargs)
    else:
        warnings.warn('Found file ls8_skies.h5', TargetSelectionUserWarning)

    if not os.path.exists('tmass_skies.h5'):
        log.info('Procesing twomass_psc.')
        get_sky_table(database, 'catalogdb.twomass_psc', 'tmass_skies.h5',
                      dec_column='decl', mag_column='h_m', mag_threshold=12,
                      downsample=downsample_list, tiles=tiles, **kwargs)
    else:
        warnings.warn('Found file tmass_skies.h5', TargetSelectionUserWarning)

    if not os.path.exists('tycho2_skies.h5'):
        log.info('Procesing tycho2.')
        get_sky_table(database, 'catalogdb.tycho2', 'tycho2_skies.h5',
                      ra_column='ramdeg', dec_column='demdeg',
                      mag_column='vtmag', mag_threshold=12,
                      downsample=downsample_list, tiles=tiles, **kwargs)
    else:
        warnings.warn('Found file tycho2_skies.h5', TargetSelectionUserWarning)

    if not os.path.exists('tmass_xsc_skies.h5'):
        log.info('Procesing twomass_xsc.')
        get_sky_table(database, 'catalogdb.twomass_xsc', 'tmass_xsc_skies.h5',
                      dec_column='decl', mag_column='h_m_k20fe', mag_threshold=14,
                      downsample=downsample_list, tiles=tiles, **kwargs)
    else:
        warnings.warn('Found file tmass_xsc_skies.h5', TargetSelectionUserWarning)

    skies = None
    col_order = []

    for file_ in ['gaia_skies.h5', 'ls8_skies.h5', 'tmass_skies.h5',
                  'tycho2_skies.h5', 'tmass_xsc_skies.h5']:

        table_name = file_[0:-9]

        table = pandas.read_hdf(file_).drop_duplicates()
        table.rename(columns={'sep_neighbour': f'sep_neighbour_{table_name}',
                              'mag_neighbour': f'mag_neighbour_{table_name}'},
                     inplace=True)

        table.loc[:, f'{table_name}_sky'] = True

        if skies is None:
            skies = table
        else:
            skies = skies.combine_first(table)

        skies.fillna({f'{table_name}_sky': False}, inplace=True)

        col_order += [f'{table_name}_sky',
                      f'sep_neighbour_{table_name}',
                      f'mag_neighbour_{table_name}']

    skies = skies.loc[:, ['ra', 'dec', 'down_pix', 'tile_32'] + col_order]
    skies.to_hdf('skies.h5', 'data')
