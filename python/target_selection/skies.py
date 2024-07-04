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

import enlighten
import healpy
import numpy
import pandas
import peewee
from astropy.coordinates import SkyCoord, match_coordinates_sky
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from mocpy import MOC

from target_selection import log
from target_selection.exceptions import TargetSelectionError, TargetSelectionUserWarning


warnings.filterwarnings("ignore", ".*invalid value encountered in power.*")


_known_flux_zpts = {
    "nMgy": 22.5,
    "Jy": 8.9,
    "mJy": 16.4,
    "uJy": 23.9,
}


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

    assert nside_in != nside_out, "nside_in cannot be equal to nside_out."

    k_in = numpy.log2(nside_in)
    k_out = numpy.log2(nside_out)

    assert k_in.is_integer() and k_out.is_integer(), "nside_in or nside_out are not power of 2."

    pixels = numpy.atleast_1d(pixels).astype(int)
    assert pixels.ndim == 1, "dimension of input pixels is invalid."

    # Npix = 12 * nside**2
    is_all = numpy.all(pixels <= 12 * nside_in**2)
    assert is_all, "some pixel indices are greater than the maximum allowed for nside_in."

    if k_in > k_out:
        degraded = numpy.right_shift(pixels, 2 * int(k_in - k_out))

        return degraded[0] if len(degraded) == 1 else degraded

    else:
        prograded = numpy.zeros((len(pixels), int(2 ** (2 * (k_out - k_in)))), dtype=int)
        prograded += pixels[numpy.newaxis].T
        for ii in range(int(k_out - k_in))[::-1]:
            prograded = numpy.left_shift(prograded, 2)

            repeats = numpy.repeat([0, 1, 2, 3], (4**ii) or 1)
            tiles = numpy.tile(repeats, prograded.shape[1] // len(repeats))

            prograded += tiles

        return prograded


def downsample(
    df,
    nsample=2048,
    tile_column="tile_32",
    tile_nside=32,
    candidate_nside=32768,
    downsample_nside=256,
    downsample_data=None,
    seed=None,
):
    """Downsamples valid sky positions for a tile."""

    df["selected"] = False

    if "down_pix" not in df:
        df.loc[:, "down_pix"] = nested_regrade(df.index, candidate_nside, downsample_nside)

    if downsample_data is not None:
        cond = True
        if tile_column in df:
            tile = df.iloc[0][tile_column]
            cond = cond & (downsample_data[tile_column] == tile)
        if "valid" in df:
            cond = cond & downsample_data.valid
        down_tile = downsample_data.loc[cond]
        selected = df.reindex(down_tile.index).dropna()
        selected = selected.loc[selected.valid]
        df.loc[selected.index, "selected"] = True
        if len(selected) == nsample:
            return df

    # Step one: get as many valid sky positions as possible in random order.
    n_selected = df.selected.sum()
    n_valid_left = (df.valid & ~df.selected).sum()

    # Check if we are done.
    if n_selected >= nsample:
        return df
    elif (~df.selected).sum() == 0:
        return df

    # Calculate how many skies per downsample pixel to get.
    k_tile = int(numpy.log2(tile_nside))
    k_downsample = int(numpy.log2(downsample_nside))
    n_downsample_pix = 4 ** (k_downsample - k_tile)
    n_skies_per_pix = nsample // n_downsample_pix
    # reset index line here is new and it is meant to avoid problem in sample
    # with invalid ids while using weights
    # line with set_index here is not part of the original code and it is
    # meant to avoid problem in sample with invalid ids using weights
    if n_valid_left > 0:
        not_selected = df[(~df.selected) & df.valid]
        not_selected.reset_index(inplace=True)
        weights = not_selected.sep_neighbour / not_selected.sep_neighbour.sum()
        valid_selected = not_selected.groupby("down_pix").sample(
            n=n_skies_per_pix,
            replace=True,
            random_state=seed,
            weights=weights,
        )
        valid_selected.drop_duplicates(inplace=True)
        valid_selected.set_index("pix_32768", inplace=True)
        df.loc[valid_selected.index, "selected"] = True

    if df.selected.sum() >= nsample or (~df.selected).sum() == 0:
        return df

    # Step two: complete downsample pixels without enough skies with the
    # invalid skies that have the largest separation.

    assigned = df.groupby("down_pix").apply(lambda x: len(x[x.selected]))
    n_pix_missing = assigned[assigned < n_skies_per_pix]

    invalid_sorted = (
        df.loc[df.down_pix.isin(n_pix_missing.index) & ~df.selected]
        .groupby("down_pix")
        .apply(lambda x: x.sort_values(["valid", "sep_neighbour"], ascending=False))
    )

    if invalid_sorted.size == 0:
        return df

    if "down_pix" in invalid_sorted.index.names:
        invalid_sorted.reset_index("down_pix", drop=True, inplace=True)

    invalid_selected = invalid_sorted.groupby("down_pix").apply(
        lambda x: x.head(n_skies_per_pix - assigned.loc[x.iloc[0].down_pix])
    )

    if invalid_selected.size == 0:
        return df

    if "down_pix" in invalid_selected.index.names:
        invalid_selected.reset_index("down_pix", drop=True, inplace=True)

    df.loc[invalid_selected.index, "selected"] = True

    return df


def _process_tile(
    tile,
    database_params=None,
    candidate_nside=None,
    tile_nside=None,
    query=None,
    min_separation=None,
    is_flux=False,
    flux_unit="nMgy",
    mag_threshold=None,
    scale_a=0.2,
    scale_b=1,
    nsample=None,
    downsample_data=None,
    downsample_nside=None,
    calculate_min_separation=True,
    seed=None,
):
    """Processes a tile from catalogue data."""

    db = peewee.PostgresqlDatabase(**database_params)
    db.connect()

    targets = pandas.read_sql(query.format(tile=tile), db)

    db.close()

    has_mag = "mag" in targets
    has_radius = "radius" in targets

    if "mag" in targets and not mag_threshold:
        raise TargetSelectionError("mag_threshold required if mag_column is set.")

    if len(targets) == 0:
        return False

    # For each target we calculate the corresponding pixel in candidate_nside
    # resolution.
    cpcol = f"pix_{candidate_nside}"

    targets[cpcol] = healpy.ang2pix(
        candidate_nside, targets.ra.array, targets.dec.array, nest=True, lonlat=True
    )

    # Get all the pixels in candidate_nside resolution that correspond to this
    # tile. Create a mask with all of them that are not occupied by a target.
    all_pix = nested_regrade(tile, tile_nside, candidate_nside)[0]
    mask_not_target = ~numpy.isin(all_pix, targets[cpcol])

    # Create a DF with the candidate pixels and calculate their RAs and Decs.
    candidates = pandas.DataFrame(all_pix[mask_not_target], columns=[cpcol])
    candidates["ra"], candidates["dec"] = healpy.pix2ang(
        candidate_nside,
        candidates[cpcol],
        nest=True,
        lonlat=True,
    )

    candidates["valid"] = True

    # If we are using the magnitudes to correct the minimum separation, we
    # first convert the column to magnitudes if it's flux. Then for each pixel
    # we mask out all the candidate pixels that lie within the corrected
    # separation.

    if is_flux:
        mask_flux = targets.mag > 0

        try:
            zpt = _known_flux_zpts[flux_unit]
        except BaseException:
            raise Exception(f"Unknown flux_unit: {flux_unit}")

        targets.loc[mask_flux, "mag"] = zpt - 2.5 * numpy.log10(targets[mask_flux].mag)
        targets.loc[targets.mag <= 0, "mag"] = numpy.nan

    vectors = healpy.pixelfunc.ang2vec(targets.ra, targets.dec, lonlat=True)

    all_masked = []
    row = 0
    for _, target in targets.iterrows():
        if has_radius:
            sep_corr = max(min_separation, target.radius)
        elif has_mag and mag_threshold is not None:
            min_sep_corr = numpy.power(mag_threshold - target.mag, scale_b) / scale_a
            sep_corr = min_separation + min_sep_corr
        else:
            sep_corr = min_separation

        masked = healpy.query_disc(
            candidate_nside,
            vec=vectors[row],
            radius=numpy.deg2rad(sep_corr / 3600.0),
            inclusive=True,
            fact=4,
            nest=True,
        )

        all_masked += masked.tolist()

        row += 1

    candidates.loc[candidates[cpcol].isin(all_masked), "valid"] = False

    # Determine the closest neighbout for each candidate pixel.
    if calculate_min_separation:
        c1 = SkyCoord(ra=candidates.ra.to_numpy(), dec=candidates.dec.to_numpy(), unit="deg")
        c2 = SkyCoord(ra=targets.ra.to_numpy(), dec=targets.dec.to_numpy(), unit="deg")

        target_idx, sep, _ = match_coordinates_sky(c1, c2, nthneighbor=1)

        matched_target = targets.iloc[target_idx, :]
        sep_arcsec = sep.value * 3600.0

        # Add column with the separation.
        candidates["sep_neighbour"] = sep_arcsec

        if has_mag:
            candidates["mag_neighbour"] = matched_target.loc[:, "mag"].to_numpy()

    candidates.loc[:, f"tile_{tile_nside}"] = tile

    # If we are not downsampling, set the index and return.
    if downsample is False or downsample is None:
        if len(candidates) > 0:
            candidates.set_index(cpcol, inplace=True)
        else:
            return False
        return candidates

    # If downsampling, regrade each candidate pixel value to the nside we'll
    # use for downsampling.
    candidates.loc[:, "down_pix"] = nested_regrade(
        candidates[cpcol],
        candidate_nside,
        downsample_nside,
    )

    if len(candidates) > 0:
        candidates.set_index(cpcol, inplace=True)
    else:
        return False

    if downsample_data is None and nsample is None:
        return candidates
    else:
        if downsample_data is not None:
            assert nsample is not None, "downsample_data requires nsample to be defined."
        downsampled = downsample(
            candidates,
            tile_column=f"tile_{tile_nside}",
            tile_nside=tile_nside,
            nsample=nsample,
            candidate_nside=candidate_nside,
            downsample_nside=downsample_nside,
            downsample_data=downsample_data,
            seed=seed,
        )
        return downsampled[downsampled.selected]


def get_sky_table(
    database,
    table,
    output,
    tiles=None,
    tile_nside=32,
    candidate_nside=32768,
    min_separation=10,
    ra_column="ra",
    dec_column="dec",
    mag_column=None,
    is_flux=False,
    radius_column=None,
    flux_unit="nMgy",
    scale_a=0.2,
    scale_b=1.0,
    mag_threshold=None,
    calculate_min_separation=True,
    nsample=2048,
    downsample_data=None,
    downsample_nside=256,
    n_cpus=1,
    seed=None,
):
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
      candidates. For each target with magnitude less than ``mag_threshold``
      all pixels within a ``min_separation`` are rejected. The remaining pixels
      are considered valid skies. Alternatively, if ``mag_column`` and
      ``mag_threshold`` are defined, the minimum separation to valid pixels is
      corrected using the expression
      :math:`s^* = s + \dfrac{(m_{thr}-m)^{\beta}}{a}` where
      :math:`s` is the minimum separation, :math:`m_{thr}` is the magnitude
      threshold, :math:`a=0.2` and :math:`\beta=1.0` are factors that
      control the relationship between the star's magnitude and the exclusion
      radius.

    - If ``nsample``, only that number of skies are returned for each tile.
      The tile is divided in pixels of nside ``downsample_nside`` (which must
      be smaller than ``candidate_nside`` but larger than ``tile_nside``) and
      for each downsample pixel ``int(downsample / downsample_npix) + 1`` skies
      are selected. If not enough valid positions can be selected in this way,
      each downsample pixel quota is completed with the invalid positions that
      have a larger separation to their nearest neighbour.

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
    tile_nside : int
        The HEALPix nside to use to tile the all-sky catalogue.
    candidate_nside : int
        The HEALPix nside used to identify candidate pixels in each tile.
        Candidates are then checked to confirm that their closest neighbour
        is at least ``min_separation`` arcsec away.
    min_separation : int
        The minimum separation, in arcsec, between skies and their closest
        neighbour in the catalogue.
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
        If `True`, assumes the ``mag_column`` values are given as fluxes
        (units of flux_unit).
    flux_unit : str
        Gives the units of flux in the 'mag_column' - known values 'nMgy', 'Jy'
    mag_threshold : float
        The value below which the separation to neighbouring sources will be
        scaled.
    radius_column : str
         Name of the database column that provided the object radius
         (an alternative to mag_column useful for extended sources)
    scale_a : float
         Value of :math:`a` in the radius vs mag relationship
    scale_b : float
         Value of :math:`\beta` in the radius vs mag relationship
    calculate_min_separation : bool
        If `True`, calculates the separation to the nearest neighbour for
        each candidate sky position.
    nsample : int or None
        The total number of skies to retrieve for each tile. If `None`,
        returns all candidate skies.
    downsample_nside : int
        The HEALPix nside used for downsampling. If ``nside``, the
        resulting valid skies will be grouped by HEALPix pixels of this
        resolution. For each pixel a random sample will be drawn so that the
        total number of skies selected matches ``nsample``.
    downsample_data : pandas.DataFrame
        A data frame with previously selected skies that will be used to
        downsample the sky candidates. This is useful when trying to create a
        sky catalogue from multiple tables to ensure that the selected sky
        positions match across the various tables. If not enough valid skies
        can be selected from the sample in the data frame, they will be
        completed up to ``nsample``.
    n_cpus : int
        Number of CPUs to use for multiprocessing.
    seed : int
        The random state seed.

    Returns
    -------
    skies : pandas.DataFrame
        The list of selected skies.

    """

    assert database.connected, "database is not connected."

    columns = (
        f"healpix_ang2ipix_nest("
        f"{tile_nside}, {ra_column}, {dec_column}) "
        f"AS tile_{tile_nside}, "
        f"{ra_column} AS ra, {dec_column} AS dec"
    )

    if mag_column and mag_threshold:
        columns += f", {mag_column} AS mag"

    if radius_column:
        columns += f", {radius_column} as radius"

    if tiles is None:
        tiles = numpy.arange(healpy.nside2npix(tile_nside))

    if downsample_data is not None:
        downsample_data = downsample_data.loc[:, ["tile_32", "valid"]]
        downsample_data = downsample_data.loc[downsample_data.valid]

    query = (
        f"SELECT {columns} FROM {table} "
        f"WHERE healpix_ang2ipix_nest("
        f"{tile_nside}, {ra_column}, {dec_column}) = {{tile}};"
    )

    pbar = enlighten.Counter(total=len(tiles), desc=output, unit="tiles")
    pbar.refresh()

    database_params = database.connection_params
    database_params.update({"database": database.dbname})

    process_tile = partial(
        _process_tile,
        database_params=database_params,
        query=query,
        candidate_nside=candidate_nside,
        tile_nside=tile_nside,
        min_separation=min_separation,
        is_flux=is_flux,
        flux_unit=flux_unit,
        scale_a=scale_a,
        scale_b=scale_b,
        mag_threshold=mag_threshold,
        nsample=nsample,
        downsample_data=downsample_data,
        calculate_min_separation=calculate_min_separation,
        downsample_nside=downsample_nside,
        seed=seed,
    )

    all_skies = None

    with multiprocessing.Pool(n_cpus) as pool:
        for tile_skies in pool.imap_unordered(process_tile, tiles, chunksize=5):
            if tile_skies is not False and len(tile_skies) > 0:
                if all_skies is None:
                    all_skies = tile_skies
                else:
                    all_skies = all_skies.append(tile_skies)

            pbar.update()

    if all_skies is None:
        return all_skies

    # Not sure why this is needed but it seems downsampling sometimes
    # converts the "valid" column to object.
    all_skies.valid = all_skies.valid.astype(bool)

    # Downcast some columns
    all_skies = all_skies.astype({f"tile_{tile_nside}": numpy.int16})

    if "sep_neighbour" in all_skies:
        all_skies = all_skies.astype({"sep_neighbour": numpy.float32})
    if "mag_neighbour" in all_skies:
        all_skies = all_skies.astype({"mag_neighbour": numpy.float32})

    if downsample_nside is not None:
        all_skies = all_skies.astype({"down_pix": numpy.int32})

    all_skies.to_hdf(output, "data")

    return all_skies


def plot_sky_density(file_or_data, nside, pix_column=None, nside_plot=32, **kwargs):
    """Plots the number of skies as a HEALPix map.

    Parameters
    ----------
    file_or_data : str or ~pandas.DataFrame
        A HDF5 file with the sky catalogue or a Pandas data frame with the
        data. It is assumed that the file or data frame contains at least a
        column ``pix_X``, where ``X`` is the order corresponding to ``nside``.
    nside : int
        The HEALPix nside for the skies.
    pix_column : str or None
        The column that contains the HEALPix pixels. If `None`, the index of
        the data frame will be used.
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
        data = file_or_data.copy()

    if not pix_column:
        pix_column = data.index.name
        data.reset_index(inplace=True)

    data[f"pix_{nside_plot}"] = nested_regrade(data[pix_column], nside, nside_plot)
    count = data.groupby(f"pix_{nside_plot}").count()

    hmap = numpy.arange(healpy.nside2npix(nside_plot), dtype=numpy.float32)
    hmap[:] = healpy.UNSEEN
    hmap[count.index] = count[pix_column]

    figure = plt.figure()
    healpy.mollview(hmap, fig=figure.number, nest=True, **kwargs)

    return figure


def plot_skies(
    file_or_data, ra, dec, radius=1.5, targets=None, show_sky_buffer=False, buffer_radius=10.0
):
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

    FOV = Ellipse(
        (ra, dec),
        radius * 2 / cos_factor,
        radius * 2,
        linewidth=2,
        color="None",
        ec="k",
    )
    ax.add_patch(FOV)

    # Plot skies
    centre = SkyCoord(ra=ra, dec=dec, unit="deg")
    sky_coords = SkyCoord(ra=data.ra.to_numpy(), dec=data.dec.to_numpy(), unit="deg")
    skies = data[sky_coords.separation(centre).value < radius]

    ax.scatter(skies.ra, skies.dec, marker="o", s=0.2, color="b", zorder=10, label="Skies")

    if show_sky_buffer:
        for _, sky in skies.iterrows():
            buffer = Ellipse(
                (sky.ra, sky.dec),
                buffer_radius / 3600.0 * 2,
                buffer_radius / 3600.0 * 2,
                color="y",
                ec="None",
                alpha=0.2,
                zorder=0,
            )
            ax.add_patch(buffer)

    # Plot targets
    if targets is not None:
        targets = numpy.array(targets)
        target_coords = SkyCoord(ra=targets[:, 0], dec=targets[:, 1], unit="deg")
        targets = targets[target_coords.separation(centre).value < radius]

        ax.scatter(
            targets[:, 0],
            targets[:, 1],
            marker="x",
            s=0.2,
            color="r",
            zorder=20,
            label="Targets",
        )

    ax.legend(loc="upper right")

    ax.set_xlim(ra - (radius - 0.1) / cos_factor, ra + (radius + 0.1) / cos_factor)
    ax.set_ylim(dec - radius - 0.1, dec + radius + 0.1)

    ax.set_xlabel("Right Ascension [deg]")
    ax.set_ylabel("Declination [deg]")

    ax.set_aspect(1 / cos_factor)

    fig.tight_layout()

    return fig


# This is the main function to use with catalogdb as catalog
def create_sky_catalogue(database, tiles=None, **kwargs):
    """A script to generate a combined sky catalogue from multiple sources."""

    default_mag_threshold = 14.0
    default_min_separation = 5.0
    default_param_a = 0.15
    default_param_b = 1.5
    nsample = 2048

    if not os.path.exists("tmass_skies.h5"):
        log.info("Procesing twomass_psc.")
        get_sky_table(
            database,
            "catalogdb.twomass_psc",
            "tmass_skies.h5",
            dec_column="decl",
            mag_column="h_m",
            mag_threshold=default_mag_threshold,
            min_separation=default_min_separation,
            scale_a=default_param_a,
            scale_b=default_param_b,
            nsample=nsample,
            tiles=tiles,
            **kwargs,
        )
    else:
        warnings.warn("Found file tmass_skies.h5", TargetSelectionUserWarning)

    # We use 2MASS as the source for the downsampled candidates.
    tmass = pandas.read_hdf("tmass_skies.h5")

    if not os.path.exists("gaia_skies.h5"):
        log.info("Procesing gaia_dr2_source.")
        get_sky_table(
            database,
            "catalogdb.gaia_dr2_source",
            "gaia_skies.h5",
            mag_column="phot_g_mean_mag",
            mag_threshold=default_mag_threshold,
            min_separation=default_min_separation,
            scale_a=default_param_a,
            scale_b=default_param_b,
            nsample=nsample,
            downsample_data=tmass,
            tiles=tiles,
            **kwargs,
        )
    else:
        warnings.warn("Found file gaia_skies.h5", TargetSelectionUserWarning)

    if not os.path.exists("ps1dr2_skies.h5"):
        log.info("Procesing ps1dr2.")
        get_sky_table(
            database,
            "catalogdb.panstarrs1",
            "ps1dr2_skies.h5",
            mag_column="r_stk_psf_flux",
            is_flux=True,
            flux_unit="Jy",
            mag_threshold=default_mag_threshold,
            min_separation=default_min_separation,
            scale_a=default_param_a,
            scale_b=default_param_b,
            nsample=nsample,
            downsample_data=tmass,
            tiles=tiles,
            **kwargs,
        )
    else:
        warnings.warn("Found file ps1dr2_skies.h5", TargetSelectionUserWarning)

    if not os.path.exists("ls8_skies.h5"):
        log.info("Procesing legacy_survey_dr8.")
        get_sky_table(
            database,
            "catalogdb.legacy_survey_dr8",
            "ls8_skies.h5",
            mag_column="flux_r",
            is_flux=True,
            flux_unit="nMgy",
            mag_threshold=default_mag_threshold,
            min_separation=default_min_separation,
            scale_a=default_param_a,
            scale_b=default_param_b,
            nsample=nsample,
            downsample_data=tmass,
            tiles=tiles,
            **kwargs,
        )
    else:
        warnings.warn("Found file ls8_skies.h5", TargetSelectionUserWarning)

    if not os.path.exists("tycho2_skies.h5"):
        log.info("Procesing tycho2.")
        get_sky_table(
            database,
            "catalogdb.tycho2",
            "tycho2_skies.h5",
            ra_column="ramdeg",
            dec_column="demdeg",
            mag_column="vtmag",
            mag_threshold=default_mag_threshold,
            min_separation=default_min_separation,
            scale_a=default_param_a,
            scale_b=default_param_b,
            nsample=nsample,
            downsample_data=tmass,
            tiles=tiles,
            **kwargs,
        )
    else:
        warnings.warn("Found file tycho2_skies.h5", TargetSelectionUserWarning)

    if not os.path.exists("tmass_xsc_skies.h5"):
        log.info("Procesing twomass_xsc.")
        get_sky_table(
            database,
            "catalogdb.twomass_xsc",
            "tmass_xsc_skies.h5",
            dec_column="decl",
            radius_column="r_ext",
            # mag_column='h_m_k20fe', mag_threshold=14,
            nsample=nsample,
            downsample_data=tmass,
            tiles=tiles,
            **kwargs,
        )
    else:
        warnings.warn("Found file tmass_xsc_skies.h5", TargetSelectionUserWarning)

    skies = None
    col_order = []

    for file_ in [
        "gaia_skies.h5",
        "ls8_skies.h5",
        "ps1dr2_skies.h5",
        "tmass_skies.h5",
        "tycho2_skies.h5",
        "tmass_xsc_skies.h5",
    ]:
        table_name = file_[0:-9]

        if not os.path.exists(file_):
            warnings.warn(f"File {file_} does not exists.", TargetSelectionUserWarning)
            continue

        table = pandas.read_hdf(file_).drop_duplicates()
        table.rename(
            columns={
                "sep_neighbour": f"sep_neighbour_{table_name}",
                "mag_neighbour": f"mag_neighbour_{table_name}",
                "valid": f"valid_{table_name}",
                "selected": f"selected_{table_name}",
            },
            inplace=True,
        )

        if skies is None:
            skies = table
        else:
            skies = skies.combine_first(table)

        for col in [
            f"valid_{table_name}",
            f"selected_{table_name}",
            f"sep_neighbour_{table_name}",
            f"mag_neighbour_{table_name}",
        ]:
            if col in skies:
                col_order.append(col)

    skies = skies.loc[:, ["ra", "dec", "down_pix", "tile_32"] + col_order]
    for col in skies:
        if col.startswith("valid_") or col.startswith("selected_"):
            skies[col].fillna(False, inplace=True)

    # Now do some masking based on healpixels that are flagged by a
    # conservative veto_mask. This is intended to catch a few cases
    # where the nearest neighbour star is not the one that is
    # contributing the most flux at a candidate sky location

    # Now get the list of all healpixels that are within circles around
    # very bright stars:
    # veto_mask = create_veto_mask(database=database, nside=32768,
    #                              moc_filename='veto_from_tycho2.moc',
    #                              table='catalogdb.tycho2',
    #                              ra_column='ramdeg', dec_column='demdeg',
    #                              mag_column='vtmag', mag_threshold=12,
    #                              min_separation=15.0, param_a=0.15,
    #                              param_b=1.5, overwrite=True)

    # # Now flag any skies that have a 'pix_32768' that is in the veto_mask
    # skies.loc[:, 'tycho2_veto'] = False
    # skies.loc[skies.index.isin(veto_mask), 'tycho2_veto'] = True

    skies.tile_32 = skies.tile_32.astype(int)
    skies.down_pix = skies.down_pix.astype(int)
    skies.to_hdf("skies.h5", "data")


def create_veto_mask(
    database,
    nside=32768,
    moc_filename=None,
    overwrite=False,
    debug_limit=None,
    table="gaia_dr2_source",
    ra_column="ra",
    dec_column="dec",
    mag_column="phot_g_mean_mag",
    mag_threshold=12.0,
    min_separation=15.0,
    param_a=0.15,
    param_b=1.5,
):
    """Generate a list of healpixels that must be avoided because
    they lie within the near-zones of bright stars (and/or galaxies TBD).

    I have a hunch that generating this list once and then testing skies
    against it will be more efficient (and reliable) than checking candidate
    skies against a list of stars. This avoids potential pitfalls of nearest
    neighbour method (e.g. when second nearest neighbour is brighter than
    nearest neighbour).

    This doesn't take too long to run, so result does not need to be preserved
    long term. However, writing out a MOC file is a convenient way to visualise
    what has been done.

    Parameters
    ----------
    database : ~sdssdb.connection.PeeweeDatabaseConnection
        A valid database connection.
    nside : int
        HEALPix resolution of the returned healpix pixel list
    moc_filename : str
        Path to the MOC file to write (or None).
    overwrite: bool
        Whether to clobber the MOC file
    debug_limit: int
        Max number of stars to return in database query - debug purposes only
    ra_column : str
        The name of the column in ``table`` that contains the Right Ascension
        coordinates, in degrees.
    dec_column : str
        The name of the column in ``table`` that contains the Declination
        coordinates, in degrees.
    mag_column : str
        The name of the column in ``table`` with the magnitude to be used to
        scale ``min_separation``.
    param_a: float
        A parameter that controls the how the radius scales with magnitude
    param_b: float
        A parameter that controls the how the radius scales with magnitude

    Returns
    -------
    mask : ~numpy.ndarray
        A (numpy) array of healpixel indices (resolution ``nside``) that
        fall within the mask.

    """

    as2rad = numpy.pi / (180.0 * 3600.0)
    hpx_order = healpy.nside2order(nside)
    pixarea = healpy.nside2pixarea(nside)

    # get the list of gaia dr2 stars brighter than G=12 from the database

    query = (
        f"SELECT {ra_column},{dec_column},{mag_column} from "
        f"{table} WHERE {mag_column} < {mag_threshold} AND "
        f"{ra_column} IS NOT NULL AND {dec_column} IS NOT NULL"
    )
    if debug_limit is not None:
        query = query + f"limit {debug_limit} "

    targets = pandas.read_sql(query, database)
    print(
        f"Working on {len(targets):,} bright stars "
        f"({mag_column} < {mag_threshold}) from {table}"
    )

    # compute coords on unit sphere
    vector = healpy.pixelfunc.ang2vec(targets[ra_column], targets[dec_column], lonlat=True)

    # compute mag-dependent exclusion radii
    corr = numpy.power(mag_threshold - targets[mag_column], param_b) / param_a
    radius = as2rad * min_separation + corr

    ipix_list = []
    for v, r in zip(vector, radius):
        i = healpy.query_disc(nside, vec=v, radius=r, inclusive=True, fact=4, nest=True)
        if len(i) > 0:
            ipix_list.extend(list(i))

    # we only one copy of each masked pixel:
    ipix = numpy.unique(ipix_list)
    npix = len(ipix)
    print(f"Result: {npix:,} masked pixels (NSIDE={nside}), area={npix*pixarea:.4f} sqdeg")

    if moc_filename is not None:
        m = MOC.from_healpix_cells(ipix=ipix, depth=numpy.repeat(hpx_order, len(ipix)))
        m.write(moc_filename, format="fits", overwrite=overwrite)

    return ipix
