#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-04-11
# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import contextlib
import io
import time

import pandas
from peewee import SQL, fn


class Timer:
    """Convenience context manager to time events.

    Modified from https://bit.ly/3ebdp3y.

    """

    def __enter__(self):
        self.start = time.time()
        self.end = None
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start

    @property
    def elapsed(self):

        if self.end:
            return self.interval

        return time.time() - self.start


def sql_apply_pm(ra_field, dec_field, pmra_field, pmdec_field,
                 epoch_delta, is_pmra_cos=True):
    """Constructs a SQL expression for applying proper motions to RA/Dec.

    Parameters
    ----------
    ra_field
        The Peewee field or value representing Right Ascension.
    dec_field
        The Peewee field or value representing Declination.
    pmra_field
        The Peewee field or value representing the proper motion in RA.
    pmdec_field
        The Peewee field or value representing the proper motion in Dec.
    epoch_delta
        A value or expression with the delta epoch, in years.
    is_pmra_cos : bool
        Whether the ``pmra_field`` includes the correction for the cosine of
        the declination

    Returns
    -------
    ra_dec : `~peewee:Expression`
        A tuple with the `~peewee:Expression` instances for the proper motion
        corrected RA and Dec.

    """

    pmra_field = fn.coalesce(pmra_field, 0) / 1000 / 3600
    pmdec_field = fn.coalesce(pmdec_field, 0) / 1000 / 3600

    if is_pmra_cos:
        cos_dec = fn.cos(fn.radians(fn.coalesce(dec_field, 0)))
        ra_corr = ra_field + pmra_field / cos_dec * fn.coalesce(epoch_delta, 0)
    else:
        ra_corr = ra_field + pmra_field * fn.coalesce(epoch_delta, 0)

    dec_corr = dec_field + pmdec_field * fn.coalesce(epoch_delta, 0)

    return ra_corr, dec_corr


def sql_iauname(ra_field, dec_field, prefix='SDSS J'):
    """Constructs a SQL expression for the IAU name from RA/Dec.

    Parameters
    ----------
    ra_field
        The Peewee field or value representing Right Ascension.
    dec_field
        The Peewee field or value representing Declination.
    prefix : str
        The prefix to add before the calculated IAU name.

    Returns
    -------
    iauname : `~peewee:Expression`
        The `~peewee:Expression` instance to construct the IAU name.

    """

    def dd2dms(dd, round=2, sign=False):

        degs = fn.trunc(dd).cast('INT')
        mins = fn.trunc((fn.abs(dd) - fn.abs(degs)) * 60).cast('INT')
        secs = (fn.round(((((fn.abs(dd) - fn.abs(degs)) * 60) - mins) * 60)
                         .cast('NUMERIC'), round).cast('REAL'))

        degs_fmt = '00' if sign is False else 'SG00'
        secs_fmt = '00.' + '0' * round

        return fn.trim(fn.to_char(degs, degs_fmt)).concat(
            fn.trim(fn.to_char(mins, '00'))).concat(
                fn.trim(fn.to_char(secs, secs_fmt)))

    return SQL('\'' + prefix + '\'').concat(
        dd2dms(ra_field / 15)).concat(
            dd2dms(dec_field, round=1, sign=True))


def copy_pandas(df, database, table_name, schema=None, columns=None):
    """Inserts a Pandas data frame using COPY.

    Parameters
    ----------
    df : ~pandas.DataFrame
        The Pandas data frame to copy or the inputs to be used to create one.
    database
        A database connection.
    table_name : str
        The name of the table into which to copy the data.
    schema : str
        The schema in which the table lives.
    columns : list
        The list of column names to copy.

    """

    if not isinstance(df, pandas.DataFrame):
        df = pandas.DataFrame(data=df)

    stream = io.StringIO()
    df.to_csv(stream, index=False, header=False)
    stream.seek(0)

    cursor = database.cursor()

    full_table_name = schema + '.' + table_name if schema else table_name

    with database.atomic():
        cursor.copy_from(stream, full_table_name, columns=columns, sep=',')

    return df


def get_epoch(xmodel):
    """Returns the epoch for an `.XMatchModel` in Julian years."""

    xmatch = xmodel._meta.xmatch
    fields = xmodel._meta.fields

    # If epoch == 0, make it null. This helps with q3c functions.
    if xmatch.epoch:
        epoch = xmatch.epoch
    else:
        epoch = fn.nullif(fields[xmatch.epoch_column], 0)

    if xmatch.epoch_format == 'jd':
        epoch = 2000 + (epoch - 2451545.0) / 365.25

    return epoch


@contextlib.contextmanager
def set_config_parameter(database, parameter, new_value, reset=True, log=None):
    """Temporarily a database configuration parameter."""

    new_value = new_value.upper()

    try:

        orig_value = database.execute_sql(f'show {parameter}').fetchone()[0].upper()
        value_changed = orig_value != new_value

        if value_changed:
            database.execute_sql(f'SET {parameter} = {new_value};')
            if log:
                log.debug(f'{parameter} is {new_value}.')

        yield

    finally:

        if reset and value_changed:
            database.execute_sql(f'SET {parameter} = {orig_value};')
            if log:
                log.debug(f'{parameter} reset to {orig_value}.')
