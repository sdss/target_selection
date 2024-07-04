#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-04-11
# @Filename: utils.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import contextlib
import hashlib
import io
import time

import click
import pandas
from peewee import SQL, DoesNotExist, fn


__all__ = [
    "Timer",
    "sql_apply_pm",
    "sql_iauname",
    "copy_pandas",
    "get_epoch",
    "set_config_parameter",
    "remove_version",
    "vacuum_outputs",
    "get_configuration_values",
    "vacuum_table",
]


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


def sql_apply_pm(ra_field, dec_field, pmra_field, pmdec_field, epoch_delta, is_pmra_cos=True):
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

    pmra_field = fn.coalesce(pmra_field, 0) / 1000.0 / 3600.0
    pmdec_field = fn.coalesce(pmdec_field, 0) / 1000.0 / 3600.0
    epoch = fn.coalesce(epoch_delta, 0)

    if is_pmra_cos:
        cos_dec = fn.cos(fn.radians(fn.coalesce(dec_field, 0)))
        ra_corr = ra_field + pmra_field / cos_dec * epoch
    else:
        ra_corr = ra_field + pmra_field * epoch

    dec_corr = dec_field + pmdec_field * epoch

    return ra_corr, dec_corr


def sql_iauname(ra_field, dec_field, prefix="SDSS J"):
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
        degs = fn.trunc(dd).cast("INT")
        mins = fn.trunc((fn.abs(dd) - fn.abs(degs)) * 60).cast("INT")

        secs_raw = ((((fn.abs(dd) - fn.abs(degs)) * 60) - mins) * 60).cast("NUMERIC")
        secs = fn.round(secs_raw, round).cast("REAL")

        degs_fmt = "00" if sign is False else "SG00"
        secs_fmt = "00." + "0" * round

        return (
            fn.trim(fn.to_char(degs, degs_fmt))
            .concat(fn.trim(fn.to_char(mins, "00")))
            .concat(fn.trim(fn.to_char(secs, secs_fmt)))
        )

    return (
        SQL("'" + prefix + "'")
        .concat(dd2dms(ra_field / 15))
        .concat(dd2dms(dec_field, round=1, sign=True))
    )


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

    full_table_name = schema + "." + table_name if schema else table_name

    with database.atomic():
        cursor.copy_from(stream, full_table_name, columns=columns, sep=",")

    return df


def get_epoch(xmodel):
    """Returns the epoch for an `.XMatchModel` in Julian years."""

    xmatch = xmodel._meta.xmatch
    fields = xmodel._meta.fields

    if not xmatch.epoch and not xmatch.epoch_column:
        return None

    # If epoch == 0, make it null. This helps with q3c functions.
    if xmatch.epoch:
        epoch = xmatch.epoch
    else:
        epoch = fn.nullif(fields[xmatch.epoch_column], 0)

    if xmatch.epoch_format == "jd":
        epoch = 2000 + (epoch - 2451545.0) / 365.25

    return epoch


@contextlib.contextmanager
def set_config_parameter(database, parameter, new_value, reset=True, log=None):
    """Temporarily a database configuration parameter."""

    new_value = new_value.upper()

    try:
        orig_value = database.execute_sql(f"SHOW {parameter}").fetchone()[0].upper()
        value_changed = orig_value != new_value

        if value_changed:
            database.execute_sql(f"SET {parameter} = {new_value};")
            if log:
                log.debug(f"{parameter} is {new_value}.")

        yield

    finally:
        if reset and value_changed:
            database.execute_sql(f"SET {parameter} = {orig_value};")
            if log:
                log.debug(f"{parameter} reset to {orig_value}.")


def remove_version(
    database, plan, schema="catalogdb", table="catalog", delete_version=True, vacuum=True
):
    """Removes all rows in ``table`` and ``table_to_`` that match a version."""

    models = []
    for path in database.models:
        schema, table_name = path.split(".")
        if schema != schema:
            continue
        if table_name != table and not table_name.startswith(table + "_to_"):
            continue
        model = database.models[path]
        if not model.table_exists() or model._meta.schema != schema:
            continue
        models.append(model)

    if len(models) == 0:
        raise ValueError("No table matches the input parameters.")

    print(
        f"Tables that will be truncated on {plan!r}: "
        + ", ".join(model._meta.table_name for model in models)
    )

    Catalog = database.models[schema + "." + table]
    Version = database.models[schema + ".version"]

    try:
        version_id = Version.get(plan=plan).id
    except DoesNotExist:
        raise ValueError(f"Version {plan!r} does not exist.")

    print(f"version_id={version_id}")

    n_targets = Catalog.select().where(Catalog.version_id == version_id).count()
    print(f"Number of targets found in {Catalog._meta.table_name}: {n_targets:,}")

    if not click.confirm("Do you really want to proceed?"):
        return

    for model in models:
        n_removed = model.delete().where(model.version_id == version_id).execute()
        print(f"{model._meta.table_name}: {n_removed:,} rows removed.")
        if vacuum:
            print("Vacuuming ...")
            vacuum_table(database, f"{model._meta.schema}.{model._meta.table_name}")

    md5 = hashlib.md5(plan.encode()).hexdigest()[0:16]
    for table in database.get_tables(schema=schema):
        if table.endswith(md5):
            print(f"Dropping temporary table {table}.")
            database.execute_sql(f"DROP TABLE IF EXISTS {table};")

    if delete_version:
        Version.delete().where(Version.id == version_id).execute()
        print("Removed entry in 'version'.")


def vacuum_table(database, table_name, vacuum=True, analyze=True, maintenance_work_mem="50GB"):
    """Vacuums and analyses a table."""

    statement = ("VACUUM " if vacuum else "") + ("ANALYZE " if analyze else "") + table_name

    with database.atomic():
        # Change isolation level to allow executing commands such as VACUUM.
        connection = database.connection()
        original_isolation_level = connection.isolation_level
        connection.set_isolation_level(0)

        database.execute_sql(f"SET maintenance_work_mem = {maintenance_work_mem!r}")

        database.execute_sql(statement)

        connection.set_isolation_level(original_isolation_level)


def vacuum_outputs(
    database,
    vacuum=True,
    analyze=True,
    schema="catalogdb",
    table="catalog",
    relational_tables=True,
):
    """Vacuums and analyses the output tables."""

    assert vacuum or analyze, "either vacuum or analyze need to be True."
    assert database.is_connection_usable(), "connection is not usable."

    tables = []
    for full_name in database.models:
        table_schema, table_name = full_name.split(".")
        if table_schema != schema:
            continue
        if table_name != table:
            is_relational = table_name.startswith(table + "_to_")
            if not is_relational or not relational_tables:
                continue
        tables.append(table_name)

    for table_name in tables:
        table_name = table_name if schema is None else schema + "." + table_name
        vacuum_table(database, table_name, vacuum=vacuum, analyze=analyze)


def get_configuration_values(database, parameters):
    """Returns a dictionary of datbase configuration parameter to value."""

    values = {}

    with database.atomic():
        for parameter in parameters:
            value = database.execute_sql(f"SHOW {parameter}").fetchone()[0]
            values[parameter] = value

    return values


def is_view(database, view_name, schema="public", materialized=False):
    """Determines if a view exists."""

    if not materialized:
        query = (
            f"SELECT * FROM pg_views where schemaname = {schema!r} "
            f"AND viewname = {view_name!r};"
        )
    else:
        query = (
            f"SELECT pc.* FROM pg_class pc "
            f"JOIN pg_namespace pn ON pc.relnamespace = pn.oid "
            f"WHERE pn.nspname = {schema!r} AND "
            f"pc.relname = {view_name!r};"
        )

    if len(database.execute_sql(query).fetchall()) > 0:
        return True
    else:
        return False
