#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-26
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import importlib
import logging
import sys
import warnings

import click
import peewee

from sdssdb.peewee.sdss5db import targetdb as tdb
from sdssdb.peewee.sdss5db.catalogdb import database

import target_selection as tsmod
from target_selection.exceptions import TargetSelectionError
from target_selection.xmatch import XMatchPlanner


# Disable warnings during import because the connection may not be working yet.
with warnings.catch_warnings():
    from target_selection.cartons import BaseCarton


def connect(profile=None, dbname=None, user=None, host=None, port=None):
    """Connects the database."""

    if profile:
        database.set_profile(profile)

    if dbname or user or host or port:
        database.connect(dbname=dbname, user=user, host=host, port=port)

    return database.connected


@click.group()
@click.option('--profile', '-p', type=str, default=None)
@click.option('--dbname', '-d', type=str, default=None)
@click.option('--user', '-u', type=str, default=None)
@click.option('--host', '-h', type=str, default=None)
@click.option('--port', '-P', type=int, default=None)
@click.option('--verbose', '-v', is_flag=True,
              help='outputs extra debug information')
@click.option('--save-log', type=str, default=None,
              help='saves the log to a file.')
def target_selection(profile, dbname, user, host, port, verbose, save_log):
    """Performs tasks related to target selection for SDSS-V."""

    if verbose:
        tsmod.log.set_level(logging.DEBUG)

    if save_log:
        tsmod.log.start_file_logger(save_log, rotating=False)

    if not connect(profile, dbname, user, host, port):
        raise TargetSelectionError('database is not connected.')

    # Reload the cartons now that we have a connection.
    importlib.reload(sys.modules['target_selection.cartons'])


@target_selection.command()
@click.argument('TARGETING-PLAN', nargs=1, type=str)
@click.option('--config-file', type=click.Path(exists=True, dir_okay=False),
              help='the file to read. Defaults to the internal '
                   'configuration file.')
@click.option('--overwrite', is_flag=True,
              help='drop intermediate tables if they exist')
@click.option('--keep', is_flag=True,
              help='keep intermediate tables after loading')
@click.option('--skip-query', is_flag=True,
              help='do not run the query, only load intermediate '
                   'table if it exists')
@click.option('--region', '-r', type=float, nargs=3, default=None,
              help='the region (ra, dec, radius) to query')
@click.option('--load/--no-load', is_flag=True, default=True,
              help='whether to load data into targetdb.target')
@click.option('--include', '-i', type=str,
              help='comma-separated carton names to include')
@click.option('--exclude', '-e', type=str,
              help='comma-separated carton names to exclude')
@click.option('--write-table', '-w', is_flag=True,
              help='write table of loaded targets as a FITS file')
@click.option('--allow-errors', is_flag=True,
              help='continue processing cartons if a carton fails')
def run(targeting_plan, config_file, overwrite, keep, region, load,
        skip_query, include, exclude, write_table, allow_errors):
    """Runs target selection for all cartons."""

    carton_classes = {Carton.name: Carton
                      for Carton in BaseCarton.__subclasses__()}

    if len(carton_classes) == 0:
        raise TargetSelectionError('no carton classes found.')

    try:
        config_plan = tsmod.config[targeting_plan]
    except KeyError:
        raise TargetSelectionError('cannot find configuration for plan '
                                   f'{targeting_plan}.')

    cartons = config_plan['cartons']

    if exclude:
        cartons = [carton for carton in cartons if carton not in exclude]

    if include:
        cartons = [carton for carton in cartons if carton in include]

    for carton_name in cartons:

        try:

            Carton = carton_classes[carton_name]
            carton = Carton(targeting_plan, config_file=config_file)

            tsmod.log.header = f'({carton.name}): '
            tsmod.log.info(f'running target selection for '
                           f'carton {carton.name!r}.')

            if carton.check_targets():
                raise ValueError(f'found existing targets for carton '
                                 f'{carton.name!r} with plan {carton.plan!r}.')

            if not skip_query:
                carton.run(query_region=(region or None), overwrite=overwrite)
            else:
                carton.has_run = True
                tsmod.log.debug('skipping query.')

            if load:
                carton.load()
            else:
                tsmod.log.debug('not loading data into targetdb.target.')

            if write_table:
                carton.write_table(mode='targetdb')

            if not keep:
                tsmod.log.info(f'dropping temporary table {carton.path!r}.')
                carton.drop_table()

        except Exception as ee:

            if allow_errors:
                if carton_name not in carton_classes:
                    tsmod.log.error(f'no carton class found for {carton_name}')
                else:
                    tsmod.log.error(f'errored processing carton '
                                    f'{carton_name}: {ee}')
            else:
                raise

    tsmod.log.header = ''


@target_selection.command()
@click.argument('TARGETING-PLAN', type=str)
@click.option('--tables', is_flag=True,
              help='also remove intermediate tables')
def clear(targeting_plan, tables):
    """Clear all data for a target selection plan."""

    if tables:
        Cartons = BaseCarton.__subclasses__()
        for Carton in Cartons:
            Carton(targeting_plan).drop_table()
    try:
        tdb.Version.get(plan=targeting_plan,
                        target_selection=True).delete_instance()
    except peewee.DoesNotExist:
        pass


@target_selection.command()
@click.argument('XMATCH-PLAN', type=str)
@click.option('--file', type=click.Path(exists=True, dir_okay=False),
              help='the file to read. Defaults to the internal '
                   'configuration file.')
def xmatch(xmatch_plan, file):
    """Runs catalogue cross-matching from a configuration file."""

    xmatch = XMatchPlanner.read(database, xmatch_plan, config_file=file)
    xmatch.run()


if __name__ == '__main__':

    target_selection()
