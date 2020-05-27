#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-26
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import sys
import importlib
import logging

import click
import peewee

from sdssdb.peewee.sdss5db import targetdb as tdb
from sdssdb.peewee.sdss5db.catalogdb import database

import target_selection as tsmod
from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError
from target_selection.xmatch import XMatchPlanner


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
def target_selection(profile, dbname, user, host, port, verbose):
    """Performs tasks related to target selection for SDSS-V."""

    if verbose:
        tsmod.log.set_level(logging.DEBUG)

    # if not connect(profile, dbname, user, host, port):
    #     raise TargetSelectionError('database is not connected.')


@target_selection.command()
@click.argument('TARGETING-PLAN', type=str)
@click.option('--overwrite', is_flag=True,
              help='drop intermediate tables if they exist')
@click.option('--keep', is_flag=True,
              help='keep intermediate tables after loading')
@click.option('--skip-query', is_flag=True,
              help='do not run the query, only load intermediate '
                   'table if it exists')
@click.option('--tile/--no-tile', is_flag=True, default=None,
              help='whether to run the query in chunks tiling the sky')
@click.option('--load/--no-load', is_flag=True, default=True,
              help='whether to load data into targetdb.target')
@click.option('--include', '-i', type=str,
              help='comma-separated carton names to include')
@click.option('--exclude', '-e', type=str,
              help='comma-separated carton names to exclude')
@click.option('--write-table', '-w', is_flag=True,
              help='write intermediate table as a FITS file')
@click.option('--allow-errors', is_flag=True,
              help='continue processing cartons if a carton fails')
def run(targeting_plan, overwrite, keep, tile, load,
        skip_query, include, exclude, write_table, allow_errors):
    """Runs target selection for all cartons."""

    # Reload the carton module. At this point the DB connection should be
    # available and we want to be sure all the cartons are available and
    # can be imported.
    tsmod.__fail_on_carton_import = True
    importlib.reload(sys.modules['target_selection.cartons'])

    carton_classes = {Carton.name: Carton
                      for Carton in BaseCarton.__subclasses__()}

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
            carton = Carton(targeting_plan)

            tsmod.log.header = f'({carton.name}): '
            tsmod.log.info(f'running target selection for '
                           f'carton {carton.name!r}.')

            if carton.check_targets():
                raise ValueError(f'found existing targets for carton '
                                 f'{carton.name!r} with plan {carton.plan!r}.')

            if not skip_query:
                carton.run(tile=tile, overwrite=overwrite)
            else:
                carton.has_run = True
                tsmod.log.debug('skipping query.')

            if load:
                carton.load()
            else:
                tsmod.log.debug('not loading data into targetdb.target.')

            if write_table:
                carton.write_table()

            if not keep:
                tsmod.log.info(f'dropping temporary table {carton.path!r}.')
                carton.drop_table()

        except Exception as ee:
            if allow_errors:
                tsmod.log.error(f'errored processing carton '
                                f'{carton.name}: {ee}')
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
        tdb.Plan.get(label=targeting_plan).delete_instance()
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
