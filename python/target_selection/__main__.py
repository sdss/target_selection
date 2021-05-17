#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-05-26
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import importlib
import logging
import os
import sys
import warnings
from glob import glob

import click
import peewee

from sdssdb.peewee.sdss5db import targetdb as tdb
from sdssdb.peewee.sdss5db.catalogdb import database

import target_selection as tsmod
from target_selection.cartons.tools import get_file_carton
from target_selection.exceptions import (TargetSelectionError,
                                         TargetSelectionUserWarning)
from target_selection.xmatch import XMatchPlanner


# Disable warnings during import because the connection may not be working yet.
with warnings.catch_warnings():
    from target_selection.cartons import BaseCarton


def all_subclasses(cls):
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in all_subclasses(c)])


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
        tsmod.log.start_file_logger(os.path.realpath(save_log),
                                    mode='a',
                                    rotating=False)

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
@click.option('--limit', '-l', type=int, default=None,
              help='limit number of targets in the carton')
@click.option('--load/--no-load', is_flag=True, default=True,
              help='whether to load data into targetdb.target')
@click.option('--include', '-i', type=str,
              help='comma-separated carton names to include')
@click.option('--exclude', '-e', type=str,
              help='comma-separated carton names to exclude')
@click.option('--write-table', '-w', is_flag=True,
              help='write table of loaded targets as a FITS file')
@click.option('--exclude-open-fiber', is_flag=True,
              help='do not process open fiber cartons')
def run(targeting_plan, config_file, overwrite, keep, region, load,
        skip_query, include, exclude, write_table, limit, exclude_open_fiber):
    """Runs target selection for all cartons."""

    carton_classes = {Carton.name: Carton
                      for Carton in all_subclasses(BaseCarton)}

    if len(carton_classes) == 0:
        raise TargetSelectionError('no carton classes found.')

    try:
        config_plan = tsmod.config[targeting_plan]
    except KeyError:
        raise TargetSelectionError('cannot find configuration for plan '
                                   f'{targeting_plan}.')

    carton_names = config_plan['cartons']

    if exclude:
        carton_names = [cn for cn in carton_names
                        if cn not in exclude.split(',')]

    if include:
        carton_names = include.split(',')
        for c in carton_names:
            if c not in carton_classes:
                raise ValueError(f'Carton {c} does not exist.')

    Cartons = [carton_classes[carton_name] for carton_name in carton_names]

    if not exclude_open_fiber:
        if 'open_fiber_path' not in config_plan:
            warnings.warn('Plan does not specify an "open_fiber_path".',
                          TargetSelectionUserWarning)
        else:
            open_fiber_path = os.path.expandvars(config_plan['open_fiber_path'])
            if not os.path.exists(open_fiber_path):
                warnings.warn(f'Open fiber path {open_fiber_path} does not exist.',
                              TargetSelectionUserWarning)
            else:
                open_fiber_files = glob(os.path.join(open_fiber_path, '*.fits'))
                OpenFiberCartons = [get_file_carton(fn, '', 'open_fiber', 'open_fiber')
                                    for fn in open_fiber_files]
                if include:
                    OpenFiberCartons = [OFC for OFC in OpenFiberCartons
                                        if OFC.name in include]
                Cartons += OpenFiberCartons
                tsmod.log.info(f'{len(OpenFiberCartons)} open fiber cartons selected.')
    else:
        tsmod.log.info('Excluding open fiber cartons.')

    for Carton in Cartons:

        carton = Carton(targeting_plan, config_file=config_file)

        tsmod.log.header = f'({carton.name}): '
        tsmod.log.info(f'Running target selection for '
                       f'carton {carton.name!r}.')

        if carton.check_targets() and not overwrite:
            raise ValueError(f'Found existing targets for carton '
                             f'{carton.name!r} with plan {carton.plan!r}.')

        if not skip_query:
            carton.run(query_region=(region or None),
                       overwrite=overwrite,
                       limit=limit)
            if write_table:
                carton.write_table()
        else:
            tsmod.log.debug('Skipping query.')

        if load:
            carton.load(overwrite=overwrite)
        else:
            tsmod.log.info('Not loading data into targetdb.target.')

        if write_table:
            carton.write_table(mode='targetdb')

        if not keep:
            tsmod.log.info(f'Dropping temporary table {carton.path!r}.')
            carton.drop_table()

    tsmod.log.header = ''


@target_selection.command('load-files')
@click.option('--program', type=str, default='open_fiber', show_default=True)
@click.option('--category', type=str, default='open_fiber', show_default=True)
@click.argument('PLAN', type=str, nargs=1)
@click.argument('FILES', nargs=-1, required=True)
def load_files(plan, files, program, category):
    """Loads a carton from a FITS file."""

    for fn in files:
        Carton = get_file_carton(fn, '', category, program)
        carton = Carton(plan)
        carton.run(overwrite=True)
        carton.load(overwrite=True)


@target_selection.command()
@click.argument('TARGETING-PLAN', type=str)
@click.option('--tables', is_flag=True,
              help='also remove intermediate tables')
def clear(targeting_plan, tables):
    """Clear all data for a target selection plan."""

    if tables:
        # Cartons = BaseCarton.__subclasses__()
        # find grandchildren classes as well
        Cartons = all_subclasses(BaseCarton)
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
