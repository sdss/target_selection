#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-24
# @Filename: base.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import logging

from sdssdb.peewee.sdss5db import catalogdb, database, targetdb

from .. import log


class TargetClassMeta(type):

    __require_overload__ = ['build_query']

    def __new__(meta, name, bases, dct):

        if bases:
            assert dct.get('name', None) is not None, 'name is not defined.'

            # Make sure we have overloaded the necessary methods. We do this
            # manually because we are already using one metaclass and cannot
            # use abc.
            for method in meta.__require_overload__:
                assert method in dct, f'overloading of {method!r} is required.'

        return super(TargetClassMeta, meta).__new__(meta, name, bases, dct)


class TargetClass(metaclass=TargetClassMeta):
    """A base class for target classes.

    Attributes
    ----------
    name : str
        The name of the target class (required).
    cadence : str
        The label of the candece rule for this target class.
    category : str
        The category of targets for this target class.
    survey : str
        The survey associated with this target class.
    orm : str
        The ORM library to be used, ``peewee`` or ``sqlalchemy``.

    """

    name = None
    cadence = None
    category = None
    survey = None

    orm = 'peewee'

    def __init__(self):

        if not self.orm == 'peewee':
            raise NotImplementedError('not implemented for SQLAlchemy.')

        assert database.connected, 'database is not connected.'

        self._check_targetdb()

    def log(self, message, level=logging.INFO):
        """Logs a message with a header of the current target class name."""

        message = f'[{self.name.upper()}]: {message}'
        log.log(level, message)

    def _check_targetdb(self):
        """Checks that target, cadence, and survey are present in targetdb."""

        assert self.category in targetdb.Category.select(targetdb.Category.label)

    def build_query(self):
        """Builds and returns the query.

        Returns
        -------
        The ORM query for the target class. Note that this must be the
        *un-executed* query, which will be executed in `.run`.

        """

        pass

    def run(self, model=None, tile=False, tile_size=20):
        """Executes the query and stores the results.

        Parameters
        ----------
        model
            The model for the table to which to append the results. If the
            table exists, the results with be appended.


        """

        pass
