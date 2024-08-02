#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-08-01
# @Filename: test_skies.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import pathlib

from typing import TYPE_CHECKING

import numpy
import polars
import pytest

from sdsstools import yanny

from target_selection.skies import is_valid_sky


if TYPE_CHECKING:
    from sdssdb.connection import PeeweeDatabaseConnection


@pytest.fixture()
def database():
    from sdssdb.peewee.sdss5db import database

    if not database.connected or database.profile != "tunnel_operations":
        database.set_profile("tunnel_operations")
        assert database.connected

    yield database


@pytest.fixture()
def sky_candidates():
    confSummary_path = pathlib.Path(__file__).parent / "data" / "confSummary-15130.par"
    confSummary = yanny(str(confSummary_path))["FIBERMAP"]

    not_on_target = confSummary[
        (confSummary["fiberType"] != "METROLOGY")
        & (confSummary["on_target"] == 0)
        & (confSummary["ra"] >= 0)
        & (confSummary["dec"] > -999)
    ]

    yield numpy.array(not_on_target[["ra", "dec"]].tolist(), dtype="f8")


def test_database(database: PeeweeDatabaseConnection):
    assert database.connected
    assert database.profile == "tunnel_operations"
    assert database.database == "sdss5db"


def test_is_valid_sky(database: PeeweeDatabaseConnection, sky_candidates: numpy.ndarray):
    mask, df = is_valid_sky(
        sky_candidates,
        database,
        catalogues=["gaia_dr3_source", "twomass_psc"],
        return_dataframe=True,
    )

    assert isinstance(mask, numpy.ndarray)
    assert isinstance(df, polars.DataFrame)

    assert (~mask).sum() == 30
    assert df.filter(polars.col.gaia_dr3_source.not_()).height == 21
