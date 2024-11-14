#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-03
# @Filename: test_target_selection.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations


def test_target_selection_import():
    from target_selection.cartons import BaseCarton

    assert BaseCarton is not None
