#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-26
# @Filename: guide.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from .base import TargetClass


class Guide(TargetClass):

    name = 'guide'
    category = 'guide'

    def build_query(self):
        return super().build_query()
