# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date: 2017-12-05 12:01:21


class TargetSelectionError(Exception):
    """A custom core TargetSelection exception"""

    def __init__(self, message=None):
        message = "There has been an error." if not message else message
        super(TargetSelectionError, self).__init__(message)


class TargetSelectionNotImplemented(TargetSelectionError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):
        message = "This feature is not implemented yet." if not message else message
        super(TargetSelectionNotImplemented, self).__init__(message)


class TargetSelectionMissingDependency(TargetSelectionError):
    """A custom exception for missing dependencies."""

    pass


class XMatchError(TargetSelectionError):
    """An error during cross-matching."""


class TargetSelectionWarning(Warning):
    """Base warning for TargetSelection."""


class TargetSelectionUserWarning(UserWarning, TargetSelectionWarning):
    """The primary warning class."""

    pass


class TargetSelectionDeprecationWarning(TargetSelectionUserWarning):
    """A warning for deprecated features."""

    pass


class TargetSelectionImportWarning(TargetSelectionUserWarning):
    """Warning for import problems."""

    pass
