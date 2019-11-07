# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date:   2017-12-05 12:01:21
# @Last modified by:   Brian Cherinka
# @Last Modified time: 2017-12-05 12:19:32

from __future__ import print_function, division, absolute_import


class Target_selectionError(Exception):
    """A custom core Target_selection exception"""

    def __init__(self, message=None):

        message = 'There has been an error' \
            if not message else message

        super(Target_selectionError, self).__init__(message)


class Target_selectionNotImplemented(Target_selectionError):
    """A custom exception for not yet implemented features."""

    def __init__(self, message=None):

        message = 'This feature is not implemented yet.' \
            if not message else message

        super(Target_selectionNotImplemented, self).__init__(message)


class Target_selectionAPIError(Target_selectionError):
    """A custom exception for API errors"""

    def __init__(self, message=None):
        if not message:
            message = 'Error with Http Response from Target_selection API'
        else:
            message = 'Http response error from Target_selection API. {0}'.format(message)

        super(Target_selectionAPIError, self).__init__(message)


class Target_selectionApiAuthError(Target_selectionAPIError):
    """A custom exception for API authentication errors"""
    pass


class Target_selectionMissingDependency(Target_selectionError):
    """A custom exception for missing dependencies."""
    pass


class Target_selectionWarning(Warning):
    """Base warning for Target_selection."""


class Target_selectionUserWarning(UserWarning, Target_selectionWarning):
    """The primary warning class."""
    pass


class Target_selectionSkippedTestWarning(Target_selectionUserWarning):
    """A warning for when a test is skipped."""
    pass


class Target_selectionDeprecationWarning(Target_selectionUserWarning):
    """A warning for deprecated features."""
    pass
