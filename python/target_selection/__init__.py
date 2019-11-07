# encoding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from .utils import get_config, get_logger


NAME = 'target_selection'


# Loads config
config = get_config(NAME)


# Inits the logging system. Only shell logging, and exception and warning catching.
# File logging can be started by calling log.start_file_logger(path).
log = get_logger(NAME)


__version__ = '0.1.0dev'
