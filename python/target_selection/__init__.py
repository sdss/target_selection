# encoding: utf-8

from sdsstools import get_config, get_logger, get_package_version


NAME = 'sdss-target-selection'

# Loads config
config = get_config(NAME)

# Inits the logging system. Only shell logging, and exception and warning catching.
# File logging can be started by calling log.start_file_logger(path).
log = get_logger(NAME)


__version__ = get_package_version(path=__file__, package_name=NAME)
