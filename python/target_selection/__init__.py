# encoding: utf-8

import enlighten

from sdsstools import get_config, get_logger, get_package_version


NAME = 'sdss-target-selection'

config = get_config('target_selection')
log = get_logger(NAME)

manager = enlighten.get_manager()


__version__ = get_package_version(path=__file__, package_name=NAME)
