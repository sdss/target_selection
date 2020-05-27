# encoding: utf-8

import os

import enlighten

from sdsstools import get_config, get_logger, get_package_version


NAME = 'sdss-target-selection'

config = get_config('target_selection',
                    os.path.join(os.path.dirname(__file__),
                                 'config/target_selection.yml'))

log = get_logger(NAME)

manager = enlighten.get_manager()


# Controls whether a failure importing a carton raises an exception.
__fail_on_carton_import = False


__version__ = get_package_version(path=__file__, package_name=NAME)
