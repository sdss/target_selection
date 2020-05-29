
import glob
import importlib
import inspect
import os

from .. import __fail_on_carton_import
from .base import BaseCarton


try:

    # Import cartons so that they can be discovered by
    # calling Carton.__subclasses__().

    exclusions = ['__init__.py', 'base.py']

    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    files = [file_ for file_ in glob.glob('**/*.py', recursive=True)
             if file_ not in exclusions]

    for file_ in files:
        modname = file_[0:-3].replace('/', '.')
        mod = importlib.import_module('target_selection.cartons.' + modname)
        for objname in dir(mod):
            obj = getattr(mod, objname)
            if inspect.isclass(obj) and issubclass(obj, BaseCarton):
                locals().update({objname: obj})

except:

    # Controls whether we raise an error if any of the cartons fails to import
    # This is set to False by default but gets changed to True in the CLI
    # when this module gets reloaded and we want to confirm that all the
    # cartons import correctly.
    if __fail_on_carton_import:
        raise
