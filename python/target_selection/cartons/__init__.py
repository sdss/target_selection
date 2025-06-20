import glob
import importlib
import inspect
import os
import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

from ..exceptions import TargetSelectionImportWarning
from .base import BaseCarton


warnings.filterwarnings("ignore", category=AstropyDeprecationWarning)


# Import cartons so that they can be discovered by
# calling Carton.__subclasses__().

exclusions = ["__init__.py", "base.py"]

cwd = os.getcwd()
os.chdir(os.path.dirname(os.path.realpath(__file__)))

files = [file_ for file_ in glob.glob("*.py") if file_ not in exclusions]

for file_ in files:
    try:
        modname = file_[0:-3].replace("/", ".")
        mod = importlib.import_module("target_selection.cartons." + modname)
        for objname in dir(mod):
            obj = getattr(mod, objname)
            if inspect.isclass(obj) and issubclass(obj, BaseCarton):
                locals().update({objname: obj})

    except Exception as ee:
        warnings.warn(f"cannot import file {file_}: {ee}", TargetSelectionImportWarning)

os.chdir(cwd)
