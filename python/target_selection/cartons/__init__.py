
# flake8: noqa

from .. import __fail_on_carton_import
from .base import BaseCarton
print(__fail_on_carton_import)
try:

    # Import cartons so that they can be discovered by
    # calling Carton.__subclasses__().
    # Eventually maybe automate this with a glob.
    from .guide import GuideCarton

except:

    # Controls whether we raise an error if any of the cartons fails to import
    # This is set to False by default but gets changed to True in the CLI
    # when this module gets reloaded and we want to confirm that all the
    # cartons import correctly.
    if __fail_on_carton_import:
        raise
