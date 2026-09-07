import os
from importlib.metadata        import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata        import version as _distribution_version

try:
    __version__ = _distribution_version("scope-qc")
except _PackageNotFoundError:
    __version__ = "0+unknown"

from scope.read_write          import *
from scope.classes_environment import Environment
from scope.classes_system      import System
