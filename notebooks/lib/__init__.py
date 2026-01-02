#!/usr/bin/env python3
"""
Lightweight package entrypoint, mirroring the cooltools pattern:
- lib.core: plotting/config utilities only (no data load)
- lib.read_data_basic: shared reference tables
- lib.read_data_hic: Hi-C loaders and fountain stacks
- lib.read_data_epigenetics: epigenetic stacks and annotations

Importing `lib` stays cheap; submodules are loaded lazily on access.
Prefer `from lib.core import *` or `import lib.read_data_hic` to scope imports.
"""

from importlib import import_module

import warnings
# Ignore warnings for deprecated proplot parameters
warnings.filterwarnings("ignore")

# Logging: default to INFO for data-loading visibility
import logging
import sys
logging.basicConfig(
    stream=sys.stdout, 
    # level=logging.INFO,
    format="[%(asctime)s] %(name)s:%(levelname)s: %(message)s",
    # format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def log_step(message: str) -> None:
    """Lightweight helper to report major data-loading steps."""
    logger.info(message)


_SUBMODULES = ("core", "read_data_basic", "read_data_hic", "read_data_epigenetics", "read_nonstandard")

__all__ = ["logger", "log_step"] + list(_SUBMODULES)


def __getattr__(name: str):
    if name in _SUBMODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__} has no attribute {name}")


def __dir__():
    return sorted(set(globals().keys()) | set(_SUBMODULES))
