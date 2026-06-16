"""scikit-build-core dynamic-metadata provider for the netgen-mesher wheel.

Computes the PEP 440 package version from the git history, reusing the logic in
``tests/utils.py`` (the single source of truth, also used by the CI ``--check-pip``
and ``--wait-pip`` helpers).

Referenced from ``pyproject.toml``::

    [tool.scikit-build.metadata.version]
    provider = "netgen_metadata"
    provider-path = "tests"

scikit-build-core inserts ``provider-path`` (``tests``) onto ``sys.path`` and
imports this module, so ``import utils`` resolves to ``tests/utils.py``.
"""

from __future__ import annotations

import os
import sys

# scikit-build-core only keeps provider-path on sys.path while importing this
# module, so resolve and import ``utils`` (tests/utils.py) now, at import time.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

__all__ = ["dynamic_metadata", "get_requires_for_dynamic_metadata"]


def __dir__() -> list[str]:
    return __all__


def dynamic_metadata(field: str, settings: dict | None = None) -> str:
    if field != "version":
        msg = "Only the 'version' field is supported"
        raise ValueError(msg)
    if settings:
        msg = "No inline configuration is supported"
        raise ValueError(msg)

    return utils.get_version(os.getcwd())


def get_requires_for_dynamic_metadata(_settings: dict | None = None) -> list[str]:
    # get_version() only needs git and the standard library.
    return []
