from __future__ import annotations

import importlib
from abc import ABC, abstractmethod
from typing import Any


class baseBackend(ABC):
    """Abstract base class that provides a structure to different quantum backends."""

    @abstractmethod
    def solve(self, *args, **kwargs) -> Any:
        """Solve the problem using quantum backend"""
        raise NotImplementedError

    @staticmethod
    def check_dependencies(solver_name: str, packages: dict):
        for pkg_name, import_path in packages.items():
            try:
                importlib.import_module(import_path)
            except ImportError:
                raise ImportError(
                    f"The {solver_name} solver requires the '{pkg_name}' package. "
                    f"Install it with: pip install {pkg_name}"
                ) from None
