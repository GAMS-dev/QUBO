from __future__ import annotations

from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

from abc import ABC, abstractmethod
import importlib


class baseBackend(ABC):
    """Abstract base class that provides a structure to different quantum backends."""

    @abstractmethod
    def map_solution(self, solution: pd.DataFrame, **kwargs) -> Any:
        """Maps the solution obtained from the backend back to the original problem."""
        raise NotImplementedError

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
                )
