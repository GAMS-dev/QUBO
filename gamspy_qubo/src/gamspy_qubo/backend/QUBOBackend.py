from abc import ABC, abstractmethod
from typing import Any


class QUBOBackend(ABC):
    """Abstract base class that maps an input to a predetermined output."""

    @abstractmethod
    def map_solution(self, input_data: Any) -> Any:
        """Return a predetermined output for the given input_data."""
        raise NotImplementedError
