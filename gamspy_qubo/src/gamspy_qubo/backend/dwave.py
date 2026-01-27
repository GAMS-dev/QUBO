from .QUBOBackend import QUBOBackend
from typing import Any, Dict


class DwaveBackend(QUBOBackend):
    """
    Concrete mapper that returns predetermined outputs based on a mapping dict.
    If the input_data is not found, returns the provided default value.
    """

    def __init__(self, mapping: Dict[Any, Any], default: Any = None) -> None:
        self._mapping = mapping
        self._default = default

    def map_solution(self, input_data: Any) -> Any:
        return self._mapping.get(input_data, self._default)


# Example usage:
# mapper = PredeterminedMapper({"a": 1, "b": 2}, default=0)
# result = mapper.map_solution("a")  # returns 1
# result = mapper.map_solution("z")  # returns 0
