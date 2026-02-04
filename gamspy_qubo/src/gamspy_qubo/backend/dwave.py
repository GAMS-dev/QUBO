from __future__ import annotations

import importlib
from typing import Any

import pandas as pd

from gamspy_qubo._utils import triu
from gamspy_qubo.backend import baseBackend


class DwaveBackend(baseBackend):
    """
    DWave Backend based on `dimod.BinaryQuadraticModel` and `SimulatedAnnealingSampler`.

    Args:
        input_data : A dictionary for providing the necessary input for Dwave backend
    """

    def __init__(self, input_data: dict[str, Any]) -> None:
        self.check_dependencies(
            "DWAVE", {"dwave-ocean-sdk": "dwave.samplers", "dimod": "dimod"}
        )
        self.input_data = input_data

    def solve(self, *args, **kwargs):
        """
        solve the model using Dwave `SimulatedAnnealingSampler` and map the solution back
        to original problem. Users can provide keyword arguments as expected by `SimulatedAnnealingSample.sample()`

        """

        print("\n--- Starting D-Wave (Ocean) Solve ---")
        ut_mat = triu(self.input_data["q_matrix"])
        q_vars = self.input_data["q_vars"]
        matrix_dict = {}
        rows, cols = ut_mat.nonzero()
        for i, j in zip(rows, cols):
            matrix_dict[(q_vars[i], q_vars[j])] = ut_mat[i, j]

        dimod = importlib.import_module("dimod")
        dwave_sampler = importlib.import_module("dwave.samplers")
        bqm = dimod.BinaryQuadraticModel.from_qubo(
            matrix_dict, offset=float(self.input_data["q_const"])
        )
        sampler = dwave_sampler.SimulatedAnnealingSampler()
        kwargs.setdefault("num_reads", 100)
        response = sampler.sample(bqm, **kwargs)
        best_sample = response.first.sample
        best_energy = response.first.energy

        sol = pd.DataFrame(best_sample.items(), columns=["i", "level"])
        # print(f"Best Energy Found: {best_energy}")

        return sol, best_energy
