from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import gamspy as gp

import importlib

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

    def map_solution(self, solution: pd.DataFrame, **kwargs) -> Any:
        container: gp.Container = self.input_data["container"]
        og_model: gp.Model = self.input_data["og_model"]
        orig_obj_var = getattr(og_model, "_objective_variable", None)
        if orig_obj_var is not None:
            original_obj_sym = orig_obj_var.name
        else:
            original_obj_sym = f"{og_model.name}_objective_variable"

        oldvars = container["x"].records
        oldvars.drop(["level"], inplace=True, axis=1)

        res = solution.merge(oldvars, how="right", on="j")
        vardict = container["j"].records
        separate_sym_domain = vardict["element_text"].str.split("(", expand=True)

        if (
            len(separate_sym_domain.columns) == 1
        ):  # check if all variables are flat, i.e., no domain
            vardict["symbol"] = separate_sym_domain[0]
            vardict["domain"] = None
        else:
            vardict[["symbol", "domain"]] = separate_sym_domain[[0, 1]]
            vardict["domain"] = vardict["domain"].str.rstrip(")")
            vardict["domain"] = vardict["domain"].str.strip(r"\'")

        vardict.drop(columns=["element_text"], inplace=True)
        vardict.rename({"uni": "j"}, axis=1, inplace=True)

        final = res.merge(vardict, how="right", on="j")

        for symbol in final["symbol"].unique():
            if symbol == original_obj_sym:
                og_model.container[symbol].records.loc[:, "level"] = kwargs.get(
                    "obj_val"
                )
            else:
                temp = final[final["symbol"] == symbol].reset_index(drop=True)
                temp = temp[["domain", "level", "marginal", "lower", "upper", "scale"]]
                split_labels = temp["domain"].str.split(",", expand=True)
                split_labels.columns = [
                    dom if isinstance(dom, str) else dom.name
                    for dom in og_model.container[symbol].domain
                ]
                temp = pd.concat([split_labels, temp], axis=1)
                temp.drop(["domain"], axis=1, inplace=True)
                temp[split_labels.columns] = temp[split_labels.columns].astype(
                    "category"
                )
                og_model.container[symbol].records = temp.reset_index(drop=True)

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
