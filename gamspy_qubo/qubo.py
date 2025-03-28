import gamspy as gp
import logging as lg
import warnings
import numpy as np
import pandas as pd

LOG_LEVEL_DICT = {0: lg.WARN, 1: lg.INFO, 2: lg.DEBUG}


def validate_value(value, allowed_values, param_name):
    if value not in allowed_values:
        raise ValueError(f"{param_name} must be one of {allowed_values}")
    return value


class Qubo:
    def __init__(
        self,
        model: gp.Model,
        /,
        penalty: int = 1,
        method: str = "classic",
        solver: str = "cplex",
        maxIter: int = 1,
        timeLimit: int = 10,
        log_on: int = False,
        examinerOn: bool = False,
        get_Q: bool = False,
    ):

        if model.__class__.__name__ != "Model":
            raise Exception(f"Qubo() only accepts a >gamspy.Model< object.")

        self.model: gp.Model = model
        self.modelName: gp.Model = model.name
        self.m: gp.Container = model.container
        self.sense: gp.Sense = model.sense
        self.penalty: int = penalty
        self.method: str = validate_value(
            method, allowed_values=["classic", "qpu"], param_name="method"
        )
        self.solver: str = solver
        self.maxIter: int = maxIter
        self.timeLimit: int = timeLimit
        self.log_on: int = validate_value(
            log_on, allowed_values=[0, 1, 2], param_name="log_on"
        )
        self.examinerOn: bool = examinerOn
        self.get_Q: bool = get_Q
        self.container: gp.Container = self.run_convert()

        if (log_level := LOG_LEVEL_DICT.get(self.log_on, lg.WARN)) < lg.WARN:
            lg.basicConfig(
                filename=f"{self.modelName}_reformulation.log",
                filemode="w",
                format="%(message)s",
                level=log_level,
                force=True,
            )

        warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    def run_convert(self) -> gp.Container:
        self.model.solve(
            solver="CONVERT",
            solver_options={"dumpgdx": f"{self.modelName}.gdx", "GDXQuadratic": 1},
        )

        return gp.Container(load_from=f"{self.modelName}.gdx")

    def transform(self):
        obj_eq_name = self.container["iobj"].records

        if obj_eq_name is None:
            raise Exception(
                "The objective is not defined using a scalar equation. `iobj` in gdx is empty. Quitting."
            )

        obj_var = self.container["jobj"].records  # fetch the objective variable name
        all_vars = self.container["j"].records  # fetches all variable names
        raw_a = self.container["A"].records  # A coefficients
        eq_data = self.container["e"].records  # fetches equation data

        if (
            raw_a[-raw_a["i"].isin(obj_eq_name["i"].tolist())]["value"]
            .mod(1)
            .sum(axis=0)
            > 0
        ):  # floating point coeffs in objective function is accepted
            raise Exception(
                "Reformulation with Non-Integer Coefficients not possible. Quitting."
            )

        raw_a = raw_a.pivot(index="i", columns="j", values="value").fillna(
            0
        )  # arranging in a matrix

        lg.debug("Coefficient matrix: raw_a\n" + raw_a.to_string())
        lg.info(f"Objective var  =\n{obj_var}")
        lg.info(f"All variables  =\n{all_vars}")
        lg.info(f"eq_data  =\n{eq_data}")
