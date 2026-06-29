from __future__ import annotations

import importlib

import pandas as pd

from gamspy_qubo._utils import triu
from gamspy_qubo.backend import baseBackend


class FixstarsBackend(baseBackend):
    """
    Fixstars Amplify Backend based on `amplify.Model` and solver clients.
    """

    def __init__(self, q_matrix, q_variables, q_constant, sense):
        super().__init__(q_matrix, q_variables, q_constant, sense)
        self.check_dependencies("AMPLIFY", {"amplify": "amplify"})

    def solve(self, *args, **kwargs) -> pd.DataFrame:
        """
        Solve the model using Fixstars Amplify and map the solution back
        to original problem.

        Users must provide a solver client via `client` keyword argument
        (e.g., `client=amplify.AmplifyAEClient()`). If no client is provided,
        you can pass a `token` (and optional `timeout`), and it will use
        `AmplifyAEClient` by default.

        Note: Follows the general return structure defined in `baseBackend.solve`
        """
        print("\n--- Starting Fixstars Amplify Solve ---")
        multiplier = -1 if self.sense.value == "MAX" else 1
        ut_mat = triu(self.q_matrix)
        ut_mat *= multiplier

        amplify = importlib.import_module("amplify")

        # 1. Variable array creation
        gen = amplify.VariableGenerator()
        n_vars = len(self.q_variables)
        q = gen.array("Binary", n_vars)

        # 2. Objective function formulation
        f = 0.0
        rows, cols = ut_mat.nonzero()
        for i, j in zip(rows, cols, strict=False):
            f += float(ut_mat[i, j]) * q[i] * q[j]

        f += float(self.q_constant) * multiplier

        # Build the combinatorial optimization model
        model = amplify.Model(f)

        # 3. Solver Client Configuration
        client = kwargs.get("client")
        if client is None:
            token = kwargs.get("token")
            if token:
                # Default to Fixstars Amplify Annealing Engine
                client = amplify.AmplifyAEClient()
                client.token = token
                # Default timeout to 1000 milliseconds (1 second) if not provided
                client.parameters.time_limit_ms = kwargs.get("timeout", 1000)
            else:
                raise ValueError(
                    "Fixstars Amplify requires a configured 'client' (e.g., amplify.AmplifyAEClient()) "
                    "or a 'token' string passed in kwargs."
                )

        # 4. Run the solver
        result = amplify.solve(model, client)

        # Check if the solver successfully returned solutions
        if len(result.solutions) == 0:
            print("No solution found by Fixstars Amplify.")
            return pd.DataFrame(columns=["i", "level"])

        # 5. Map solution back to original variables
        best_values = result.best.values
        q_evaluated = q.evaluate(best_values)

        sol_items = []
        for i, var_name in enumerate(self.q_variables):
            # Evaluate outputs a float array (e.g. 0.0 or 1.0), we cast to int
            sol_items.append((var_name, int(q_evaluated[i])))

        sol = pd.DataFrame(sol_items, columns=["i", "level"])

        return sol
