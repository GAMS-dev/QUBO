import gamspy as gp
import logging as log
import warnings
import numpy as np
import pandas as pd
from typing import Tuple, Optional

LOG_LEVEL_DICT = {0: log.WARN, 1: log.INFO, 2: log.DEBUG}


def validate_value(value, allowed_values, param_name):
    if value not in allowed_values:
        raise ValueError(f"{param_name} must be one of {allowed_values}")
    return value


class Qubo(gp.Model):
    def __init__(
        self,
        model: gp.Model,
        /,
        name: str = "QUBO",
        problem: str = "MIQCP",
        penalty: int = 1,
        log_on: int = 0,
    ):

        if model.__class__.__name__ != "Model":
            raise Exception(f"Qubo() only accepts a >gamspy.Model< object.")

        self._model: gp.Model = model
        self._modelName: gp.Model = model.name
        self._sense: gp.Sense = model.sense
        self.penalty: gp.Sense = penalty
        # TODO: add other methods
        # self.method: str = validate_value(
        #     method, allowed_values=["classic", "qpu"], param_name="method"
        # )
        self._container: gp.Container = self._run_convert()
        self._q_container = gp.Container()
        self._name = name
        self._problem_type = problem

        if (
            log_level := LOG_LEVEL_DICT.get(
                validate_value(log_on, allowed_values=[0, 1, 2], param_name="log_on"),
                log.WARN,
            )
        ) < log.WARN:
            log.basicConfig(
                filename=f"{self._modelName}_reformulation.log",
                filemode="w",
                format="%(message)s",
                level=log_level,
                force=True,
            )

        warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    def __str__(self) -> str:
        return (
            f"Model QUBO:\n  Problem Type: MIQCP\n  Sense:"
            f" {self._sense}\n  Equations: QUBO_objective"
        )

    def _run_convert(self) -> gp.Container:
        try:
            self._model.solve(
                solver="CONVERT",
                solver_options={"dumpgdx": f"{self._modelName}.gdx", "GDXQuadratic": 1},
            )
        except Exception as e:
            raise Exception(
                f"Error while running the >CONVERT< operation.\nMessage: {e}"
            )

        return gp.Container(load_from=f"{self._modelName}.gdx")

    def transform(self, penalty: int = None):
        if penalty is None:
            penalty = self.penalty

        obj_eq_name: pd.DataFrame = self._container["iobj"].records

        if obj_eq_name is None:
            raise Exception(
                "The objective is not defined using a scalar equation. `iobj` in gdx is empty. Quitting."
            )

        obj_var: pd.DataFrame = self._container[
            "jobj"
        ].records  # fetch the objective variable name
        all_vars: pd.DataFrame = self._container[
            "j"
        ].records  # fetches all variable names
        raw_a: pd.DataFrame = self._container["A"].records  # A coefficients
        eq_data: pd.DataFrame = self._container["e"].records  # fetches equation data

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

        if (
            eq_data["lower"].mod(1).sum(axis=0) > 0
            or eq_data["upper"].mod(1).sum(axis=0) > 0
        ):
            raise Exception(
                "Reformulation with Non-Integer RHS not possible. Quitting."
            )

        bin_vars = self._container[
            "jb"
        ].records  # fetches binary variable names, if any
        int_vars = self._container[
            "ji"
        ].records  # fetches integer variable names, if any
        bin_vars = (
            [] if bin_vars is None else bin_vars["j"].to_list()
        )  # check if any bin_vars are present
        int_vars = (
            [] if int_vars is None else int_vars["j"].to_list()
        )  # check if any int_vars are present
        obj_var = obj_var["j"].to_list()
        all_var_vals = self._container[
            "x"
        ].records  # get all variable values, viz., [level, marginal, lower, upper, scale]

        if (
            len(all_vars) - len(bin_vars) - len(int_vars) != 1
        ):  # Continuous variables are not allowed
            raise Exception("There are continuous variables. Quitting.")

        obj_eq_name = obj_eq_name["i"].to_list()

        check_quad = self._container["ANL"].records

        """
        Check if there are any fixed variables in the gdx, i.e., lb=ub=level of any variable.
        If such variables exist, separate them from the list of non-fixed vairables and treat them as constanst in the objective function.

        We also need to check if the level of variables are set and handle them separately
        """

        vars_with_lower_bounds = {
            var.j: var.lower
            for _, var in all_var_vals.iterrows()
            if (var.lower > 0) and (var.lower != var.upper)
        }  # would only contain, integer varaibles with lower bound defined
        fixed_vars = {
            var.j: var.level
            for _, var in all_var_vals.iterrows()
            if (var.level == var.lower) and (var.level == var.upper)
        }  # check for fixed variables
        fixed_and_lower_bounds = {**vars_with_lower_bounds, **fixed_vars}
        sum_fixed_obj_var_coeffs = 0

        if check_quad is not None:
            rawquad = self._container[
                "Q"
            ].records  # fetch quadratic terms from the original problem, if any.

            if len(int_vars) != 0:
                raise Exception(
                    "Quadratic Program with integer variables are not supported."
                )

            if any(check_quad["j"].isin(fixed_and_lower_bounds.keys())):
                raise Exception(
                    "Quadratic terms with non-zero variable levels are not supported at the moment."
                )

        log.debug("Coefficient matrix: raw_a\n" + raw_a.to_string())
        log.debug("\nEquation Data: eq_data\n" + eq_data.to_string())
        log.debug("\nVariable Data: all_var_vals\n" + all_var_vals.to_string())

        def var_contribution(
            A: pd.DataFrame, vars: dict, cons: Optional[list] = None
        ) -> np.array:
            """
            helper function to calculate the contribution of given variables
            in a constraint or set of constraints

            Args:
                A:      df of coefficients
                vars:   contributing variables
                cons:   participating constraints

            Returns:
                np.array of Total contribution of all variables for that constraint
            """
            cons = slice(None) if cons is None else cons
            coeffs_of_vars_in_constraint = A.loc[cons, vars.keys()].to_numpy()
            lb_var_levels = np.array(list(vars.values())).reshape((len(vars), 1))
            if coeffs_of_vars_in_constraint.size > 0:
                return coeffs_of_vars_in_constraint @ lb_var_levels

            return np.array([0])

        if (
            fixed_and_lower_bounds
        ):  # adjust the rhs of equations when level of variables > 0
            log.info(
                f"\nList of variables with lower bounds:\n{vars_with_lower_bounds}"
            )
            contribution = var_contribution(raw_a, fixed_and_lower_bounds)
            eq_data.loc[:, ["lower", "upper"]] -= contribution
            if fixed_vars:
                log.info(f"\nList of Fixed Variables:\n{fixed_vars}")
                # remove the fixed variables from computation
                bin_vars = [var for var in bin_vars if var not in fixed_vars]
                int_vars = [var for var in int_vars if var not in fixed_vars]
                sum_fixed_obj_var_coeffs += np.ndarray.item(
                    var_contribution(raw_a, fixed_vars, cons=obj_eq_name)
                )
                raw_a.drop(
                    fixed_vars, axis=1, inplace=True
                )  # dropping columns from the coefficient matrix
                fixed_var_vals = all_var_vals[all_var_vals["j"].isin(fixed_vars)].copy(
                    deep=True
                )

            log.debug(
                "\nAfter removing fixed variables and adjusting for non-zero levels: raw_a\n"
                + raw_a.to_string()
            )
            log.debug(
                "\nAfter removing fixed variables and adjusting for non-zero levels: eq_data\n"
                + eq_data.to_string()
            )

        """
        Check if there exist a row in coefficient matrix with all zero values. This can happen if all vars in a constraint are fixed.
        Such row is irrelevant for QUBO and can be dropped out of the matrix and set of constraints.
        """
        redundant_cons = list(raw_a[raw_a.apply(abs).sum(axis=1) == 0].index)
        if len(redundant_cons) > 0:
            log.info(f"\nDropping these redundant constraint: \n{redundant_cons}")
            raw_a.drop(redundant_cons, axis=0, inplace=True)
            eq_data.drop(
                eq_data[eq_data["i"].isin(redundant_cons)].index, axis=0, inplace=True
            )

        def gen_slacks(var_range: float) -> np.array:
            """
            helper function to generate slacks depending on the range of variables or rhs

            Args:
                var_range: upper bound of variable

            Returns:
                Numpy array containing slack co-efficients

            example:
                if var_range=5, then gen_slacks(5) returns [1, 2, 2]
            """
            if var_range >= 1e4:
                raise Exception(
                    "The Upper bound is greater than or equal to 1e+4, Quitting!"
                )

            power = int(np.log2(var_range)) if var_range > 0 else 0
            bounded_coef = var_range - (2**power - 1)
            D_val = [2**i for i in range(power)] + [bounded_coef]
            return np.array(D_val)

        """
        If integer variables exist, convert all integers to binary with '@' as a delimiter of variable names
        If Integer variables with lower bound exist, i.e., lb >=1 and lb!=ub, then convert binary variable for that range
        If these variables contribute to the objective function, their lower bounds are added as a constant
        """
        sum_lower_bound_of_int_vars = 0
        if len(int_vars) != 0:
            if vars_with_lower_bounds:
                sum_lower_bound_of_int_vars += np.ndarray.item(
                    var_contribution(raw_a, vars_with_lower_bounds, obj_eq_name)
                )

            int_var_vals = all_var_vals[all_var_vals["j"].isin(int_vars)]
            int_to_bin_bounds = {
                row["j"]: gen_slacks(row["upper"] - row["lower"])
                for _, row in int_var_vals.iterrows()
            }  # generate coeffs for converted binary vars
            int_bin_vals = pd.DataFrame(columns=["intName", "binName", "value"])
            for var, bin_bounds in int_to_bin_bounds.items():
                for i in range(len(bin_bounds)):  # naming the converted binary vars
                    new_row = pd.DataFrame(
                        {
                            "intName": var,
                            "binName": f"{var}@_bin{i}",
                            "value": bin_bounds[i],
                        },
                        index=[0],
                    )
                    int_bin_vals = (
                        pd.concat([int_bin_vals, new_row], ignore_index=True)
                        if not int_bin_vals.empty
                        else new_row.copy()
                    )

            binName_list = int_bin_vals[
                "binName"
            ].to_list()  # list of all converted binary variable names
            int_bin_name_map = list(
                int_bin_vals[["intName", "binName"]].itertuples(index=None, name=None)
            )
            log.info(
                "\nInteger to Binary Mapping: int_bin_vals\n" + int_bin_vals.to_string()
            )
            int_bin_vals = int_bin_vals.pivot(
                index="intName", columns="binName", values="value"
            ).fillna(
                0
            )  # mapping each binary var to its integer var component
            int_bin_vals = int_bin_vals.reindex(labels=int_vars, axis="index")
            int_bin_vals = int_bin_vals.reindex(labels=binName_list, axis="columns")
            # int_bin_vals.columns = pd.MultiIndex.from_tuples(int_bin_name_map)

            raw_a_int = raw_a[int_vars]
            raw_a_int = raw_a_int.dot(
                int_bin_vals
            )  # updating the "A" coeff matrix with the new coeffs for converted binary vars
            raw_a_rest = raw_a[obj_var + bin_vars]
            raw_a = pd.concat(
                [raw_a_rest, raw_a_int], axis="columns"
            )  # new "A" coeff matrix
            log.info("\nInteger to Binary Mapping: raw_a\n" + raw_a.to_string())
            bin_vars += binName_list  # append the list of original binary variables with the list of converted binary variables

        cons = eq_data[-eq_data["i"].isin(obj_eq_name)].reset_index(
            drop=True
        )  # fetch only the constrainsts and not the objective equation
        nvars = len(bin_vars)
        nslacks = 0
        obj_var_direction = raw_a[obj_var].loc[obj_eq_name].to_numpy()
        obj_var_coeff = raw_a[bin_vars].loc[obj_eq_name].to_numpy()
        if obj_var_direction > 0:
            obj_var_coeff = -1 * obj_var_coeff
        obj = np.zeros((nvars, nvars))
        np.fill_diagonal(obj, obj_var_coeff)

        """
        Pre-processing in case special constraints exist
        Remove the constraint from the "A" matrix and include the special penalty directly in the objective
        Doing so, reduces the number of slack variables used in the final reformulation
        special penalty case 1: sum(x_i | 1 <= i <= n) <= 1 => P*sum(x_i*x_j | i < j)
        special penalty case 2: x_i  + x_j >= 1 => P*(1 - x_i - x_j + x_i*x_j)
        """

        def check_row_entries(df: pd.DataFrame) -> pd.DataFrame:
            """
            helper function to filter DataFrame having either 0 or 1 entries in each row.
            Args:
                df: A Pandas DataFrame.

            Returns:
                Filtered DataFrame with rows having either 0 or 1
            """
            row_contains_only_0s_or_1s = df.isin([0, 1]).all(axis=1)
            return df[row_contains_only_0s_or_1s]

        # Case 1 implementation
        special_cons_case_1_lable = [
            ele.i for _, ele in cons.iterrows() if ele.upper == 1 and ele.lower != 1
        ]
        if special_cons_case_1_lable:
            case1_cons = raw_a[bin_vars].loc[special_cons_case_1_lable]
            case1_cons = check_row_entries(case1_cons.copy())
            case1_cons_index_lable = list(case1_cons.index)
            case1_penalty = case1_cons.to_numpy()
            if case1_penalty.size > 0:
                case1_penalty = (case1_penalty.T @ case1_penalty) / 2
                np.fill_diagonal(case1_penalty, np.zeros((1, len(bin_vars))))
            else:  # if there are no rows with only 0/1 entries
                case1_penalty = np.zeros((nvars, nvars))
            log.debug(f"\nSpecial constraint case 1:\n{special_cons_case_1_lable}")
        else:
            case1_cons_index_lable = []
            case1_penalty = np.zeros((nvars, nvars))

        # Case 2 implementation
        special_cons_case_2_lable = [
            ele.i for _, ele in cons.iterrows() if ele.lower == 1 and ele.upper != 1
        ]
        if special_cons_case_2_lable:
            case2_cons = raw_a[bin_vars].loc[special_cons_case_2_lable]
            case2_cons = check_row_entries(case2_cons.copy())
            case2_cons = case2_cons[case2_cons.sum(axis=1) == 2]
            case2_cons_index_lable = list(case2_cons.index)
            case2_penalty = case2_cons.to_numpy()
            if case2_penalty.size > 0:
                case2_penalty = (case2_penalty.T @ case2_penalty) / 2
                case2_diag = np.diag_indices_from(case2_penalty)
                case2_penalty[case2_diag] *= -2
            else:  # if there are no rows with two 1s in them
                case2_penalty = np.zeros((nvars, nvars))
            log.debug("\nSpecial constraint case 2:\n{special_cons_case_2_lable}")

        else:
            case2_cons_index_lable = []
            case2_penalty = np.zeros((nvars, nvars))

        final_special_cons = case1_cons_index_lable + case2_cons_index_lable
        final_special_penalty = case1_penalty + case2_penalty
        case2_penalty_offset_factor = len(case2_cons_index_lable)

        is_max = True if self._sense == gp.Sense.MAX else False
        P = -1 * penalty if is_max else penalty  # penalty term for classic solvers
        # P_qpu = penalty  # penalty term for qpu, since we always going to minimize the qubo using qpu no treatment is required

        # if self.method == "classic":
        obj += P * final_special_penalty
        # else:
        #     obj = (
        #         P_qpu * final_special_penalty - obj
        #         if is_max
        #         else P_qpu * final_special_penalty + obj
        #     )  # Here, we update the obj and fix the direction to always Minimize for 'qpu'

        cons.drop(cons[cons["i"].isin(final_special_cons)].index, axis=0, inplace=True)
        raw_a.drop(final_special_cons, axis=0, inplace=True)

        A_coeff = raw_a.loc[cons["i"], bin_vars]

        quad = None

        def fetch_quadratic_coeff(raw_df: pd.DataFrame) -> np.array:
            """
            helper function to convert the original Q matrix of the problem to a symmetric matrix

            Args:
                raw_df: Original problem Q data in a pd.DataFrame

            Returns:
                Numpy Q matrix
            """
            raw_df["value"] /= 2
            mask = raw_df["j_1"].astype(str) == raw_df["j_2"].astype(str)
            filtered_quad = raw_df.loc[mask, :].copy()
            raw_df = raw_df.loc[~mask].copy()
            diag_quad = filtered_quad.reset_index(drop=True)
            quad = raw_df.copy(deep=True)
            quad["j_1"], quad["j_2"] = raw_df["j_2"], raw_df["j_1"]
            quad = pd.concat([raw_df, quad, diag_quad], axis=0)
            quad = quad.pivot(index="j_1", columns="j_2", values="value").fillna(0)
            quad = quad.reindex(labels=bin_vars, axis="index")
            quad = quad.reindex(labels=bin_vars, axis="columns")
            return quad.to_numpy()

        if (
            check_quad is not None
        ):  # check if quadratic terms are present in the original problem
            log.debug("\nRaw Q data from GDX: Q\n" + rawquad.to_string())
            rawquad_obj = rawquad[rawquad["i_0"].isin(obj_eq_name)].copy(deep=True)
            if (
                len(rawquad_obj.index) != 0
            ):  # check if quadratic terms exist in the objective function
                rawquad_obj.drop(["i_0"], axis=1, inplace=True)
                quad = fetch_quadratic_coeff(raw_df=rawquad_obj)
                sum_fixed_obj_var_coeffs /= 2

            rawquad_cons = rawquad[
                -rawquad["i_0"].isin(obj_eq_name)
            ]  # non-linear constraints without objective equation
            if len(rawquad_cons.index) != 0:  # non-linear constraints exists
                raise Exception("There are non-linear constraints. Quitting.")

                ### Removed the support for quadratic constraints.
                # mask = rawquad_cons['j_1'].astype(str) == rawquad_cons['j_2'].astype(str)
                # filtered_quad = rawquad_cons.loc[mask, :].reset_index(drop=True)
                # if len(filtered_quad.index) == len(rawquad_cons.index): # if pair-wise quadratic terms are not present, i.e., only terms like x1*x1 and not x1*x2
                #     filtered_quad.drop(['j_2'], axis=1, inplace=True)
                #     filtered_quad['value'] /= 2
                #     filtered_quad = filtered_quad.pivot(index='i_0', columns='j_1', values='value').fillna(0)
                #     if len(filtered_quad.columns) != len(bin_vars):
                #         remain_cols = [var for var in bin_vars if var not in set(filtered_quad.columns)]
                #         filtered_quad[remain_cols] = 0
                #         filtered_quad = filtered_quad[bin_vars]
                #     A_coeff.loc[filtered_quad.index] += filtered_quad.loc[filtered_quad.index].values
                # else: # if pair-wise quadratic terms present, i.e., x1*x2. Quit
                #     raise Exception("There are non-linear constraints. Quitting.")

        if quad is not None:  # add the old quadratic terms/matrix to the new objective
            log.debug("\nUpdate Objective by adding Q: \n" + np.array2string(quad))
            obj += -1 * quad if obj_var_direction > 0 else quad
            log.debug("\nNew Q: \n" + np.array2string(obj))

        def modify_matrix(
            b_vec: np.array,
            rhs: float,
            slacks: np.array,
            A_coeff: pd.DataFrame,
            ele: pd.Series,
            nslacks: int,
        ) -> Tuple[np.array, pd.DataFrame, int]:
            """
            helper function to update the original "A" matrix of coeffs

            Args:
                b_vec: The n*1 vector
                rhs : The Right hand side of a constraint
                slacks: result of gen_slacks()
                A_coeff: "A" matrix
                ele: constraint
                nslacks: number of slacks

            Returns:
                updated b_vec, A_coeff and number of slacks
            """
            con_index = ele.i
            slack_names = [
                f"slack_{con_index}_{i}"
                for i in range(nslacks + 1, nslacks + len(slacks) + 1)
            ]
            A_coeff[slack_names] = 0
            A_coeff.loc[con_index, slack_names] = slacks
            nslacks += len(slacks)
            return np.append(b_vec, [rhs]), A_coeff, nslacks

        def get_lhs_bounds(ele: pd.DataFrame) -> Tuple[float, float]:
            """
            helper function to find the bounds of a constraint

            Args:
                ele: The coefficents of the constraint

            Returns:
                lower_bound, upper_bound
            """
            return ele[ele < 0].sum(), ele[ele > 0].sum()

        b_vec = np.array([])
        log.info("\nFinal Cons: \n" + cons.to_string())
        for _, ele in cons.iterrows():

            if ele.upper == ele.lower:  # equal-to type constraint
                rhs = ele.lower
                lhs_min_lb, lhs_max_ub = get_lhs_bounds(A_coeff.loc[ele.i])
                if (rhs - lhs_min_lb) < 0 or (lhs_max_ub - rhs) < 0:
                    raise Exception(f"Constraint is infeasible: {ele.i}")
                else:
                    b_vec = np.append(b_vec, [rhs])
                    slacks = []  # do not introduce slacks for equality type constraints

            elif ele.upper == np.inf:  # greater than type constraint
                rhs = ele.lower
                _, lhs_max_ub = get_lhs_bounds(A_coeff.loc[ele.i])
                slacks_range = lhs_max_ub - rhs
                if slacks_range > 0:
                    slacks = -1 * gen_slacks(slacks_range)
                    b_vec, A_coeff, nslacks = modify_matrix(
                        b_vec, rhs, slacks, A_coeff, ele, nslacks
                    )
                elif slacks_range == 0:
                    b_vec = np.append(b_vec, [rhs])
                    slacks = []
                else:
                    raise Exception(f"Constraint is infeasible: {ele.i}")

            else:  # less-than type constraint
                rhs = ele.upper
                lhs_min_lb, _ = get_lhs_bounds(A_coeff.loc[ele.i])
                slacks_range = rhs - lhs_min_lb
                if slacks_range > 0:
                    slacks = gen_slacks(slacks_range)
                    b_vec, A_coeff, nslacks = modify_matrix(
                        b_vec, rhs, slacks, A_coeff, ele, nslacks
                    )
                elif slacks_range == 0:
                    b_vec = np.append(b_vec, [rhs])
                    slacks = []
                else:
                    raise Exception(f"Constraint is infeasible: {ele.i}")

        logging_a_mat = A_coeff.unstack().reset_index()
        logging_a_mat = logging_a_mat[logging_a_mat[0] != 0]
        log.debug("\nFinal coefficient matrix: \n" + logging_a_mat.to_string())
        log.debug(f"\nFinal RHS: \n{b_vec}")
        log.debug(f"Constant RHS term: {b_vec.T@b_vec}")
        log.debug(f"Case 2 Offset Penalty Factor: {case2_penalty_offset_factor}")
        log.debug(
            f"Fixed Variable contribution to Objective Function: {sum_fixed_obj_var_coeffs}"
        )
        log.debug(
            f"Integer variable lower bound contribution: {sum_lower_bound_of_int_vars}"
        )

        # A matrix and b_vec are available. Now, for penalization: $(A.X - B)^{2}$  = $(A.X - B)^{T} * (A.X - B)$
        a_mat = A_coeff.to_numpy()
        nvars += nslacks  # increment the total number of variables by total number of slack variable used
        X_diag = np.zeros((nvars, nvars))
        np.fill_diagonal(X_diag, -b_vec.T @ a_mat)
        new_x = a_mat.T @ a_mat + 2 * X_diag

        newobj = np.zeros((nvars, nvars))
        newobj[: len(bin_vars), : len(bin_vars)] = (
            obj  # define the new objective: Q for the qubo
        )

        def qubo_to_ising(Q: dict, offset: float = 0.0) -> Tuple[dict, dict, float]:
            """
            This is the Qubo to Ising Reformulation. Here, the variable X in {-1,1}

            Args:
                    Q: in a form of dict, {(i,j): val}
                    offset: offset from the Qubo reformulation

            Returns:
                    h: the bias vector
                    J: the coupling matrix
                    offset: adjusted offset for the Ising model
            """

            h = {}
            J = {}
            linear_offset = 0.0
            quadratic_offset = 0.0

            for (u, v), bias in Q.items():
                if u == v:
                    if u in h:
                        h[u] += 0.5 * bias
                    else:
                        h[u] = 0.5 * bias
                    linear_offset += bias
                else:
                    if bias != 0.0:
                        J[(u, v)] = 0.25 * bias

                    if u in h:
                        h[u] += 0.25 * bias
                    else:
                        h[u] = 0.25 * bias

                    if v in h:
                        h[v] += 0.25 * bias
                    else:
                        h[v] = 0.25 * bias

                    quadratic_offset += bias

            offset += 0.5 * linear_offset + 0.25 * quadratic_offset

            return h, J, offset

        def qubo_to_maxcut(Q: np.array) -> np.array:
            """
            This is the Qubo to Maxcut Reformulation. This can be used for SDP procedures.

            Args:
                Q: a n x n symmetric numpy matrix

            Returns:
                n x 1 vector associated with the extra variable required in max cut transformation
            """

            return -1 * np.sum(Q, axis=1)

        ### Section to get the qubo for submitting it to the dwave-hybrid method
        # TODO: if self.method == "qpu": ### Plugin QuSol?
        #     Q = (
        #         P_qpu * new_x + newobj
        #     )  # Note: We are always minimizing in the case of qpu and we have already fixed the direction of obj on line 184
        #     offset = P_qpu * b_vec.T @ b_vec + P_qpu * case2_penalty_offset_factor
        #     log.debug(f"\nPenaly: {P_qpu} | Total Offset: {offset}")
        #     if Q.size > 0:
        #         Qdf = pd.DataFrame(
        #             Q, columns=list(A_coeff.columns), index=list(A_coeff.columns)
        #         )
        #         Qdf = Qdf.unstack()
        #         Qdf = Qdf.reset_index()
        #         Qdf_dict = {
        #             (row["level_0"], row["level_1"]): row[0]
        #             for _, row in Qdf.iterrows()
        #         }
        #     else:
        #         Qdf_dict = {}

        ### Section to solve the qubo using miqcp
        # TODO: if self.method in ["classic", "sdp"]:
        const = (
            P * b_vec.T @ b_vec
            + P * case2_penalty_offset_factor
            + sum_fixed_obj_var_coeffs
            + sum_lower_bound_of_int_vars
        )
        log.debug(f"\nPenalty: {P} | Total Offset: {const}\n")
        Q = newobj + P * new_x
        # sparsity = 1.0 - (np.count_nonzero(Q) / float(Q.size))
        ### For producing the qs format for QUBOWL
        # non_zero_indices = np.tril_indices_from(Q)
        # non_zero_values = Q[non_zero_indices]
        # with open('QMat_{self._modelName}.qs', 'w') as f:
        #     f.write(f"{Q.shape[0]} {len(non_zero_values)} {const}\n")
        #     for i, j, value in zip(*non_zero_indices, non_zero_values):
        #         if value != 0:
        #             f.write(f"{i+1} {j+1} {value}\n")
        Qdf = pd.DataFrame(
            Q, columns=list(A_coeff.columns), index=list(A_coeff.columns)
        )
        Qdf = Qdf.unstack()
        Qdf = Qdf.reset_index()
        Qdf = list(Qdf.itertuples(index=None, name=None))

        qconst = gp.Parameter(
            self._q_container,
            "qconst",
            None,
            const,
            description="Constant term to offset the objective value to original",
        )
        if Qdf:
            qi = gp.Set(
                self._q_container,
                "qi",
                records=A_coeff.columns,
                description="QUBO variables",
            )
            qd = gp.Parameter(
                self._q_container,
                "qd",
                [qi, qi],
                records=Qdf,
                description="Q matrix",
            )
            # TODO: add method = sdp if required
            # if self.method == "sdp":
            #     max_cut_var = qubo_to_maxcut(Q)
            #     gp.Parameter(
            #         self._q_container,
            #         "c_sun",
            #         [qi],
            #         records=list(zip(A_coeff.columns, max_cut_var)),
            #         description="extra variable for max cut transformation, ",
            #     )
            # else:
            max_cut_var = np.zeros((len(A_coeff.columns), 1))
            gp.Parameter(
                self._q_container,
                "c_sun",
                [qi],
                records=list(zip(A_coeff.columns, max_cut_var)),
                description="extra variable for max cut transformation, all set to zero since SDP method is not chosen.",
            )
        else:
            qi = gp.Set(
                self._q_container,
                "qi",
                records=["all_variables_fixed"],
                description="QUBO variables",
            )
            qd = gp.Parameter(
                self._q_container,
                "qd",
                [qi, qi],
                records=[("all_variables_fixed", "all_variables_fixed", 0)],
                description="Q matrix",
            )
            # if self.method == "sdp":
            #     raise Exception(
            #         "All variables are fixed. Q matrix is Empty. Quitting SDP SOLVE."
            #     )
        return qd, qi, qconst

    def write_gdx(self, gdxName: str = None):
        if gdxName:
            self._q_container.write(gdxName)
        else:
            self._q_container.write(f"qout_{self._modelName}.gdx")

    def get_q_matrix(self) -> np.ndarray:
        """
        Convenient method to get the Q matrix as numpy array
        """
        try:
            qd: gp.Parameter = self._q_container.getSymbols("qd")[0]
            qd = qd.records.pivot(index="qi_0", columns="qi_1", values="value")

            return qd.to_numpy()
        except KeyError as e:
            raise Exception(
                f"Symbol `qd` not yet present in the Container. Either use .transform() or .solve() to generate the symbol."
            )

    def model(self):
        try:
            qd, qi, qconst = self._q_container.getSymbols(["qd", "qi", "qconst"])
        except Exception as e:
            raise Exception(
                f"Something went from while fetching the q symbols.\nMessage: {e}"
            )

        i = gp.Alias(self._q_container, name="i", alias_with=qi)

        j = gp.Alias(self._q_container, name="j", alias_with=qi)

        x = gp.Variable(
            self._q_container,
            name="x",
            type="binary",
            domain=i,
            description="QUBO variables",
        )

        qubo_obj = gp.Sum(gp.Domain(i, j), x[i] * qd[i, j] * x[j]) + qconst

        return super().__init__(
            container=self._q_container,
            name=self._name,
            problem=self._problem_type,
            sense=self._sense,
            objective=qubo_obj,
        )

    def solve(self, *args, **kwargs):

        if "qd" not in self._q_container.data:
            self.transform()

        if f"{self._name}_objective" not in self._q_container.data:
            self.model()

        try:
            return super().solve(*args, **kwargs)
        except Exception as e:
            raise Exception(f"Something went wrong while solving QUBO.\nMessage: {e}")
