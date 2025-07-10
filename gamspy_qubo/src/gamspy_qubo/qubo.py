import re
import gamspy as gp
import logging as log
import numpy as np
import pandas as pd

from gamspy.exceptions import ValidationError, GamspyException

LOG_LEVEL_DICT = {0: log.WARN, 1: log.INFO, 2: log.DEBUG}


def validate_value(value: int, allowed_values: list, param_name: str) -> int:
    if value not in allowed_values:
        raise ValueError(f"{param_name} must be one of {allowed_values}")
    return value


class Qubo(gp.Model):
    """
    QUBO Reformulation for Integer Programs (IPs) modelled in GAMSPy.

    Parameters
    ----------
    model: gp.Model
        GAMSPy generated IP Model.
    name: str | None
        Name for the reformulated model. By default it is "QUBO"
    penalty: int | None
        Set the appropriate penalty to be used to penalize the constraints. By default it is 1
    log_on: int | None
        Enable logging information.
        Options are {0 = WARN, 1 = INFO, 2 = DEBUG}.
        By default it is 0.

    Notes
    -----
    We first use the `CONVERT` solver to generate a standardized version of the original problem.
    One can also pass keyword arguments to the `CONVERT` solver, for. e.g., `options=gp.Options(hold_fixed_variables=True)`
    """

    def __init__(
        self,
        model: gp.Model,
        /,
        name: str = "QUBO",
        penalty: int = 1,
        log_on: int = 0,
        **kwargs,
    ):
        if not isinstance(model, gp.Model):
            raise ValidationError("Qubo() only accepts a >gamspy.Model< object.")

        self._og_model: gp.Model = model
        self._og_modelName: str = model.name
        self._modelName: str = name
        self._work_dir: str = model.container.working_directory
        self._sense: gp.Sense = model.sense
        self.penalty: int = penalty
        self._container: gp.Container = self._run_convert(
            workdir=self._work_dir, **kwargs
        )
        self._q_container = gp.Container(working_directory=self._work_dir)
        self.Q: np.ndarray = None
        self.Qconst: float = None
        self._TRANSFORMATION_COMPLETE = False

        if (
            log_level := LOG_LEVEL_DICT.get(
                validate_value(log_on, allowed_values=[0, 1, 2], param_name="log_on"),
                log.WARN,
            )
        ) < log.WARN:
            log.basicConfig(
                filename=f"{self._og_modelName}_reformulation.log",
                filemode="w",
                format="%(message)s",
                level=log_level,
                force=True,
            )

    def __str__(self) -> str:
        return (
            f"Model {self._modelName}:\n  Problem Type: {self._problem_type}\n  Sense:"
            f" {self._sense}\n  Equations: {self._modelName}_objective"
        )

    def _run_convert(self, workdir, **kwargs) -> gp.Container:
        try:
            self._og_model.solve(
                solver="CONVERT",
                solver_options={
                    "dumpgdx": f"{self._og_modelName}.gdx",
                    "GDXQuadratic": 1,
                    "GDXHessian": 1,
                },
                **kwargs,
            )
        except Exception as e:
            raise GamspyException(
                f"Error while running the >CONVERT< operation.\nMessage: {e}"
            )

        return gp.Container(
            load_from=f"{self._og_modelName}.gdx",
            working_directory=workdir,
        )

    @staticmethod
    def _var_contribution(
        A: pd.DataFrame, vars: dict, cons: list | None = None
    ) -> np.ndarray:
        """
        helper function to calculate the contribution of given variables
        in a constraint or set of constraints

        Args:
            A:      df of coefficients
            vars:   contributing variables
            cons:   participating constraints

        Returns:
            np.ndarray of Total contribution of all variables for that constraint
        """
        cons = slice(None) if cons is None else cons
        coeffs_of_vars_in_constraint = A.loc[cons, vars.keys()].to_numpy()
        lb_var_levels = np.array(list(vars.values())).reshape((len(vars), 1))
        if coeffs_of_vars_in_constraint.size > 0:
            return coeffs_of_vars_in_constraint @ lb_var_levels

        return np.array([0])

    def _check_transformation(self):
        if not self._TRANSFORMATION_COMPLETE:
            raise ValidationError("Run transform first to generate the Q Matrix.")

    def transform(
        self, penalty: int = None
    ) -> tuple[gp.Parameter, gp.Set, gp.Parameter]:
        """
        Method to perform the QUBO reformulation on the provided model.

        Arguments:
            penalty: different penalty can be provided again to generate a new QUBO.

        Returns:
            qd: gp.Parameter: the Q Matrix
            qi: gp.Set: binary variables participating in the QUBO
            qconst: gp.Parameter: the offset calcluated based on the penalty provided.
        """
        if penalty is None:
            penalty = self.penalty

        self._fixed_vars_flag = False
        self._lower_bounded_vars_flag = False
        self._obj_eq_name: pd.DataFrame = self._container["iobj"].records

        if self._obj_eq_name is None:
            raise ValidationError(
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
            raw_a[-raw_a["i"].isin(self._obj_eq_name["i"].tolist())]["value"]
            .mod(1)
            .sum(axis=0)
            > 0
        ):  # floating point coeffs in objective function is accepted
            raise ValidationError(
                "Reformulation with Non-Integer Coefficients not possible. Quitting."
            )

        raw_a = raw_a.pivot(index="i", columns="j", values="value").fillna(
            0
        )  # arranging in a matrix

        if (
            eq_data["lower"].mod(1).sum(axis=0) > 0
            or eq_data["upper"].mod(1).sum(axis=0) > 0
        ):
            raise ValidationError(
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
        self._int_vars_flag = False if len(int_vars) == 0 else True
        obj_var = obj_var["j"].to_list()
        all_var_vals = (
            self._container["x"].records
        )  # get all variable values, viz., [level, marginal, lower, upper, scale]

        if (
            len(all_vars) - len(bin_vars) - len(int_vars) != 1
        ):  # Continuous variables are not allowed
            raise ValidationError("There are continuous variables. Quitting.")

        self._obj_eq_name = self._obj_eq_name["i"].to_list()

        check_quad = self._container["ANL"].records

        """
        Check if there are any fixed variables in the gdx, i.e., lb=ub=level of any variable.
        If such variables exist, separate them from the list of non-fixed vairables and treat them as constanst in the objective function.

        We also need to check if the level of variables are set and handle them separately
        """

        self._vars_with_lower_bounds = {
            var.j: var.lower
            for _, var in all_var_vals.iterrows()
            if (var.lower > 0) and (var.lower != var.upper)
        }  # would only contain, integer varaibles with lower bound defined
        fixed_vars = {
            var.j: var.level
            for _, var in all_var_vals.iterrows()
            if (var.level == var.lower) and (var.level == var.upper)
        }  # check for fixed variables
        fixed_and_lower_bounds = {**self._vars_with_lower_bounds, **fixed_vars}
        sum_fixed_obj_var_coeffs = 0

        if check_quad is not None:
            rawquad = self._container[
                "Q"
            ].records  # fetch quadratic terms from the original problem, if any.

            if self._int_vars_flag:
                raise ValidationError(
                    "Quadratic Program with integer variables are not supported."
                )

            if any(check_quad["j"].isin(fixed_and_lower_bounds.keys())):
                raise ValidationError(
                    "Quadratic terms with non-zero variable levels are not supported at the moment."
                )

        log.debug("Coefficient matrix: raw_a\n" + raw_a.to_string())
        log.debug("\nEquation Data: eq_data\n" + eq_data.to_string())
        log.debug("\nVariable Data: all_var_vals\n" + all_var_vals.to_string())

        if (
            fixed_and_lower_bounds
        ):  # adjust the rhs of equations when level of variables > 0
            log.info(
                f"\nList of variables with lower bounds:\n{self._vars_with_lower_bounds}"
            )
            self._lower_bounded_vars_flag = True
            contribution = self._var_contribution(raw_a, fixed_and_lower_bounds)
            eq_data.loc[:, ["lower", "upper"]] -= contribution
            if fixed_vars:
                self._fixed_vars_flag = True
                log.info(f"\nList of Fixed Variables:\n{fixed_vars}")
                # remove the fixed variables from computation
                bin_vars = [var for var in bin_vars if var not in fixed_vars]
                int_vars = [var for var in int_vars if var not in fixed_vars]
                sum_fixed_obj_var_coeffs += np.ndarray.item(
                    self._var_contribution(raw_a, fixed_vars, cons=self._obj_eq_name)
                )
                raw_a.drop(
                    fixed_vars, axis=1, inplace=True
                )  # dropping columns from the coefficient matrix
                self._fixed_var_vals: pd.DataFrame = all_var_vals[
                    all_var_vals["j"].isin(fixed_vars)
                ].copy(deep=True)

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

        def gen_slacks(var_range: float) -> np.ndarray:
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
                raise ValidationError(
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
        if self._int_vars_flag:
            if self._vars_with_lower_bounds:
                sum_lower_bound_of_int_vars += np.ndarray.item(
                    self._var_contribution(
                        raw_a, self._vars_with_lower_bounds, self._obj_eq_name
                    )
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

            self._binName_list = int_bin_vals[
                "binName"
            ].to_list()  # list of all converted binary variable names
            log.info(
                "\nInteger to Binary Mapping: int_bin_vals\n" + int_bin_vals.to_string()
            )
            int_bin_vals = int_bin_vals.pivot(
                index="intName", columns="binName", values="value"
            ).fillna(0)  # mapping each binary var to its integer var component
            int_bin_vals = int_bin_vals.reindex(labels=int_vars, axis="index")
            int_bin_vals = int_bin_vals.reindex(
                labels=self._binName_list, axis="columns"
            )
            self._int_bin_vals = int_bin_vals

            raw_a_int = raw_a[int_vars]
            raw_a_int = raw_a_int.dot(
                int_bin_vals
            )  # updating the "A" coeff matrix with the new coeffs for converted binary vars
            raw_a_rest = raw_a[obj_var + bin_vars]
            raw_a = pd.concat(
                [raw_a_rest, raw_a_int], axis="columns"
            )  # new "A" coeff matrix
            log.info("\nInteger to Binary Mapping: raw_a\n" + raw_a.to_string())
            bin_vars += self._binName_list  # append the list of original binary variables with the list of converted binary variables

        cons = eq_data[-eq_data["i"].isin(self._obj_eq_name)].reset_index(
            drop=True
        )  # fetch only the constrainsts and not the objective equation
        nvars = len(bin_vars)
        nslacks = 0
        self._obj_var_direction = raw_a[obj_var].loc[self._obj_eq_name].to_numpy()
        obj_var_coeff = raw_a[bin_vars].loc[self._obj_eq_name].to_numpy()
        if self._obj_var_direction > 0:
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
            log.debug(f"\nSpecial constraint case 2:\n{special_cons_case_2_lable}")

        else:
            case2_cons_index_lable = []
            case2_penalty = np.zeros((nvars, nvars))

        final_special_cons = case1_cons_index_lable + case2_cons_index_lable
        final_special_penalty = case1_penalty + case2_penalty
        case2_penalty_offset_factor = len(case2_cons_index_lable)

        is_max = True if self._sense == gp.Sense.MAX else False
        P = -1 * penalty if is_max else penalty  # penalty term for classic solvers
        obj += P * final_special_penalty

        cons.drop(cons[cons["i"].isin(final_special_cons)].index, axis=0, inplace=True)
        raw_a.drop(final_special_cons, axis=0, inplace=True)

        A_coeff = raw_a.loc[cons["i"], bin_vars]

        quad = None

        def fetch_quadratic_coeff(raw_df: pd.DataFrame) -> np.ndarray:
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

        self._quad_val = 0
        if (
            check_quad is not None
        ):  # check if quadratic terms are present in the original problem
            log.debug("\nRaw Q data from GDX: Q\n" + rawquad.to_string())
            rawquad_obj = rawquad[rawquad["i_0"].isin(self._obj_eq_name)].copy(
                deep=True
            )
            if (
                len(rawquad_obj.index) != 0
            ):  # check if quadratic terms exist in the objective function
                rawquad_obj.drop(["i_0"], axis=1, inplace=True)
                quad = fetch_quadratic_coeff(raw_df=rawquad_obj)
                self._quad_val = quad
                sum_fixed_obj_var_coeffs /= 2

            rawquad_cons = rawquad[
                -rawquad["i_0"].isin(self._obj_eq_name)
            ]  # non-linear constraints without objective equation
            if len(rawquad_cons.index) != 0:  # non-linear constraints exists
                raise ValidationError("There are non-linear constraints. Quitting.")
                ### Removed the support for quadratic constraints.

        if quad is not None:  # add the old quadratic terms/matrix to the new objective
            log.debug("\nUpdate Objective by adding Q: \n" + np.array2string(quad))
            obj += -1 * quad if self._obj_var_direction > 0 else quad
            log.debug("\nNew Q: \n" + np.array2string(obj))

        def modify_matrix(
            b_vec: np.ndarray,
            rhs: float,
            slacks: np.ndarray,
            A_coeff: pd.DataFrame,
            ele: pd.Series,
            nslacks: int,
        ) -> tuple[np.ndarray, pd.DataFrame, int]:
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

        def get_lhs_bounds(ele: pd.DataFrame) -> tuple[float, float]:
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
                    raise ValidationError(f"Constraint is infeasible: {ele.i}")
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
                    raise ValidationError(f"Constraint is infeasible: {ele.i}")

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
                    raise ValidationError(f"Constraint is infeasible: {ele.i}")

        logging_a_mat = A_coeff.unstack().reset_index()
        logging_a_mat = logging_a_mat[logging_a_mat[0] != 0]
        log.debug("\nFinal coefficient matrix: \n" + logging_a_mat.to_string())
        log.debug(f"\nFinal RHS: \n{b_vec}")
        log.debug(f"Constant RHS term: {b_vec.T @ b_vec}")
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

        self.Qconst = (
            P * b_vec.T @ b_vec
            + P * case2_penalty_offset_factor
            + sum_fixed_obj_var_coeffs
            + sum_lower_bound_of_int_vars
        )
        log.debug(f"\nPenalty: {P} | Total Offset: {self.Qconst}\n")
        self.Q = newobj + P * new_x
        Qdf = pd.DataFrame(
            self.Q, columns=list(A_coeff.columns), index=list(A_coeff.columns)
        )
        Qdf = Qdf.unstack()
        Qdf = Qdf.reset_index()
        Qdf = list(Qdf.itertuples(index=None, name=None))

        qconst = gp.Parameter(
            self._q_container,
            "qconst",
            None,
            self.Qconst,
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

        self._TRANSFORMATION_COMPLETE = True
        return qd, qi, qconst

    def write_gdx(self, gdxName: str = None):
        if gdxName:
            self._q_container.write(gdxName)
        else:
            self._q_container.write(f"qout_{self._modelName}.gdx")

    def write_qubowl(self) -> None:
        """
        Convenient method to produce QUBO in QUBOWL ready file format (.qs)

        Returns:
                A QMat_<modelName>.qs file accepted by QUBOWL.
        """
        self._check_transformation()

        non_zero_indices = np.tril_indices_from(self.Q)
        non_zero_values = self.Q[non_zero_indices]
        with open(f"QMat_{self._modelName}.qs", "w") as fp:
            fp.write(f"{self.Q.shape[0]} {len(non_zero_values)} {self.Qconst}\n")
            for i, j, value in zip(*non_zero_indices, non_zero_values):
                if value != 0:
                    fp.write(f"{i + 1} {j + 1} {value}\n")

    def _model(self) -> gp.Model:
        try:
            qd, qi, qconst = self._q_container.getSymbols(["qd", "qi", "qconst"])
        except Exception as e:
            raise GamspyException(
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
            name=self._modelName,
            problem="MIQCP",
            sense=self._sense,
            objective=qubo_obj,
        )

    def solve(self, *args, **kwargs) -> pd.DataFrame:
        if not self._TRANSFORMATION_COMPLETE:
            self.transform()

        if f"{self._modelName}_objective" not in self._q_container.data:
            self._model()

        try:
            solved = super().solve(*args, **kwargs)
            self._map_solution()
            return solved
        except Exception as e:
            raise GamspyException(
                f"Something went wrong while solving QUBO.\nMessage: {e}"
            )

    def _map_solution(self) -> None:
        """
        This function maps the QUBO solution to the original Problem
        """
        solveStatus = super().solve_status

        if solveStatus.value in [2, 3, 5, 8]:
            # Continue mapping incumbant solution if solve_status is one of *Interrupt.
            pass

        elif solveStatus.value != 1:
            raise GamspyException("Solver did not yield NormalCompletion.")

        obj_var_coeff: pd.DataFrame = self._q_container[
            f"{self._modelName}_objective_variable"
        ].records
        obj_var = self._container["jobj"].records["j"].values[0]

        all_vars = self._container["j"].records
        original_obj_sym = self._og_model._objective_variable.name

        rem_syms = all_vars[all_vars["uni"] != obj_var]["uni"].to_list()
        optimized_vals = self._q_container["x"].records

        if self._fixed_vars_flag:
            self._fixed_var_vals.rename({"j": "i"}, axis=1, inplace=True)
            optimized_vals = pd.concat(
                [optimized_vals, self._fixed_var_vals], ignore_index=True
            )

        if self._int_vars_flag:
            # check if integer variable exist. If yes, combine and merge the solution of converted binary variables to their integer representation
            int_bin_vals_unstack = self._int_bin_vals.unstack().reset_index()
            int_bin_vals_unstack.drop(
                int_bin_vals_unstack[int_bin_vals_unstack[0] == 0].index, inplace=True
            )
            bin_to_int_vals = pd.merge(
                int_bin_vals_unstack,
                optimized_vals,
                how="left",
                left_on="binName",
                right_on="i",
            )
            bin_to_int_vals["final_level"] = (
                bin_to_int_vals[0] * bin_to_int_vals["level"]
            )
            bin_to_int_vals = (
                bin_to_int_vals.groupby("intName")["final_level"].sum().reset_index()
            )

            original_int_vals = self._container["x"].records
            original_int_vals = original_int_vals[
                original_int_vals["j"].isin(bin_to_int_vals["intName"])
            ].copy(deep=True)
            original_int_vals = pd.merge(
                original_int_vals,
                bin_to_int_vals,
                left_on="j",
                right_on="intName",
                how="left",
            )
            original_int_vals.drop(["level", "intName"], axis=1, inplace=True)
            original_int_vals.rename(
                {"j": "i", "final_level": "level"}, axis=1, inplace=True
            )
            original_int_vals = original_int_vals[
                ["i", "level", "marginal", "lower", "upper", "scale"]
            ]

            optimized_vals.drop(
                optimized_vals[optimized_vals.i.isin(self._binName_list)].index,
                inplace=True,
            )
            optimized_vals = pd.concat(
                [optimized_vals, original_int_vals], ignore_index=True
            )
            if self._vars_with_lower_bounds:
                optimized_vals.loc[
                    optimized_vals["i"].isin(self._vars_with_lower_bounds.keys()),
                    "level",
                ] += list(self._vars_with_lower_bounds.values())

        mapper = {}
        for (
            _,
            ele,
        ) in all_vars.iterrows():  # extract the domain and symbols from the gdx
            domain = re.findall(r"\((.*?)\)", ele["element_text"])
            var_sym = re.findall(r"(.+)(?=\()", ele["element_text"])
            if len(domain) > 0:
                domain = domain[0].strip(r"\'")
                if var_sym[0] not in mapper:
                    mapper[var_sym[0]] = {ele["uni"]: domain}
                else:
                    mapper[var_sym[0]][ele["uni"]] = domain

        for vars, ele in mapper.items():
            newsol = optimized_vals[optimized_vals["i"].isin(ele.keys())].copy(
                deep=True
            )
            newsol["i"] = newsol["i"].map(ele)
            newsol.rename(columns={"i": "QUBO_label"}, inplace=True)

            newsol = newsol[
                ["QUBO_label", "level", "marginal", "lower", "upper", "scale"]
            ]
            split_labels = newsol["QUBO_label"].str.split(",", expand=True)
            split_labels.columns = [
                dom if isinstance(dom, str) else dom.name
                for dom in self._og_model.container[vars].domain
            ]
            newsol = pd.concat([split_labels, newsol], axis=1)
            newsol.drop(["QUBO_label"], axis=1, inplace=True)
            newsol[split_labels.columns] = newsol[split_labels.columns].astype(
                "category"
            )
            self._og_model.container[vars].records = newsol.reset_index(drop=True)

        """
        The code below is required to calculate the contribution of variables towards the objective using the new levels obtained from the QUBO solve.
        Since the QUBO solve returns a different level for the objective variable when the optimal solution is not returned, for example, it includes the penalty for every constraint not satisfied.
        """
        x_l = optimized_vals[optimized_vals["i"].isin(rem_syms)]["level"].to_numpy()
        orig_syms_w_new_levels = (
            optimized_vals[optimized_vals["i"].isin(rem_syms)]
            .set_index("i")["level"]
            .to_dict()
        )
        # at the moment `quad` only contains contribution from the objective row.
        quad_contribution = (
            x_l.T @ self._quad_val @ x_l
            if isinstance(self._quad_val, np.ndarray)
            else 0
        )
        orig_jacobian: pd.DataFrame = self._container["A"].records  # A coefficients
        orig_jacobian = orig_jacobian.pivot(
            index="i", columns="j", values="value"
        ).fillna(0)  # arranging in a matrix
        linear_contribution = self._var_contribution(
            orig_jacobian, orig_syms_w_new_levels, cons=self._obj_eq_name
        ).flatten()[0]
        total_objective_contribution = linear_contribution + quad_contribution
        total_objective_contribution = (
            -1 * total_objective_contribution
            if self._obj_var_direction > 0
            else total_objective_contribution
        )
        obj_var_coeff.loc[0, "level"] = total_objective_contribution
        self._og_model.container[original_obj_sym].records = obj_var_coeff.reset_index(
            drop=True
        )

    @staticmethod
    def check_convexity(Q: np.ndarray) -> str:
        """
        Convenience method to check the convexity of the QUBO using eigenvalues
        """
        eigenvalues = np.linalg.eigvals(Q)

        if np.all(eigenvalues > 0):
            return "Function is strictly convex."
        elif np.all(eigenvalues >= 0):
            return "Function is convex."
        else:
            return "Function is not convex."

    @staticmethod
    def qubo_to_ising(Q: dict, offset: float = 0.0) -> tuple[dict, dict, float]:
        """
        This is the Qubo to Ising Reformulation. Here, the variable X in {-1,1}

        Args:
                Q: in a form of dict, {(i,j): val}
                offset: offset from the Qubo reformulation

        Returns:
                h: the bias vector \
                J: the coupling matrix\
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

    @staticmethod
    def qubo_to_maxcut(Q: np.ndarray) -> np.ndarray:
        """
        This is the Qubo to Maxcut Reformulation. This can be used for SDP procedures.

        Args:
            Q: a n x n symmetric numpy matrix

        Returns:
            n x 1 vector associated with the extra variable required in max cut transformation
        """

        return -1 * np.sum(Q, axis=1)

    @property
    def qubo(self):
        return self.Q
