$offEolCom
$offListing
$eolcom #

$setArgs modelName modelType direction obj penalty method_name solver_name maxIter_num timeLimit_num num_of_threads log_on examiner_on

$ifThenE.defmethod sameas('%method_name%','') 
    $$set method classic # classic method by default 
$elseIfE.defmethod sameas('%method_name%','classic')
    $$set method classic
$elseIfE.defmethod sameas('%method_name%','qpu')
    $$set method qpu
$else.defmethod
    $$abort 'Wrong method chosen, [qpu, classic]'
$endIf.defmethod

$ifThenE.defsolver sameas('%solver_name%','') 
    $$set solver cplex # default: cplex
$else.defsolver
    $$set solver %solver_name%
$endIf.defsolver

$ifThenE.defmaxIter sameas('%maxIter_num%','') 
    $$set maxIter 1 # 1 by default 
$else.defmaxIter
    $$eval maxIter %maxIter_num%
$endIf.defmaxIter

$ifThenE.defTimeLimit sameas('%timeLimit_num%','') 
    $$set timeLimit 10 # 10 by default 
$else.defTimeLimit
    $$eval timeLimit %timeLimit_num%
$endIf.defTimeLimit

$ifThenE.defnumThreads sameas('%num_of_threads%','') 
    $$set num_threads 1 # 1 by default
$else.defnumThreads
    $$eval num_threads %num_of_threads%
$endIf.defnumThreads

$ifThenE.defExaminer sameas('%examiner_on%','1') 
    $$set check_examiner 1
    $$echo examineGamsPoint 1 > examiner.opt
$else.defExaminer
    $$eval check_examiner 0
$endIf.defExaminer

$onEcho > convert.opt
dumpgdx %modelName%.gdx
GDXQuadratic 1
$offEcho

option %modelType%=convert;

%modelName%.optfile = 1;

Solve %modelName% use %modelType% %direction% %obj%; # dumping the problem data in a gdx file

* QUBO Reformulations
EmbeddedCode Python: 
import pandas as pd
from gams import transfer as gt
import numpy as np
from math import log2, ceil
import re
import warnings
import logging

log_level_dict = {"0": logging.WARN, "1": logging.INFO, "2": logging.DEBUG}
log_level = log_level_dict.get('%log_on%', logging.WARN)

if log_level < 30:
        logging.basicConfig(filename='%modelName%_reformulation.log', filemode='w', format='%(message)s', level=log_level, force=True)

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

is_max = "x" in "%direction%".lower()
gdxFile = r"%modelName%.gdx"
gdx_data = gt.Container(gdxFile)
obj_eq_name = gdx_data['iobj'].records

if obj_eq_name is None:
    raise Exception("The objective is not defined using a scalar equation. `iobj` in gdx is empty. Quitting.")

obj_var = gdx_data['jobj'].records # fetch the objective variable name
all_vars = gdx_data['j'].records # fetches all variable names
rawA = gdx_data['A'].records # A coefficients 
eq_data = gdx_data['e'].records # fetches equation data

if rawA[-rawA['i'].isin(obj_eq_name['i'].tolist())]['value'].mod(1).sum(axis=0) > 0: # floating point coeffs in objective function is accepted
    raise Exception("Reformulation with Non-Integer Coefficients not possible. Quitting.")

rawA = rawA.pivot(index="i", columns="j", values="value").fillna(0) # arranging in a matrix

if eq_data['lower'].mod(1).sum(axis=0) > 0 or eq_data['upper'].mod(1).sum(axis=0) > 0:
    raise Exception("Reformulation with Non-Integer RHS not possible. Quitting.")

bin_vars = gdx_data['jb'].records # fetches binary variable names, if any
int_vars = gdx_data['ji'].records # fetches integer variable names, if any
bin_vars = [] if bin_vars is None else bin_vars['j'].to_list() # check if any bin_vars are present
int_vars = [] if int_vars is None else int_vars['j'].to_list() # check if any int_vars are present
obj_var = obj_var['j'].to_list()
all_var_vals = gdx_data['x'].records # get all variable values, viz., [level, marginal, lower, upper, scale]

if len(all_vars) - len(bin_vars) - len(int_vars) != 1: # Continuous variables are not allowed
    raise Exception("There are continuous variables. Quitting.")

obj_eq_name = obj_eq_name['i'].to_list()

checkQuad = gdx_data['ANL'].records
   
"""
Check if there are any fixed variables in the gdx, i.e., lb=ub=level of any variable.
If such variables exist, separate them from the list of non-fixed vairables and treat them as constanst in the objective function.

We also need to check if the level of variables are set and handle them separately
"""

vars_with_lower_bounds = {var.j: var.lower for _, var in all_var_vals.iterrows() if (var.lower > 0) and (var.lower != var.upper)} # would only contain, integer varaibles with lower bound defined
fixed_vars = {var.j: var.level for _, var in all_var_vals.iterrows() if (var.level == var.lower) and (var.level == var.upper)} # check for fixed variables
fixed_and_lower_bounds = {**vars_with_lower_bounds, **fixed_vars}
sum_fixed_obj_var_coeffs = 0

if checkQuad is not None:
    rawquad = gdx_data['Q'].records # fetch quadratic terms from the original problem, if any.

    if len(int_vars) != 0:
        raise Exception("Quadratic Program with integer variables are not supported.")

    if any(checkQuad['j'].isin(fixed_and_lower_bounds.keys())):
        raise Exception("Quadratic terms with non-zero variable levels are not supported at the moment.")

logging.debug("Coefficient matrix: rawA\n"+rawA.to_string())
logging.debug("\nEquation Data: eq_data\n"+eq_data.to_string())
logging.debug("\nVariable Data: all_var_vals\n"+all_var_vals.to_string())

def varContribution(A, vars, cons=slice(None)):
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
    coeffs_of_vars_in_constraint = A.loc[cons, vars.keys()].to_numpy()
    lb_var_levels = np.array(list(vars.values())).reshape((len(vars), 1))
    if coeffs_of_vars_in_constraint.size > 0:
        return coeffs_of_vars_in_constraint@lb_var_levels
    else:
        return np.array([0])

if fixed_and_lower_bounds: # adjust the rhs of equations when level of variables > 0
    logging.info(f"\nList of variables with lower bounds:\n{vars_with_lower_bounds}")
    contribution = varContribution(rawA, fixed_and_lower_bounds)
    eq_data.loc[:,['lower', 'upper']] -= contribution
    if fixed_vars:
        logging.info(f"\nList of Fixed Variables:\n{fixed_vars}")
        # remove the fixed variables from computation
        bin_vars = [var for var in bin_vars if var not in fixed_vars.keys()]
        int_vars = [var for var in int_vars if var not in fixed_vars.keys()]
        sum_fixed_obj_var_coeffs += np.ndarray.item(varContribution(rawA, fixed_vars, cons=obj_eq_name))
        rawA.drop(fixed_vars, axis=1, inplace=True) # dropping columns from the coefficient matrix
        fixed_var_vals = all_var_vals[all_var_vals['j'].isin(fixed_vars)].copy(deep=True)

    logging.debug("\nAfter removing fixed variables and adjusting for non-zero levels: rawA\n"+rawA.to_string())
    logging.debug("\nAfter removing fixed variables and adjusting for non-zero levels: eq_data\n"+eq_data.to_string())  

"""
Check if there exist a row in coefficient matrix with all zero values. This can happen if all vars in a constraint are fixed.
Such row is irrelevant for QUBO and can be dropped out of the matrix and set of constraints.
"""
redundant_cons = list(rawA[rawA.apply(abs).sum(axis=1)==0].index)
if len(redundant_cons) > 0:
    logging.info(f"\nDropping these redundant constraint: \n{redundant_cons}")
    rawA.drop(redundant_cons, axis=0, inplace=True)
    eq_data.drop(eq_data[eq_data['i'].isin(redundant_cons)].index, axis=0, inplace=True)


def genslacks(var_range):
    """
    helper function to generate slacks depending on the range of variables or rhs

    Args:
        var_range: upper bound of variable

    Returns: 
        Numpy array containing slack co-efficients

    example: 
        if var_range=5, then genslacks(5) returns [1, 2, 2]
    """
    if var_range >= 1e+4:
        raise Exception("The Upper bound is greater than or equal to 1e+4, Quitting!")
    
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
        sum_lower_bound_of_int_vars += np.ndarray.item(varContribution(rawA, vars_with_lower_bounds, obj_eq_name))

    int_var_vals = all_var_vals[all_var_vals['j'].isin(int_vars)]
    int_to_bin_bounds = {row['j']: genslacks(row['upper']-row['lower']) for _,row in int_var_vals.iterrows()} # generate coeffs for converted binary vars
    int_bin_vals = pd.DataFrame(columns=['intName', 'binName', 'value'])
    for var, bin_bounds in int_to_bin_bounds.items():
        for i in range(len(bin_bounds)): # naming the converted binary vars
            int_bin_vals = pd.concat([int_bin_vals, pd.DataFrame({'intName': var, 'binName': f"{var}@_bin{i}", 'value': bin_bounds[i]}, index=[0])], ignore_index=True)
    
    binName_list = int_bin_vals['binName'].to_list() # list of all converted binary variable names 
    int_bin_name_map = list(int_bin_vals[['intName', 'binName']].itertuples(index=None, name=None))
    logging.info("\nInteger to Binary Mapping: int_bin_vals\n"+int_bin_vals.to_string())
    int_bin_vals = int_bin_vals.pivot(index='intName', columns='binName', values='value').fillna(0) # mapping each binary var to its integer var component
    int_bin_vals = int_bin_vals.reindex(labels=int_vars, axis='index')
    int_bin_vals = int_bin_vals.reindex(labels=binName_list, axis='columns')
    # int_bin_vals.columns = pd.MultiIndex.from_tuples(int_bin_name_map)

    rawA_int = rawA[int_vars]
    rawA_int = rawA_int.dot(int_bin_vals) # updating the "A" coeff matrix with the new coeffs for converted binary vars
    rawA_rest = rawA[obj_var+bin_vars]
    rawA = pd.concat([rawA_rest, rawA_int], axis='columns') # new "A" coeff matrix
    logging.info("\nInteger to Binary Mapping: rawA\n"+rawA.to_string())
    bin_vars += binName_list # append the list of original binary variables with the list of converted binary variables 

cons = eq_data[-eq_data['i'].isin(obj_eq_name)].reset_index(drop=True) # fetch only the constrainsts and not the objective equation
nvars = len(bin_vars)
nslacks = 0
obj_var_direction = rawA[obj_var].loc[obj_eq_name].to_numpy()
obj_var_coeff = rawA[bin_vars].loc[obj_eq_name].to_numpy()
if obj_var_direction > 0:
    obj_var_coeff = -1*obj_var_coeff
obj = np.zeros((nvars, nvars))
np.fill_diagonal(obj, obj_var_coeff)

"""
Pre-processing in case special constraints exist
Remove the constraint from the "A" matrix and include the special penalty directly in the objective
Doing so, reduces the number of slack variables used in the final reformulation
special penalty case 1: sum(x_i | 1 <= i <= n) <= 1 => P*sum(x_i*x_j | i < j)
special penalty case 2: x_i  + x_j >= 1 => P*(1 - x_i - x_j + x_i*x_j)
"""

def check_row_entries(df):
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
special_cons_case_1_lable = [ele.i for _, ele in cons.iterrows() if ele.upper == 1 and ele.lower != 1]
if special_cons_case_1_lable:
    case1_cons = rawA[bin_vars].loc[special_cons_case_1_lable]
    case1_cons = check_row_entries(case1_cons.copy())
    case1_cons_index_lable = list(case1_cons.index)
    case1_penalty = case1_cons.to_numpy()
    if case1_penalty.size > 0:    
        case1_penalty = (case1_penalty.T@case1_penalty)/2
        np.fill_diagonal(case1_penalty, np.zeros((1, len(bin_vars))))
    else: # if there are no rows with only 0/1 entries
        case1_penalty = np.zeros((nvars, nvars))
    logging.debug(f"\nSpecial constraint case 1:\n{special_cons_case_1_lable}")
else:
    case1_cons_index_lable = []
    case1_penalty = np.zeros((nvars, nvars))

# Case 2 implementation
special_cons_case_2_lable = [ele.i for _, ele in cons.iterrows() if ele.lower == 1 and ele.upper != 1]
if special_cons_case_2_lable:
    case2_cons = rawA[bin_vars].loc[special_cons_case_2_lable]
    case2_cons = check_row_entries(case2_cons.copy())
    case2_cons = case2_cons[case2_cons.sum(axis=1)==2]
    case2_cons_index_lable = list(case2_cons.index)
    case2_penalty = case2_cons.to_numpy()
    if case2_penalty.size > 0: 
        case2_penalty = (case2_penalty.T@case2_penalty)/2
        case2_diag = np.diag_indices_from(case2_penalty)
        case2_penalty[case2_diag] *= -2
    else: # if there are no rows with two 1s in them
        case2_penalty = np.zeros((nvars, nvars))
    logging.debug("\nSpecial constraint case 2:\n{special_cons_case_2_lable}")

else:
    case2_cons_index_lable = []
    case2_penalty = np.zeros((nvars, nvars))

final_special_cons = case1_cons_index_lable + case2_cons_index_lable
final_special_penalty = case1_penalty + case2_penalty
case2_penalty_offset_factor = len(case2_cons_index_lable)

P = -1 * %penalty% if is_max else %penalty% # penalty term for classic solvers
P_qpu = %penalty% # penalty term for qpu, since we always going to minimize the qubo using qpu no treatment is required

if '%method%' == 'classic':
    obj += P*final_special_penalty
else:
    obj = P_qpu*final_special_penalty - obj if is_max else P_qpu*final_special_penalty + obj # Here, we update the obj and fix the direction to always Minimize for 'qpu'

cons.drop(cons[cons["i"].isin(final_special_cons)].index, axis=0, inplace=True)
rawA.drop(final_special_cons, axis=0, inplace=True)

A_coeff = rawA.loc[cons['i'], bin_vars]

quad = None

def fetchQuadraticCoeff(rawdf):
    """
    helper function to convert the original Q matrix of the problem to a symmetric matrix

    Args:
        rawdf: Original problem Q data in a pd.DataFrame

    Returns: 
        Numpy Q matrix
    """
    rawdf['value'] /= 2
    mask = rawdf['j_1'].astype(str) == rawdf['j_2'].astype(str)
    filtered_quad = rawdf.loc[mask, :].copy()
    rawdf = rawdf.loc[~mask].copy()
    diag_quad = filtered_quad.reset_index(drop=True)
    quad = rawdf.copy(deep=True)
    quad['j_1'], quad['j_2'] = rawdf['j_2'], rawdf['j_1']
    quad  = pd.concat([rawdf, quad, diag_quad], axis=0)
    quad = quad.pivot(index="j_1", columns="j_2", values="value").fillna(0)
    quad = quad.reindex(labels=bin_vars, axis='index')
    quad = quad.reindex(labels=bin_vars, axis='columns')
    return quad.to_numpy()
    

if checkQuad is not None: # check if quadratic terms are present in the original problem
    logging.debug("\nRaw Q data from GDX: Q\n"+rawquad.to_string())
    rawquad_obj = rawquad[rawquad['i_0'].isin(obj_eq_name)].copy(deep=True)
    if len(rawquad_obj.index) != 0:  # check if quadratic terms exist in the objective function
        rawquad_obj.drop(['i_0'],axis=1,inplace=True)
        quad = fetchQuadraticCoeff(rawdf=rawquad_obj)
        sum_fixed_obj_var_coeffs /= 2
    
    rawquad_cons = rawquad[-rawquad['i_0'].isin(obj_eq_name)] # non-linear constraints without objective equation
    if len(rawquad_cons.index) != 0: # non-linear constraints exists 
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

if quad is not None: # add the old quadratic terms/matrix to the new objective
    logging.debug("\nUpdate Objective by adding Q: \n"+np.array2string(quad))
    obj += -1*quad if obj_var_direction > 0 else quad
    logging.debug("\nNew Q: \n"+np.array2string(obj))


def modify_matrix(Bvec, rhs, slacks, A_coeff, ele, nslacks):
    """
    helper function to update the original "A" matrix of coeffs

    Args:   
        Bvec: The n*1 vector 
        rhs : The Right hand side of a constraint
        slacks: result of genslacks()
        A_coeff: "A" matrix
        ele: constraint
        nslacks: number of slacks

    Returns: 
        updated Bvec, A_coeff and number of slacks 
    """
    con_index = ele.i
    slack_names = [f'slack_{con_index}_{i}' for i in range(nslacks+1, nslacks + len(slacks) + 1)]
    A_coeff[slack_names] = 0
    A_coeff.loc[con_index, slack_names] = slacks
    nslacks += len(slacks)
    return np.append(Bvec, [rhs]), A_coeff, nslacks

def getLhsBounds(ele):
    """
    helper function to find the bounds of a constraint

    Args:
        ele: The coefficents of the constraint

    Returns:
        lower_bound, upper_bound 
    """
    return ele[ele < 0].sum(), ele[ele > 0].sum()

Bvec = np.array([])
logging.info(f"\nFinal Cons: \n"+cons.to_string())
for _, ele in cons.iterrows():

    if ele.upper == ele.lower: # equal-to type constraint
        rhs = ele.lower
        lhs_min_lb, lhs_max_ub = getLhsBounds(A_coeff.loc[ele.i])
        if (rhs - lhs_min_lb) < 0 or (lhs_max_ub - rhs) < 0:
            raise Exception(f"Constraint is infeasible: {ele.i}")
        else:
            Bvec = np.append(Bvec, [rhs])
            slacks = [] # do not introduce slacks for equality type constraints

    elif ele.upper == np.inf: # greater than type constraint
        rhs = ele.lower
        _, lhs_max_ub = getLhsBounds(A_coeff.loc[ele.i])
        slacks_range = lhs_max_ub - rhs
        if slacks_range > 0:
            slacks = -1*genslacks(slacks_range)
            Bvec, A_coeff, nslacks = modify_matrix(Bvec, rhs, slacks, A_coeff, ele, nslacks)
        elif slacks_range == 0:
            Bvec = np.append(Bvec, [rhs])
            slacks = []
        else:
            raise Exception(f"Constraint is infeasible: {ele.i}")

    else: # less-than type constraint
        rhs = ele.upper
        lhs_min_lb, _ = getLhsBounds(A_coeff.loc[ele.i])
        slacks_range = rhs - lhs_min_lb
        if slacks_range > 0:
            slacks = genslacks(slacks_range)
            Bvec, A_coeff, nslacks = modify_matrix(Bvec, rhs, slacks, A_coeff, ele, nslacks)
        elif slacks_range == 0:
            Bvec = np.append(Bvec, [rhs])
            slacks = []
        else:
            raise Exception(f"Constraint is infeasible: {ele.i}")


logging_a_mat = A_coeff.unstack().reset_index()
logging_a_mat = logging_a_mat[logging_a_mat[0]!= 0]
logging.debug(f"\nFinal coefficient matrix: \n"+logging_a_mat.to_string())
logging.debug(f"\nFinal RHS: \n{Bvec}")
logging.debug(f"Constant RHS term: {Bvec.T@Bvec}")
logging.debug(f"Case 2 Offset Penalty Factor: {case2_penalty_offset_factor}")
logging.debug(f"Fixed Variable contribution to Objective Function: {sum_fixed_obj_var_coeffs}")
logging.debug(f"Integer variable lower bound contribution: {sum_lower_bound_of_int_vars}")

# A matrix and Bvec are available. Now, for penalization: $(A.X - B)^{2}$  = $(A.X - B)^{T} * (A.X - B)$
Amat = A_coeff.to_numpy()
nvars += nslacks # increment the total number of variables by total number of slack variable used
X_diag = np.zeros((nvars, nvars))
np.fill_diagonal(X_diag,-Bvec.T@Amat)
newX = Amat.T@Amat + 2*X_diag

newobj = np.zeros((nvars,nvars))
newobj[:len(bin_vars), :len(bin_vars)] = obj # define the new objective: Q for the qubo


def qubo_to_ising(Q, offset=0.0):
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
                h[u] += .5 * bias
            else:
                h[u] = .5 * bias
            linear_offset += bias
        else:
            if bias != 0.0:
                J[(u, v)] = .25 * bias

            if u in h:
                h[u] += .25 * bias
            else:
                h[u] = .25 * bias

            if v in h:
                h[v] += .25 * bias
            else:
                h[v] = .25 * bias

            quadratic_offset += bias

    offset += .5 * linear_offset + .25 * quadratic_offset
    
    return h, J, offset


def qubo_to_maxcut(Q):
    """
    This is the Qubo to Maxcut Reformulation. This can be used for SDP procedures.
    
    Args: 
        Q: a n x n symmetric numpy matrix

    Returns:
        n x 1 vector associated with the extra variable required in max cut transformation
    """

    return -1*np.sum(Q, axis=1)

### Section to get the qubo for submitting it to the dwave-hybrid method
if '%method%' == 'qpu':
    Q = P_qpu*newX + newobj # Note: We are always minimizing in the case of qpu and we have already fixed the direction of obj on line 184
    offset = P_qpu*Bvec.T@Bvec + P_qpu*case2_penalty_offset_factor
    logging.debug(f"\nPenaly: {P_qpu} | Total Offset: {offset}")
    if Q.size > 0:    
        Qdf = pd.DataFrame(Q, columns=list(A_coeff.columns), index=list(A_coeff.columns))
        Qdf = Qdf.unstack()
        Qdf = Qdf.reset_index()
        Qdf_dict = {(row['level_0'], row['level_1']): row[0] for _,row in Qdf.iterrows()}
    else:
        Qdf_dict = {}

### Section to solve the qubo using miqcp
elif '%method%' in ['classic', 'sdp']:
    const = P*Bvec.T@Bvec + P*case2_penalty_offset_factor + sum_fixed_obj_var_coeffs + sum_lower_bound_of_int_vars
    logging.debug(f"\nPenalty: {P} | Total Offset: {const}\n")
    Q = newobj + P*newX
    #sparsity = 1.0 - (np.count_nonzero(Q) / float(Q.size))
    #np.savetxt(fr'%modelName%_p{P}_sprs{sparsity}.csv', Q, delimiter=",")
    # logging.debug("\nFinal Q:\n"+np.array2string(Q))
    ### For producing the qs format for QUBOWL
    # non_zero_indices = np.tril_indices_from(Q)
    # non_zero_values = Q[non_zero_indices]
    # with open('QMat_%modelName%.qs', 'w') as f:
    #     f.write(f"{Q.shape[0]} {len(non_zero_values)} {const}\n")
    #     for i, j, value in zip(*non_zero_indices, non_zero_values):
    #         if value != 0:
    #             f.write(f"{i+1} {j+1} {value}\n")
    Qdf = pd.DataFrame(Q, columns=list(A_coeff.columns), index=list(A_coeff.columns))
    Qdf = Qdf.unstack()
    Qdf = Qdf.reset_index()
    Qdf = list(Qdf.itertuples(index=None, name=None))
    
    m = gt.Container()
    qconst = gt.Parameter(m, "qconst", None, const, description="Constant term to offset the objective value to original")
    if Qdf:
        qi = gt.Set(m, "qi", records=A_coeff.columns, description="QUBO variables")
        qd = gt.Parameter(m, "qd", [qi, qi], records=Qdf, description="Q matrix")
        if '%method%' == 'sdp':
            max_cut_var = qubo_to_maxcut(Q)
            gt.Parameter(m, "c_sun", [qi], records=list(zip(A_coeff.columns, max_cut_var)), description="extra variable for max cut transformation, ")
        else:
            max_cut_var = np.zeros((len(A_coeff.columns),1))
            gt.Parameter(m, "c_sun", [qi], records=list(zip(A_coeff.columns, max_cut_var)), description="extra variable for max cut transformation, all set to zero since SDP method is not chosen.")
    else:
        qi = gt.Set(m, "qi", records=['all_variables_fixed'], description="QUBO variables")
        qd = gt.Parameter(m, "qd", [qi, qi], records=[('all_variables_fixed', 'all_variables_fixed', 0)], description="Q matrix")
        if '%method%' == 'sdp':
            raise Exception("All variables are fixed. Q matrix is Empty. Quitting SDP SOLVE.")
    m.write('qout_%modelName%.gdx')
endEmbeddedCode

abort$execError 'Error occured in Reformulation!';

* Solve the QUBO using a miqcp
$ifThenE.method sameas('%method%','classic')
embeddedCode GAMS: lo=%gams.lo%
set qi;

Alias (qi,i,j);

Parameter qd(i,j);

Scalar qconst;

Parameter c_sun(i);

$gdxload 'qout_%modelName%.gdx' qconst qi qd c_sun

Equation defQubo;

Binary variable x(i);

Free variable z;

defQubo.. z =E= sum((i,j), x(i)*qd(i,j)*x(j)) + qconst;

Model qubo /all/;

option miqcp=%solver%, threads=%num_threads%;

qubo.reslim = %timeLimit%;

Solve qubo %direction% z using miqcp;

execute_unload 'qout_%modelName%.gdx';
endembeddedCode
    
* Map the QUBO's solution to the original problem
EmbeddedCode Python:
qout = gt.Container('qout_%modelName%.gdx')
obj_var_coeff = qout['z'].records
obj_var = gdx_data['jobj'].records['j'].values[0]

all_vars = gdx_data['j'].records
original_obj_sym = all_vars[all_vars['uni'] == obj_var]['element_text'].values[0]
rem_syms = all_vars[all_vars['uni'] != obj_var]['uni'].to_list()
optimized_vals = qout['x'].records

if fixed_vars:
    fixed_var_vals.rename({'j': 'i'}, axis=1, inplace=True)
    optimized_vals = pd.concat([optimized_vals, fixed_var_vals], ignore_index=True)

if len(int_vars) != 0: # check if integer variable exist. If yes, combine and merge the solution of converted binary variables to their integer representation
    int_bin_vals_unstack = int_bin_vals.unstack().reset_index()
    int_bin_vals_unstack.drop(int_bin_vals_unstack[int_bin_vals_unstack[0] == 0].index, inplace=True)
    bin_to_int_vals = pd.merge(int_bin_vals_unstack, optimized_vals, how="left", left_on="binName", right_on="i")
    bin_to_int_vals['final_level'] = bin_to_int_vals[0]*bin_to_int_vals['level']
    bin_to_int_vals = bin_to_int_vals.groupby('intName')['final_level'].sum().reset_index()
    
    original_int_vals = gdx_data['x'].records
    original_int_vals = original_int_vals[original_int_vals['j'].isin(bin_to_int_vals['intName'])].copy(deep=True)
    original_int_vals = pd.merge(original_int_vals, bin_to_int_vals, left_on="j", right_on="intName", how="left")
    original_int_vals.drop(['level', 'intName'], axis=1, inplace=True)
    original_int_vals.rename({'j':'i', 'final_level': 'level'}, axis=1, inplace=True)
    original_int_vals = original_int_vals[["i", "level", "marginal", "lower", "upper", "scale"]]

    optimized_vals.drop(optimized_vals[optimized_vals.i.isin(binName_list)].index, inplace=True)
    optimized_vals = pd.concat([optimized_vals, original_int_vals], ignore_index=True)
    if vars_with_lower_bounds:
        optimized_vals.loc[optimized_vals["i"].isin(vars_with_lower_bounds.keys()), 'level'] += list(vars_with_lower_bounds.values())


mapper = {}
for index, ele in all_vars.iterrows(): # extract the domain and symbols from the gdx
    domain  = re.findall(r'\((.*?)\)', ele['element_text'])
    var_sym = re.findall(r'(.+)(?=\()', ele['element_text'])
    if len(domain) > 0:
        domain = domain[0].strip(r"\'")
        if var_sym[0] not in mapper.keys():
            mapper[var_sym[0]] = {ele['uni'] : domain}
        else:
            mapper[var_sym[0]][ele['uni']] = domain

for vars, ele in mapper.items():
    newsol = optimized_vals.copy(deep=True)
    newsol = newsol[newsol['i'].isin(ele.keys())]
    newsol.index = newsol['i']
    newsol.drop(['i'], axis=1, inplace=True)
    newsol['oldlabel'] = None
    
    for key, val in ele.items():
        newsol.loc[key, 'oldlabel'] = val
    newsol = newsol[['oldlabel', 'level', 'marginal', 'lower','upper','scale']]
    split_labels = newsol['oldlabel'].str.split(',', expand=True)
    new_columns = ['newlabel{}'.format(i+1) for i in range(split_labels.shape[1])]
    split_labels.columns = new_columns
    newsol = pd.concat([split_labels, newsol], axis=1)
    newsol.drop(['oldlabel'], axis=1, inplace=True)
    gams.set(vars, list(newsol.itertuples(index=None,name=None)))

if len(mapper) == 0: # support for flat gms file
    opt_list = optimized_vals.copy(deep=True)
    opt_list.drop(['i'], axis=1, inplace=True)
    flat_var_mapping = all_vars.set_index('uni')['element_text'].to_dict()
    for key, val in flat_var_mapping.items():
        if val != original_obj_sym:
            gams.set(val, list(opt_list.iloc[optimized_vals[optimized_vals['i']==key].index].itertuples(index=None,name=None)))

"""
The code below is required to calculate the contribution of variables towards the objective using the new levels obtained from the QUBO solve.
Since the QUBO solve returns a different level for the objective variable when the optimal solution is not returned, for example, it includes the penalty for every constraint not satisfied.
"""
x_l = optimized_vals[optimized_vals['i'].isin(rem_syms)]['level'].to_numpy()
orig_syms_w_new_levels = optimized_vals[optimized_vals['i'].isin(rem_syms)].set_index('i')['level'].to_dict()
quad_contribution = x_l.T@quad@x_l if checkQuad is not None else 0 #at the moment `quad` only contains contribution from the objective row.
orig_jacobian = gdx_data['A'].records # A coefficients 
orig_jacobian = orig_jacobian.pivot(index="i", columns="j", values="value").fillna(0) # arranging in a matrix
linear_contribution = varContribution(orig_jacobian, orig_syms_w_new_levels, cons=obj_eq_name).flatten()[0]
total_objective_contribution = linear_contribution + quad_contribution
total_objective_contribution = -1*total_objective_contribution if obj_var_direction > 0 else total_objective_contribution    
obj_var_coeff.loc[0, 'level'] = total_objective_contribution

gams.set(original_obj_sym, list(obj_var_coeff.itertuples(index=None, name=None)))
endEmbeddedCode

* Solve the QUBO on the QPU
$elseIfE.method sameas('%method%','qpu')
embeddedCode Python:
try:
    from dwave.system import LeapHybridSampler
except:
    raise Exception("Package 'Dwave.system' is not available")

try:
    from dimod import BinaryQuadraticModel
except:
    raise Exception("Package 'Dimod' is not available")

print("Submitting problem to QPU")
solver = LeapHybridSampler.default_solver
sampler = LeapHybridSampler(solver=solver)
from datetime import datetime

cur_date = datetime.now().strftime("%d.%m.%Y")

if Q.size > 0:# check if QMatrix exist, if not skip model generation
    bqm = BinaryQuadraticModel.from_qubo(Qdf_dict, offset=offset)
    max_iter = %maxIter%
    all_sol = pd.DataFrame()
    sampler_label_base = f'QUBO-%modelName%-{cur_date}'

    for counter in range(max_iter): # For repeated sampling, a default option is not avaiable for the chosen sampler
        print("Solve number:", counter+1)
        answer = sampler.sample(bqm, label=f'{sampler_label_base}-{counter+1}', time_limit=%timeLimit%)
        new_df = answer.to_pandas_dataframe()
        all_sol = pd.concat([all_sol, new_df], ignore_index=True)
        
    best_val = all_sol['energy'].min()
    best_sol = all_sol[all_sol['energy']==best_val].reset_index(drop=True)
    best_sol.drop(['energy', 'num_occurrences'], axis=1, inplace=True)
    best_sol = best_sol.loc[0].to_dict()
        
    print(f"Best Objective Value: {best_val}")
    print(f"All solutions:\n{all_sol}")

    sample = pd.DataFrame(best_sol.items(), columns=['j', 'level'])
else:
    warnings.warn("QMatrix size is zero. Skipping model creation.")
    sample = pd.DataFrame(columns=['j', 'level'])
    best_val = 0

if fixed_vars:
    sample = pd.concat([sample, fixed_var_vals[['j', 'level']]], ignore_index=True)

if len(int_vars) != 0: # check if integer variable exist. If yes, combine and merge the solution of converted binary variables to their integer representation
    int_bin_vals_unstack = int_bin_vals.unstack().reset_index()
    int_bin_vals_unstack.drop(int_bin_vals_unstack[int_bin_vals_unstack[0] == 0].index, inplace=True)
    
    bin_to_int_vals = pd.merge(int_bin_vals_unstack, sample, how="left", left_on="binName", right_on="j")
    bin_to_int_vals['final_level'] = bin_to_int_vals[0]*bin_to_int_vals['level']
    bin_to_int_vals = bin_to_int_vals.groupby('intName')['final_level'].sum().reset_index()
            
    original_int_vals = gdx_data['x'].records
    original_int_vals = original_int_vals[original_int_vals['j'].isin(bin_to_int_vals['intName'])].copy(deep=True)
    original_int_vals = pd.merge(original_int_vals, bin_to_int_vals, left_on="j", right_on="intName", how="left")
    original_int_vals.drop(['level', 'intName', 'marginal','lower','upper','scale'], axis=1, inplace=True)
    original_int_vals.rename({'final_level': 'level'}, axis=1, inplace=True)
    
    if vars_with_lower_bounds:
        original_int_vals.loc[original_int_vals["j"].isin(vars_with_lower_bounds.keys()), 'level'] += list(vars_with_lower_bounds.values())

    sample.drop(sample[sample.j.isin(binName_list)].index, inplace=True)
    sample = pd.concat([sample, original_int_vals], ignore_index=True)

def map_vars(gdx_data, solution, obj_val):

    """
    helper function to map the solution returned from the qpu to the original problem and set the respective gams symbols
    """

    oldvars = gdx_data['x'].records
    oldvars.drop(['level'], inplace=True, axis=1)
   
    res = solution.merge(oldvars, how='right', on='j')
    vardict = gdx_data['j'].records
    separate_sym_domain = vardict['element_text'].str.split('(', expand=True)

    if len(separate_sym_domain.columns) == 1: # check if all variables are flat, i.e., no domain
        vardict['symbol'] = separate_sym_domain[0]
        vardict['domain'] = None
    else:
        vardict[['symbol','domain']] = separate_sym_domain[[0,1]]
        vardict['domain'] = vardict['domain'].str.rstrip(')')
        vardict['domain'] = vardict['domain'].str.strip(r"\'")
    
    vardict.drop(columns=['element_text'], inplace=True)
    vardict.rename({"uni": "j"}, axis=1, inplace=True)
    
    final = res.merge(vardict, how='right', on='j')
    final.loc[0, ['level', 'marginal']].to_list()
    final.loc[final['symbol']=='%obj%', 'level'] = obj_val
    
    for symbol in final['symbol'].unique():
        new_vals = []
        temp = final[final['symbol']==symbol].reset_index(drop=True)
        temp = temp[['domain','level', 'marginal', 'lower', 'upper', 'scale']]
        for i in range(len(temp)):
            domain = temp.loc[i, 'domain'].split(",") if type(temp.loc[i, 'domain']) == str else []
            new_vals += [(*domain, *temp.loc[i, ['level', 'marginal', 'lower', 'upper', 'scale']].to_list())]
        gams.set(symbol, new_vals, mergeType=MergeType.REPLACE)

if is_max:
    best_val = (sum_fixed_obj_var_coeffs + sum_lower_bound_of_int_vars) - best_val
else:
    best_val += (sum_fixed_obj_var_coeffs + sum_lower_bound_of_int_vars)

map_vars(gdx_data, sample, best_val)
endEmbeddedCode
$endIf.method

$ifThenE.run_examiner %check_examiner%

option solver=examiner;
%modelName%.optfile = 1;
solve %modelName% %direction% %obj% use %modelType%;

$endIf.run_examiner

$onListing