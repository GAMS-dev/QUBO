# QUBOs

A QUBO, Quadratic Unconstrained Binary Optimization, is a type of optimization problem where the aim is to find the best combination of binary choices (0/1) to maximize or minimize a quadratic objective function. It involves no constraints except that the variables are binary.

The QUBO model's significance in combinatorial optimization is heightened by its equivalence to the Ising model, which is prominent in physics. Consequently, the broad range of optimization problems solved effectively by state-ofthe-art QUBO solution methods are joined by an important domain of problems arising in physics applications.

## QUBO Reformulation

The reformulation is derived from, A Tutorial on Formulating and Using QUBO Models. Glover, F., Kochenberger, G., & Du, Y. [2019](https://arxiv.org/abs/1811.11538). Following are some examples included to test the qubo reformulation.

1. setPacking.gms, Set Packing Problem (max)
2. O1program.gms, General 0/1 Problem (max)
3. QAP.gms, Quadratic Assignment Problem (min)
4. setPartition.gms, Set Partitioning Problem (min)
5. QKP.gms, Quadratic Knapsack Problem (max)
6. generalIP.gms, a general integer problem (max)
7. qplib_5881.gms, a flat/scalar gms file (max)
8. knights.gms, A Max problem from GAMS modlib


## Required Packages

1. gamsapi[transfer], [link](https://www.gams.com/latest/docs/API_PY_GETTING_STARTED.html#PY_PIP_INSTALL_BDIST)
2. dwave-system, required when solving on [Dwave's](https://docs.ocean.dwavesys.com/projects/system/en/latest/installation.html) Hybrid QPU otherwise optional.

## Input

Once the problem is defined, it can be solved through the qubo reformulation by including the `qubo_solve.gms` using $batinclude.
The `qubo_solve.gms` file should be in the same location where the main gms file is located. If not, the location of the file must be specified by either including it in the $batinclude statement, for e.g., `$batinclude 'location\of\the\file\qubo_solve.gms'` or by setting the command line parameter, `-IDIR`.

The `qubo_solve.gms` requires the following 5 positional arguments. Since they are positional arguments they must be in the exact order as mentioned below.

1. modelName
2. modelType
3. objective (max/min)
4. objectiveVariable
5. Penalty factor for the constraints

Following is the list of optional `-key=val` pair arguments, some of which are method specific.

6. method, [qpu, classic] (default: classic)
7. solver, choice of miqcp solver (default: cplex | effective only if `-method=classic`).
8. maxIter, Number of times the problem is solved on the QPU (default: 1 | effective only if `-method=qpu`)
9. timeLimit, Time limit for 1 iteration on QPU or TimeLimit for a classical solve (default: 10)
10. numThreads, Number of threads to be used in case of a classical solve (default: min(8,num_of_cores) | effective only if `-method=classic`)
11. logOn, Creates a log for the reformulation [0, 1, 2] (default: 0, don't create a log)
12. examinerOn, [0, 1] (default: 0) The quality of returned qubo solution w.r.t the original problem can be checked through the use of `examiner` [tool](https://www.gams.com/latest/docs/S_EXAMINER.html).

Note: Generating the API key and setting up the Python-Dwave Environment is considered to be available when chosen method of solving is `qpu`.

## How to run

- Download GAMS from https://www.gams.com/download/
- Install GAMS
- Run the main gms file with the desired options by including them in the main file through the `$batinclude` statement. For e.g., `$batinclude qubo_solve.gms setPacking MIP max z 6 -solver=cplex -timeLimit=60 -numThreads=2 -logOn=2`
  - from GAMS Studio: Open the main problem file in GAMS studio. If qubo_solve.gms is not in the same directory as the main problem file, enter `-IDIR=<path//to//qubo_solve.gms>` in the [parameter editor](https://www.gams.com/latest/docs/T_STUDIO.html#STUDIO_TOOLBAR) and hit the run button (or press F9)
  - from the command line
    ```
    gams '.\QAP.gms' -IDIR=<path//to//qubo_solve.gms>
    ```

## Output

The script generates two gdx files. One for the standard problem which is saved as `modelName.gdx` and another for the reformulted model, saved as `qout_modeName.gdx`. A successful run will then return the level of binary variables and the objective variable.

## Limitations

Note: In order to binarize the right hand side of each constraint, the reformulation needs to generate slack variables depending on the RHS. The number of binary variables to be generated for a number $x$ is $2^n$ where n = $log_2(x)$. This with the fact that there can be many constraints, can make the reformulation computationally challenging.

1. The reformulation does not work with continuous variables with the exception of the free objective variable.
2. When binarizing integer variables they must have an upper bound less than or equal to 1e4.
3. The coefficients for the variables in the constraints must be integers.
4. The use of quadratic terms is limited to binary variables where levels are not equal to 1.
5. Although the objective function can have quadratic terms, the constraints cannot have any quadratic terms. This is due to the penalization step. The step requires the constraint to be squared which results in a polynomial of degree greater than 2. This becomes a problem with 3 or more interacting variables. For example, $(xy + yz)^2 = (xy)^2 + (yz)^2 + 2xy^2z$. Here, the term $2xy^2z$ is problematic since it cannot be reduced to bilinear terms.
6. The reformulation expects the objective function to be defined via the use of a scalar equation, i.e., the symbol `iobj` must be nonzero in the gdx file created by [Convert](https://www.gams.com/latest/docs/S_CONVERT.html).

The reformulation will throw appropriate exceptions if the limitations are not statisfied.

## Testing

The file `test_qubo_solve.gms` tests the correctness of qubo_solve in certain scenarios. The test file should be in the same location as `qubo_solve.gms`. There is a test for checking the correctness of the reformulation and solution obtained from the Dwave QPU. This is not enabled by default. In order to enable this test one should run the file with the command line option `--TESTDWAVE=yes`. It follows that the required python packages are already present in the python environment defined by `GMSPYTHONLIB`.

## Choosing the right Penalty

A penalty value that is too large can impede the solution process as the penalty terms overwhelm the original objective function information, making it difficult to distinguish the quality of one solution from another. On the other hand, a penalty value that is too small jeopardizes the search for feasible solutions. Generally, there is a ‘Goldilocks region’ of considerable size that contains penalty values that work well. A little preliminary thought about the model can yield a ballpark estimate of the original objective function value. Taking P to be some percentage (75% to 150%) of this estimate is often a good place to start. In the end, solutions generated can always be checked for feasibility, leading to changes in penalties and further rounds of the solution process as needed to zero in on an acceptable solution.