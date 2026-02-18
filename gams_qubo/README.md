# QUBO Reformulation in GAMS

This is the GAMS implementation for the QUBO reformulation tool. Following are some examples included to test the QUBO reformulation.

1. [setPacking.gms](./examples/setPacking.gms), Set Packing Problem (max)
2. [O1program.gms](./examples/01program.gms), General 0/1 Problem (max)
3. [QAP.gms](./examples/QAP.gms), Quadratic Assignment Problem (min)
4. [setPartition.gms](./examples/setPartition.gms), Set Partitioning Problem (min)
5. [QKP.gms](./examples/QKP.gms), Quadratic Knapsack Problem (max)
6. [generalIP.gms](./examples/generalIP.gms), a general integer problem (max)
7. [qplib_5881.gms](./examples/qplib_5881.gms), a flat/scalar gms file (max)
8. [knights.gms](./examples/knights.gms), A Max problem from GAMS modlib
9. [Q01.gms](./examples/Q01.gms), Quadratic 0/1 Problem (min)
10. [flightGate.gms](./examples/flightGate.gms), Flight gate assignment Problem (min)
11. [maxColorSubgraphs.gms](./examples/maxColorSubgraph.gms), Maximum colorable subgraph problem (min)
12. [tsp.gms](./examples/tsp.gms), Traveling Salesman Problem (min)


## Required Packages

1. gamsapi[transfer], [link](https://www.gams.com/latest/docs/API_PY_GETTING_STARTED.html#PY_PIP_INSTALL_BDIST)
2. dwave-system, required when solving on [Dwave's](https://docs.ocean.dwavesys.com/projects/system/en/latest/installation.html) Hybrid QPU otherwise optional.

## Input

Once the problem is defined, it can be solved through the QUBO reformulation by including the `qubo_solve.gms` using $batinclude.
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
12. examinerOn, [0, 1] (default: 0) The quality of returned QUBO solution w.r.t the original problem can be checked through the use of `examiner` [tool](https://www.gams.com/latest/docs/S_EXAMINER.html).
13. getQ, [y, n] (default: n) Specify if only the Q-matrix is exported to a CSV file. This will produce a CSV file with the following filename schema, `modelName_p(penalty)_c(total_offset).csv`

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

The script generates two gdx files. One for the standard problem which is saved as `modelName.gdx` and another for the reformulated model, saved as `qout_modeName.gdx`. A successful run will then return the level of binary variables and the objective variable.

## Testing

The file `test_qubo_solve.gms` tests the correctness of qubo_solve in certain scenarios. The test file should be in the same location as `qubo_solve.gms`. There is a test for checking the correctness of the reformulation and solution obtained from the Dwave QPU. This is not enabled by default. In order to enable this test one should run the file with the command line option `--TESTDWAVE=yes`. It follows that the required python packages are already present in the python environment defined by `GMSPYTHONLIB`.


## Acknowledgments
The QUBO reformulation tool has been developed under the financial support of:

- ProvideQ (BMWi project, ID: 01MQ22006D)
- QuSol (BMFTR project, ID: 13N17172)