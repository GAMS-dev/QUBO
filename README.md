# QUBOs

A QUBO, Quadratic Unconstrained Binary Optimization, is a type of optimization problem where the aim is to find the best combination of binary choices (0/1) to maximize or minimize a quadratic objective function. It involves no constraints except that the variables are binary.

The QUBO model's significance in combinatorial optimization is heightened by its equivalence to the Ising model, which is prominent in physics. Consequently, the broad range of optimization problems solved effectively by state-of-the-art QUBO solution methods are joined by an important domain of problems arising in physics applications.

**Note:** This repository contains two versions of the QUBO tool. The GAMS implementation is in the gams_qubo folder. The Python-based GAMSPy implementation is in the gamspy_qubo folder. The GAMSPy implementation will be the primary focus of all future updates and active maintenance.

## QUBO Reformulation

The reformulation is derived from, A Tutorial on Formulating and Using QUBO Models. Glover, F., Kochenberger, G., & Du, Y. [2019](https://arxiv.org/abs/1811.11538). 

## Quickstart

### 1. Clone the repository into your directory

```bash
git clone git@git.gams.com:devel/qubo.git
cd qubo/gamspy_qubo
```

### 2. Install and Project Setup

This project uses `uv` for environment management. You can install the base package or include specific dependency groups for development or quantum backend access.

**Base installation:**

```bash
uv sync
```

**Installing specific dependency groups:**
You can include optional groups using the `--group` flag:

* **dev**: Tools for testing and linting (mypy, ruff, pytest, etc.).
* **dwave**: Includes `dwave-ocean-sdk` for solving on D-Wave hardware.
* **kipu**: Includes `planqk-service-sdk` for Kipu Quantum integration.

```bash
# Example: Install with dev and dwave support
uv sync --group dev --group dwave
```
You can also install all dependencies at once with `uv sync --all-group`.

#### Fallback method (using `requirements.txt`) if you are not using uv:

**macOS / Linux:**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

**Windows:**

```bash
python -m venv .venv
.venv\Scripts\activate
```

Then, install using requirements.txt:

```bash
pip install -r requirements.txt
```

This will not install the dependencies for the quantum backends.

## Project Structure

```
gamspy_qubo/
├── src/
│   └── gamspy_qubo/
│       ├── __init__.py
│       └── qubo.py
├── examples/
│   └── ...
├── tests/
│   └── ...
├── pyproject.toml
├── requirements.txt
├── README.md
└── .gitignore
```

## Using the Package

Once installed, you can use it like:

```python
from gamspy_qubo import Qubo
```

Or run example scripts:

```bash
uv run examples/qap.py
```

## Examples
Following are some examples included to test the QUBO reformulation.

1. [setPacking.py](./examples/setPacking.py), Set Packing Problem (max)
2. [O1program.py](./examples/01program.py), General 0/1 Problem (max)
3. [QAP.py](./examples/QAP.py), Quadratic Assignment Problem (min)
4. [setPartition.py](./examples/setPartition.py), Set Partitioning Problem (min)
5. [QKP.py](./examples/QKP.py), Quadratic Knapsack Problem (max)
6. [generalIP.py](./examples/generalIP.py), a general integer problem (max)
7. [qplib_5881.py](./examples/qplib_5881.py), a flat/scalar py file (max)
8. [knights.py](./examples/knights.py), A Max problem from GAMS modlib
9. [Q01.py](./examples/Q01.py), Quadratic 0/1 Problem (min)
10. [flightGate.py](./examples/flightGate.py), Flight gate assignment Problem (min)
11. [maxColorSubgraphs.py](./examples/maxColorSubgraph.py), Maximum colorable subgraph problem (min)
12. [tsp.py](./examples/tsp.py), Traveling Salesman Problem (min)


## Input

Once your GAMSPy problem is defined and a `gp.Model` object is instantiated, you can solve it through the QUBO reformulation by wrapping your model in the `Qubo` class. 

```python
from gamspy_qubo import Qubo

# 1. Wrap your existing GAMSPy model
qubo_model = Qubo(
    model=demo_model,  # Your gp.Model instance
    name="QUBO",       # (Optional) Name for the reformulated model
    penalty=10,        # (Optional) Penalty factor for constraints (default: 1)
)

# 2. Solve the QUBO model
# By default, solver="cplex". You can also specify quantum backends.
q.solve(solver="dwave", num_reads=1000)

```
*(Note: Generating the API key and setting up the required python environments for `dwave` or `kipu` must be done prior to using them as a solver backend).*

### Keyword Arguments for `.solve()`

* **`solver`**: Choice of solver backend. Supported classical solvers are standard GAMSPy MIQCP solvers (default: `"cplex"`). Supported quantum backends sofar are `"dwave"` and `"kipu"`.
* **`**kwargs`**: Any additional backend-specific arguments. For example, if using `"dwave"`, you can pass arguments like `num_reads=1000`. If using a classic solver, arguments are passed natively to GAMSPy's solve method.


## Output

The GAMSPy integration automatically translates the QUBO solution back to your problem's original scope. Upon a successful run, your initial model variables are directly populated with the new levels. 

## Limitations

Note: In order to binarize the right hand side of each constraint, the reformulation needs to generate slack variables depending on the RHS. The number of binary variables to be generated for a number $x$ is $2^n$ where n = $log_2(x)$. This with the fact that there can be many constraints, can make the reformulation computationally challenging.

1. The reformulation does not work with continuous variables with the exception of the free objective variable.
2. When binarizing integer variables they must have an upper bound less than or equal to 1e4.
3. The coefficients for the variables in the constraints must be integers.
4. The use of quadratic terms is limited to binary variables where levels are not equal to 1.
5. Although the objective function can have quadratic terms, the constraints cannot have any quadratic terms. This is due to the penalization step. The step requires the constraint to be squared which results in a polynomial of degree greater than 2. This becomes a problem with 3 or more interacting variables. For example, $(xy + yz)^2 = (xy)^2 + (yz)^2 + 2xy^2z$. Here, the term $2xy^2z$ is problematic since it cannot be reduced to bilinear terms.
6. The reformulation expects the objective function to be defined via the use of a scalar equation, i.e., the symbol `iobj` must be nonzero in the gdx file created by [Convert](https://www.gams.com/latest/docs/S_CONVERT.html).

The reformulation will throw appropriate exceptions if the limitations are not satisfied.

## Choosing the right Penalty

A penalty value that is too large can impede the solution process as the penalty terms overwhelm the original objective function information, making it difficult to distinguish the quality of one solution from another. On the other hand, a penalty value that is too small jeopardizes the search for feasible solutions. Generally, there is a ‘Goldilocks region’ of considerable size that contains penalty values that work well. A little preliminary thought about the model can yield a ballpark estimate of the original objective function value. Taking P to be some percentage (75% to 150%) of this estimate is often a good place to start. In the end, solutions generated can always be checked for feasibility, leading to changes in penalties and further rounds of the solution process as needed to zero in on an acceptable solution.

## Run Tests

There are tests available that can be run using the `pytest` package. You can install it via,

```bash
uv sync --group dev
```
The tests can then be run via `pytest tests\test_qubo.py`

## GAMS Version

The original tool was developed for the GAMS modeling language. While this repository now focuses on the GAMSPy (Python) implementation, the classic GAMS version remains available in the GAMS subfolder.

## Acknowledgments
The QUBO reformulation tool has been developed under the financial support of:

- ProvideQ (BMWi project, ID: 01MQ22006D)
- QuSol (BMFTR project, ID: 13N17172)