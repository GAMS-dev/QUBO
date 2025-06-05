"""
This MIP model finds the maximum number of knights that can be
placed on a board. Two different formulations are presented.
The second formulation is 'tight' and may perform better with certain
MIP codes. Once we found the max number of knights, we solve a series
of MIPs to find ALL solutions.

We will use lags (relative positions) to describe the allowed moves.
The labels H and V indicate horizontal and vertical moves as shown
below:

                 0 0
                0   0
                  X
                0   0
                 0 0


Dudeney, H E, Amusements in Mathematics. Dover, New York, 1970.

Keywords: mixed integer linear programming, maximum knights problem, mathematics
"""

import sys
import gamspy as gp
import pandas as pd
from qubo import Qubo

m = gp.Container(working_directory="./workdir")

i = gp.Set(
    m,
    name="i",
    records=[f"{i}" for i in range(1, 9)],
    description="size of board",
)

n = gp.Set(
    m,
    name="n",
    records=[f"m{i}" for i in range(1, 9)],
    description="number of possible moves",
)

j = gp.Alias(m, name="j", alias_with=i)
k = gp.Alias(m, name="k", alias_with=i)


moves_data = pd.DataFrame(
    [
        ("H", "m1", -2),
        ("H", "m2", -2),
        ("H", "m3", -1),
        ("H", "m4", -1),
        ("H", "m5", 1),
        ("H", "m6", 1),
        ("H", "m7", 2),
        ("H", "m8", 2),
        ("V", "m1", -1),
        ("V", "m2", 1),
        ("V", "m3", -2),
        ("V", "m4", 2),
        ("V", "m5", -2),
        ("V", "m6", 2),
        ("V", "m7", -1),
        ("V", "m8", 1),
    ]
)

move = gp.Parameter(
    m,
    name="move",
    domain=["*", n],
    records=moves_data,
    description="all possible knight moves",
)

total = gp.Variable(m, name="total")
x = gp.Variable(m, name="x", type="binary", domain=[i, j])

deftotal = gp.Equation(m, name="deftotal", description="total knights on board")
defmove = gp.Equation(m, name="defmove", domain=[i, j], description="move restrictions")
defmovex = gp.Equation(
    m, name="defmovex", domain=[n, i, j], description="move restrictions"
)

deftotal[...] = total == gp.Sum(gp.Domain(i, j), x[i, j])
defmove[i, j] = gp.Sum(n, x[i + move["H", n], j + move["V", n]]) <= gp.Card(i) * (
    1 - x[i, j]
)
defmovex[n, i, j] = x[i + move["H", n], j + move["V", n]] <= 1 - x[i, j]

knight = gp.Model(
    m,
    name="knight",
    problem="MIP",
    equations=[deftotal, defmove],
    sense=gp.Sense.MAX,
    objective=total,
)

knightx = gp.Model(
    m,
    name="knightx",
    problem="MIP",
    equations=[deftotal, defmovex],
    sense=gp.Sense.MAX,
    objective=total,
)


knight.solve(
    solver="CPLEX",
    output=sys.stdout,
    options=gp.Options(
        relative_optimality_gap=0, absolute_optimality_gap=0.999, time_limit=60
    ),
)

q = Qubo(knight, penalty=10)

q.solve(
    options=gp.Options(
        relative_optimality_gap=0, absolute_optimality_gap=0.999, time_limit=60
    )
)

print(f"Original Objective Variable:\n{knight._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable total:\n{total.records}")
