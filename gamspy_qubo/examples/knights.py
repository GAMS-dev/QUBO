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
import numpy as np
from gamspy_qubo import Qubo

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
d = gp.Set(m, name="d", records=["H", "V"])

move = gp.Parameter(
    m,
    name="move",
    domain=[d, n],
    description="all possible knight moves",
)

move.setRecords(
    np.array(
        [
            # m1 m2 m3 m4 m5 m6 m7 m8
            [-2, -2, -1, -1, +1, +1, +2, +2],  # H
            [-1, +1, -2, +2, -2, +2, -1, +1],  # V
        ]
    )
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

print(x.pivot())

q = Qubo(knight, penalty=10)

q.solve(
    options=gp.Options(
        relative_optimality_gap=0, absolute_optimality_gap=0.999, time_limit=60
    )
)

print(f"Original Objective Variable:\n{knight._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable total:\n{total.records}")
