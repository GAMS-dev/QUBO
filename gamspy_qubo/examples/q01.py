"""This example handles case where there are quadratic constraints and
objective function.

It should be noted that the quadratic constraints do not contain
pair-wise quadratic terms. This is a limitation of QUBO reformulation.
The quadratic terms otherwise, for e.g., (x_1)^2 can be treated as x_1
since x_1 is binary.
"""

import sys
import gamspy as gp
import pandas as pd
from gamspy_qubo.src.gamspy_qubo.qubo import Qubo

m = gp.Container(working_directory="./workdir")

i = gp.Set(m, name="i", records=[1, 2, 3, 4, 5])

j = gp.Alias(m, name="j", alias_with=i)

c = gp.Parameter(
    m, name="c", domain=i, records=[[1, 6], [2, 4], [3, 8], [4, 5], [5, 5]]
)
a1 = gp.Parameter(
    m, name="a1", domain=i, records=[[1, 2], [2, 2], [3, 4], [4, 4], [5, 2]]
)
a2 = gp.Parameter(
    m, name="a2", domain=i, records=[[1, 1], [2, 2], [3, 2], [4, 1], [5, 2]]
)
a3 = gp.Parameter(
    m, name="a3", domain=i, records=[[1, 3], [2, 3], [3, 3], [4, 4], [5, 4]]
)

x = gp.Variable(m, name="x", type="binary", domain=i)

z = gp.Variable(m, name="z")

obj = gp.Equation(m, name="obj")
c1 = gp.Equation(m, name="c1")
c2 = gp.Equation(m, name="c2")
c3 = gp.Equation(m, name="c3")

obj[...] = gp.Sum(i, x[i] * c[i] * x[i]) == z
c1[...] = gp.Sum(i, a1[i] * x[i]) <= 12
c2[...] = (
    gp.Sum(i.where[gp.Ord(i) <= 3], a2[i] * x[i])
    - gp.Sum(i.where[gp.Ord(i) > 3], a2[i] * x[i])
    >= 5
)
c3[...] = gp.Sum(i.where[gp.Ord(i) <= 3], a3[i] * x[i]) <= 10

quadZeroOne = gp.Model(
    m,
    name="quadZeroOne",
    problem=gp.Problem.MIQCP,
    equations=m.getEquations(),
    sense=gp.Sense.MIN,
    objective=z,
)

q = Qubo(quadZeroOne, penalty=10)

q.solve()

print(f"Original Objective Variable:\n{quadZeroOne._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable z:\n{z.records}")
