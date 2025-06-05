import sys
import gamspy as gp
from qubo import Qubo

m = gp.Container(working_directory="./workdir")

i = gp.Set(m, name="i", records=[1, 2, 3, 4])

j = gp.Alias(m, name="j", alias_with=i)

uc = gp.Parameter(
    m,
    name="uc",
    domain=[i, j],
    records=[
        [1, 1, 2],
        [2, 2, 5],
        [3, 3, 2],
        [4, 4, 4],
        [1, 2, 8],
        [1, 3, 6],
        [1, 4, 10],
        [2, 3, 2],
        [2, 4, 6],
        [3, 4, 4],
    ],
)

pc = gp.Parameter(m, name="pc", domain=[i], records=[[1, 8], [2, 6], [3, 5], [4, 3]])

x = gp.Variable(m, name="x", type="binary", domain=i)
z = gp.Variable(m, name="z")

obj = gp.Equation(m, name="obj")
c1 = gp.Equation(m, name="c1")

obj[...] = gp.Sum(gp.Domain(i, j), x[i] * uc[i, j] * x[j]) == z
c1[...] = gp.Sum(i, pc[i] * x[i]) <= 16

qkp = gp.Model(
    m,
    name="qkp",
    equations=m.getEquations(),
    problem="MIQCP",
    sense=gp.Sense.MAX,
    objective=z,
)

q = Qubo(qkp, penalty=10)

q.solve()

print(f"Original Objective Variable:\n{qkp._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable z:\n{z.records}")
