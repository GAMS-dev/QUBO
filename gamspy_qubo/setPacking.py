import sys
import gamspy as gp
from qubo import Qubo

m = gp.Container(working_directory="./workdir")

i = gp.Set(m, name="i", records=[f"{i}" for i in range(1, 5)])

x = gp.Variable(m, name="x", type="binary", domain=i)

z = gp.Variable(m, name="z")

obj = gp.Equation(m, name="obj")
c1 = gp.Equation(m, name="c1")
c2 = gp.Equation(m, name="c2")
c3 = gp.Equation(m, name="c3")

obj[...] = gp.Sum(i, x[i]) == z

c1[...] = gp.Sum(i.where[~i.sameAs("2")], x[i]) <= 1

c2[...] = gp.Sum(i.where[gp.Ord(i) <= 2], x[i]) <= 1

c3[...] = x["2"] + x["3"] == 1

setPacking = gp.Model(
    m,
    name="setPacking",
    problem="MIP",
    equations=m.getEquations(),
    sense=gp.Sense.MAX,
    objective=z,
)

# setPacking.solve(solver="CPLEX", output=sys.stdout)
# print(f"{setPacking.objective_value = }")

q = Qubo(setPacking, penalty=10)

q.solve(solver="CPLEX")

print(f"Original Objective Variable:\n{setPacking._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable z:\n{z.records}")
