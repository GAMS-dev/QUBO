import gamspy as gp

from gamspy_qubo import Qubo

m = gp.Container()

i = gp.Set(m, name="i", records=[f"b{i}" for i in range(1, 5)])

cost = gp.Parameter(
    m, name="cost", domain=i, records=[["b1", 2], ["b2", 4], ["b3", 7], ["b4", 9]]
)

x = gp.Variable(m, name="x", type="integer", domain=i)

x.up[i] = 9
x.fx["b3"] = 5

y = gp.Variable(m, name="y", type="binary", domain=i)

z = gp.Variable(m, name="z")

newX = gp.Variable(m, name="newX", type="integer", domain=i)

newX.up[i] = 5

obj = gp.Equation(m, name="obj")
c1 = gp.Equation(m, name="c1")
c2 = gp.Equation(m, name="c2")
c3 = gp.Equation(m, name="c3")

obj[...] = (
    gp.Sum(i, cost[i] * x[i]) + gp.Sum(i, cost[i] * y[i]) - gp.Sum(i, newX[i]) == z
)

c1[...] = gp.Sum(i, x[i]) <= 20
c2[...] = gp.Sum(i, y[i]) <= 3
c3[...] = gp.Sum(i, newX[i]) >= 3

demo_model = gp.Model(
    m,
    name="demo_model",
    problem="MIP",
    equations=m.getEquations(),
    sense=gp.Sense.MAX,
    objective=z,
)

q = Qubo(demo_model, penalty=10, backend="dwave")
q.solve(num_reads=1000)

print(f"Original Objective Variable:\n{demo_model._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable z:\n{z.records}")
