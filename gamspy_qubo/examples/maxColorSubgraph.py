import sys
import gamspy as gp
from gamspy_qubo import Qubo

m = gp.Container(working_directory="./workdir")

nodes = gp.Set(m, name="nodes", records=["a", "b", "c", "d"])

num_clr = gp.Set(m, name="num_clr", records=[0, 1, 2])

n = gp.Alias(m, name="n", alias_with=nodes)
n1 = gp.Alias(m, name="n1", alias_with=nodes)
n2 = gp.Alias(m, name="n2", alias_with=nodes)

c = gp.Alias(m, name="c", alias_with=num_clr)
c1 = gp.Alias(m, name="c", alias_with=num_clr)
c2 = gp.Alias(m, name="c", alias_with=num_clr)

edges = gp.Set(
    m,
    name="edges",
    domain=[n1, n2],
    records=[["a", "b"], ["b", "c"], ["b", "d"], ["c", "d"]],
)

x = gp.Variable(
    m,
    name="x",
    type="binary",
    domain=[n, c],
    description="1 if node n is colored with color c and 0 otherwise",
)

total_cost = gp.Variable(m, name="total_cost")

cost_fn = gp.Equation(m, name="cost_fn", description="Objective function")

eq_get_one_clr = gp.Equation(
    m, name="eq_get_one_clr", domain=n, description="Each node gets only one color"
)

cost_fn[...] = total_cost == gp.Sum(gp.Domain(edges[n1, n2], c), x[n1, c] * x[n2, c])

eq_get_one_clr[n] = gp.Sum(c, x[n, c]) == 1

mcs = gp.Model(
    m,
    name="mcs",
    problem=gp.Problem.MIQCP,
    equations=m.getEquations(),
    sense=gp.Sense.MIN,
    objective=total_cost,
)

# mcs.solve(solver="CPLEX", output=sys.stdout)
# print(f"{mcs.objective_value = }")

q = Qubo(mcs, penalty=10)

q.solve(solver="CPLEX", options=gp.Options(time_limit=60))

print(f"Original Objective Variable:\n{mcs._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable total_cost:\n{total_cost.records}")
