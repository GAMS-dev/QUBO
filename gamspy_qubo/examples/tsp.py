import sys
import gamspy as gp
from gamspy_qubo import Qubo

m = gp.Container(working_directory="./workdir")

nodes = gp.Set(m, name="nodes", records=["a", "b", "c", "d"])
position = gp.Set(m, name="position", records=[0, 1, 2, 3, 4])

n = gp.Alias(m, name="n", alias_with=nodes)
n1 = gp.Alias(m, name="n1", alias_with=nodes)
n2 = gp.Alias(m, name="n2", alias_with=nodes)

i = gp.Alias(m, name="i", alias_with=position)
i1 = gp.Alias(m, name="i1", alias_with=position)

edges = gp.Set(m, name="edges", domain=[n1, n2], description="possible edges")

weights = gp.Parameter(
    m,
    name="weights",
    domain=[n1, n2],
    records=[
        ["a", "b", 8],
        ["a", "c", 1],
        ["a", "d", 1],
        ["b", "a", 4],
        ["b", "c", 1],
        ["b", "d", 6],
        ["c", "a", 2],
        ["c", "b", 4],
        ["c", "d", 9],
        ["d", "a", 7],
        ["d", "b", 1],
        ["d", "c", 5],
    ],
)

edges[n1, n2] = weights[n1, n2]

total_cost = gp.Variable(m, name="total_cost")

x = gp.Variable(m, name="x", type="binary", domain=[n, i])

x.fx[n, i].where[n.first & i.first] = 1
x.fx[n, i].where[n.first & i.last] = 1
x.fx[n, i].where[n.first & ~i.first & ~i.last] = 0
x.fx[n, i].where[~n.first & i.first] = 0
x.fx[n, i].where[~n.first & i.last] = 0

cost_fn = gp.Equation(m, name="cost_fn", description="Objective function")
eq_exact_pos = gp.Equation(
    m, name="eq_exact_pos", domain=i, description="each position is used exactly once"
)
eq_exact_node = gp.Equation(
    m, name="eq_exact_node", domain=n, description="visit each node exactly once"
)

cost_fn[...] = total_cost == gp.Sum(
    gp.Domain(edges[n1, n2], i).where[~i.last],
    weights[n1, n2] * x[n1, i] * x[n2, i + 1],
)

eq_exact_pos[i].where[~i.first & ~i.last] = gp.Sum(n.where[~n.first], x[n, i]) == 1

eq_exact_node[n].where[~n.first] = gp.Sum(i.where[~i.first & ~i.last], x[n, i]) == 1

tsp = gp.Model(
    m,
    name="tsp",
    problem="MIQCP",
    equations=m.getEquations(),
    sense=gp.Sense.MIN,
    objective=total_cost,
)

# tsp.solve(
#     solver="CPLEX",
#     output=sys.stdout,
#     options=gp.Options(hold_fixed_variables=True),
# )

q = Qubo(tsp, penalty=10, options=gp.Options(hold_fixed_variables=True))

q.solve()

print(f"Original Objective Variable:\n{tsp._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable total_cost:\n{total_cost.records}")
