import sys
import gamspy as gp
import random as rp
import pandas as pd
from qubo import Qubo


m = gp.Container(working_directory="./workdir")

facility = gp.Set(
    m, name="facility", domain=["*"], records=["chicago", "boston", "denver"]
)


i = gp.Alias(m, name="i", alias_with=facility)

j = gp.Alias(m, name="j", alias_with=facility)

k = gp.Alias(m, name="k", alias_with=facility)

l = gp.Alias(m, name="l", alias_with=facility)

flow = gp.Parameter(
    m,
    name="flow",
    domain=[i, j],
)

dist = gp.Parameter(
    m,
    name="dist",
    domain=[k, l],
)

uflow = pd.DataFrame(
    [
        ("chicago", "boston", 5),
        ("chicago", "denver", 2),
        ("boston", "denver", 3),
    ],
    columns=["i", "j", "value"],
)

udist = pd.DataFrame(
    [
        ("chicago", "boston", 8),
        ("chicago", "denver", 15),
        ("boston", "denver", 13),
    ],
    columns=["k", "l", "value"],
)

_swap = uflow.rename(columns={"i": "j", "j": "i"})
_swap = pd.concat([uflow, _swap], ignore_index=True)
_swap = _swap.astype({'i': 'category', 'j': 'category', 'value': 'float'})
flow.records = _swap

_swap = udist.rename(columns={"k": "l", "l": "k"})
_swap = pd.concat([udist, _swap], ignore_index=True)
_swap = _swap.astype({'k': 'category', 'l': 'category', 'value': 'float'})
dist.records = _swap

x = gp.Variable(m, name="x", domain=[i, j], type="binary")

c1 = gp.Equation(
    m,
    name="c1",
    domain=j
)
c1[j] = gp.Sum(i, x[i,j]) == 1


c2 = gp.Equation(
    m,
    name="c2",
    domain=i
)
c2[i] = gp.Sum(j, x[i,j]) == 1

obj = gp.Sum(gp.Domain(i,j,k,l), flow[i,j]*x[i,k]*x[j,l]*dist[k,l])


qap = gp.Model(
    m,
    name="qap",
    problem="MIQCP",
    equations=m.getEquations(),
    sense=gp.Sense.MIN,
    objective=obj
)

# qap.solve(solver="CPLEX")

q = Qubo(qap, penalty=200)
print(q.solve(solver="CPLEX"))
print(f"\n{x.records}")
