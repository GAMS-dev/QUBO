import os

import gamspy as gp
import pandas as pd

from gamspy_qubo import Qubo

m = gp.Container()

i = gp.Set(m, name="i", records=[f"b{i}" for i in range(1, 7)])

c_data = pd.DataFrame(
    [
        ("b1", 3),
        ("b2", 2),
        ("b3", 1),
        ("b4", 1),
        ("b5", 3),
        ("b6", 2),
    ]
)

c = gp.Parameter(m, name="c", domain=i, records=c_data)

x = gp.Variable(m, name="x", type="binary", domain=i)

z = gp.Variable(m, name="z")

obj = gp.Equation(m, name="obj")
c1 = gp.Equation(m, name="c1")
c2 = gp.Equation(m, name="c2")
c3 = gp.Equation(m, name="c3")
c4 = gp.Equation(m, name="c4")

obj[...] = gp.Sum(i, c[i] * x[i]) == z
c1[...] = x["b1"] + x["b3"] + x["b6"] == 1
c2[...] = x["b2"] + x["b3"] + x["b5"] + x["b6"] == 1
c3[...] = x["b3"] + x["b4"] + x["b5"] == 1
c4[...] = x["b1"] + x["b2"] + x["b4"] + x["b6"] == 1

setPartition = gp.Model(
    m,
    name="setPartition",
    problem="MIP",
    equations=m.getEquations(),
    sense=gp.Sense.MIN,
    objective=z,
)

q = Qubo(setPartition, penalty=10, backend="kipu")
ACCESS_KEY_ID = os.environ["PLANQ_APP_ACCESS_KEY_ID"]
SECRET_ACCESS_KEY = os.environ["PLANQ_APP_SECRET_ACCESS_KEY"]

### shape = (6,6)
q.solve(secret_access_key=SECRET_ACCESS_KEY, access_key_id=ACCESS_KEY_ID)

print(f"Original Objective Variable:\n{setPartition._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable z:\n{z.records}")
