import sys
import gamspy as gp
import random as rp
from qubo import Qubo

rp.seed(42)

m = gp.Container(working_directory="./workdir")


i = gp.Set(m, name="i", domain=["*"], records=[f"i{i}" for i in range(1, 6)])


def generate_parameters(
    m: gp.Container, num: int, prefix: str, domain: gp.Set
) -> gp.Container:
    for i in range(1, num + 1):
        name = f"{prefix}{i}"
        gp.Parameter(
            m,
            name=name,
            domain=domain,
            records=[
                (ele, rp.randint(1, 20))
                for ele, _ in domain.records.itertuples(index=False)
            ],
        )

    return m


a1, a2, a3, a4 = generate_parameters(m, num=4, prefix="a", domain=i).getParameters()

x = gp.Variable(m, name="x", type="binary", domain=i)

obj = gp.Sum(i, a1[i] * x[i])

gp.Equation(m, name="e1", definition=gp.Sum(i, a2[i] * x[i]) <= 25)
gp.Equation(m, name="e2", definition=gp.Sum(i, a3[i] * x[i]) == 21)
gp.Equation(m, name="e3", definition=gp.Sum(i, a4[i] * x[i]) >= 10)

zero_one = gp.Model(
    m,
    name="zero_one",
    problem="MIP",
    equations=m.getEquations(),
    sense=gp.Sense.MAX,
    objective=obj,
)

# zero_one.solve(solver="CPLEX", output=sys.stdout)
# print(f"{zero_one.objective_value = }")

q = Qubo(zero_one, penalty=10)

# qd, qi, qconst = q.transform(penalty=10)
# h, J, const = q.qubo_to_ising(qd.toDict())

# q.solve(solver="CPLEX", options=gp.Options(equation_listing_limit=1))
q.solve(solver="CPLEX")
q.map_solution()

print(f"Original Objective Variable:\n{zero_one._objective_variable.records}")
print(f"Variable x:\n{x.records}")

# q_mat = q.get_q_matrix()
# print(q_mat.shape)
