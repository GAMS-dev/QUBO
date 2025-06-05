import sys
import gamspy as gp
from qubo import Qubo

m = gp.Container(working_directory="./workdir")

f = gp.Set(m, name="f", records=["F1", "F2", "F3"], description="set of flights")
g = gp.Set(m, name="g", records=["G1", "G2"], description="set of gates")

i = gp.Alias(m, name="i", alias_with=f)
j = gp.Alias(m, name="j", alias_with=f)
k = gp.Alias(m, name="k", alias_with=g)
l = gp.Alias(m, name="l", alias_with=g)

buffer_time = gp.Parameter(
    m,
    name="buffer_time",
    records=10,
    description="the buffer time between two flights at the same gate",
)

arr_time = gp.Parameter(
    m,
    name="arr_time",
    domain=g,
    records=[["G1", 2], ["G2", 2.5]],
    description="the time it takes for a passenger to get from gate G to baggage claim",
)
dep_time = gp.Parameter(
    m,
    name="dep_time",
    domain=g,
    records=[["G1", 3], ["G2", 2]],
    description="the time it takes for a passenger to get from check-in to gate G",
)
trnsfr_time = gp.Parameter(
    m,
    name="trnsfr_time",
    domain=[k, l],
    records=[["G1", "G1", 5], ["G1", "G2", 7], ["G2", "G1", 9], ["G2", "G2", 5]],
    description="the time it takes to get from gate l to gate m",
)
passengers_in = gp.Parameter(
    m,
    name="passengers_in",
    domain=f,
    records=[["F1", 50], ["F3", 100]],
    description="the number of passengers from flight F which arrive at the airport",
)
passengers_out = gp.Parameter(
    m,
    name="passengers_out",
    domain=f,
    records=[["F1", 10], ["F2", 25], ["F3", 65]],
    description="'the number of passengers which depart from the airport on flight F",
)
passenger_trnsfr = gp.Parameter(
    m,
    name="passenger_trnsfr",
    domain=[i, j],
    records=[
        ["F1", "F1", 10],
        ["F1", "F3", 20],
        ["F2", "F1", 20],
        ["F2", "F2", 10],
        ["F3", "F1", 20],
        ["F3", "F3", 10],
    ],
    description="the number of lay over passengers which arrive on flight i and depart with flight j",
)
time_in = gp.Parameter(
    m,
    name="time_in",
    domain=f,
    records=[["F2", 20], ["F3", 35]],
    description="the arrival time of flight F",
)
time_out = gp.Parameter(
    m,
    name="time_out",
    domain=f,
    records=[["F1", 16], ["F2", 30], ["F3", 50]],
    description="the departure time of flight F",
)

x = gp.Variable(
    m,
    name="x",
    type="binary",
    domain=[f, g],
    description="1 iff flight F is assigned to gate G, 0 otherwise",
)

total_cost = gp.Variable(m, name="total_cost")

cost_fn = gp.Equation(m, name="cost_fn", description="Objective Function")
eq_use_one_gate = gp.Equation(
    m, name="eq_use_one_gate", domain=i, description="Allow only one gate per flight"
)
eq_restrict_arrival = gp.Equation(
    m,
    name="eq_restrict_arrival",
    domain=[i, j, k],
    description="NO arrival before departure at the same gate",
)
eq_restrict_arrival_linear = gp.Equation(
    m,
    name="eq_restrict_arrival_linear",
    domain=[i, j, k],
    description="linear alternative for its quadratic counterpart",
)

p = gp.Set(m, name="p", domain=[i, j])

p[i, j].where[
    (time_in[i] < time_in[j]) & (time_in[j] < time_out[i] + buffer_time)
] = True

cost_fn[...] = total_cost == gp.Sum(
    gp.Domain(i, k),
    (passengers_out[i] * dep_time[k] + passengers_in[i] * arr_time[k]) * x[i, k],
) + gp.Sum(
    gp.Domain(i, j, k, l),
    passenger_trnsfr[i, j] * trnsfr_time[k, l] * x[i, k] * x[j, l],
)

eq_use_one_gate[i] = gp.Sum(k, x[i,k]) == 1
eq_restrict_arrival[i,j,k].where[p[i,j]] = x[i,k]*x[j,k] == 0
eq_restrict_arrival_linear[i, j, k].where[p[i, j]] = x[i,k] + x[j,k] <= 1


fga = gp.Model(
    m,
    name="fga",
    problem="MIQCP",
    equations=[
        eq for eq in m.getEquations() if eq.name != "eq_restrict_arrival_linear"
    ],
    sense=gp.Sense.MIN,
    objective=total_cost,
)

fgal = gp.Model(
    m,
    name="fga",
    problem="MIQCP",
    equations=[eq for eq in m.getEquations() if eq.name != "eq_restrict_arrival"],
    sense=gp.Sense.MIN,
    objective=total_cost,
)

# fgal.solve(solver="CPLEX", output=sys.stdout)
# print(f"{fgal.objective_value = }")

# set penalty 650 !!This comes from the Scalar - pen_one_gate
q = Qubo(fgal, penalty=650)

q.solve(solver="CPLEX", options=gp.Options(time_limit=60, threads=4))

print(f"Original Objective Variable:\n{fgal._objective_variable.records}")
print(f"Variable x:\n{x.records}")
print(f"Variable total_cost:\n{total_cost.records}")
