from __future__ import annotations

import numpy as np
import pytest

import gamspy as gp
from gamspy import Container, ModelStatus
from gamspy_qubo import Qubo

from gamspy.exceptions import ValidationError


@pytest.fixture
def data():
    m = Container()
    i = gp.Set(m, "i", "*", records=range(1, 6))
    x = gp.Variable(m, "x", domain=[i], type="binary")
    c1 = gp.Equation(m, "c1")
    obj = gp.Equation(m, "obj")
    z = gp.Variable(m, "objective")

    yield m, i, x, c1, obj, z
    m.close()


def test_qubo_bad_init():
    pytest.raises(ValidationError, Qubo, "not_gp_Model_instance", 1)


def test_qubo_continous_variable(data):
    m, _, _, c1, obj, z = data

    real_x = gp.Variable(m, "real_X")
    c1[...] = real_x == 1

    obj[...] = real_x == z

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)

    with pytest.raises(ValidationError) as e:
        qubo.transform()
    assert "There are continuous variables. Quitting." in str(e.value)


def test_qubo_real_coefficients(data):
    m, _, _, c1, obj, z = data

    bin_x = gp.Variable(m, "real_X", type="binary")
    c1[...] = 2.3 * bin_x == 1
    obj[...] = bin_x == z

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)

    with pytest.raises(ValidationError) as e:
        qubo.transform()
    assert "Reformulation with Non-Integer Coefficients not possible. Quitting." in str(
        e.value
    )


def test_qubo_nonscalar_objective(data):
    m, i, x, c1, _, z = data

    ind_obj = gp.Equation(m, "ind_eqn", domain=[i])
    ind_obj[...] = x[...] == z

    c1[...] = gp.Sum(i, x) == 1

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, ind_obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)

    with pytest.raises(ValidationError) as e:
        qubo.transform()
    assert (
        "The objective is not defined using a scalar equation. `iobj` in gdx is empty. Quitting."
        in str(e.value)
    )


def test_qubo_nonscalar_objective(data):
    m, i, x, c1, _, z = data

    ind_obj = gp.Equation(m, "ind_eqn", domain=[i])
    ind_obj[...] = x[...] == z

    c1[...] = gp.Sum(i, x) == 1

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, ind_obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)

    with pytest.raises(ValidationError) as e:
        qubo.transform()
    assert (
        "The objective is not defined using a scalar equation. `iobj` in gdx is empty. Quitting."
        in str(e.value)
    )


def test_qubo_infeasible_constraint(data):
    m, i, x, c1, obj, z = data

    obj[...] = gp.Sum(i, x) == z
    c1[...] = gp.Sum(i, x) == 10

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)
    with pytest.raises(ValidationError) as e1:
        qubo.transform()
    assert (
        "Constraint is infeasible: e1"
        in str(e1.value)
    )

    c1[...] = gp.Sum(i, x) <= -10

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)
    with pytest.raises(ValidationError) as e2:
        qubo.transform()
    assert (
        "Constraint is infeasible: e1"
        in str(e2.value)
    )

    c1[...] = gp.Sum(i, x) >= 10

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)
    with pytest.raises(ValidationError) as e3:
        qubo.transform()
    assert (
        "Constraint is infeasible: e1"
        in str(e3.value)
    )

def test_qubo_upper_bound_variable(data):
    m, i, _, c1, obj, z = data

    int_x = gp.Variable(m, "int_x", type="integer", domain=[i])
    int_x.up[...] = 1e5

    obj[...] = gp.Sum(i, int_x) == z
    c1[...] = gp.Sum(i, int_x) <= 50

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIP",
        equations=[c1, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)
    with pytest.raises(ValidationError) as e1:
        qubo.transform()
    assert (
        "The Upper bound is greater than or equal to 1e+4, Quitting!"
        in str(e1.value)
    )

def test_qubo_nonlinear_constraint(data):
    m, i, x, _, obj, z = data

    j = gp.Alias(m, "j", alias_with=i)

    obj[...] = gp.Sum(i, x) == z

    nonlinear_cons = gp.Equation(m, "c2")   
    nonlinear_cons[...] = gp.Sum((i,j), x[i]*x[j] ) >= 2

    test1 = gp.Model(
        m,
        name="test1",
        problem="MIQCP",
        equations=[nonlinear_cons, obj],
        sense=gp.Sense.MAX,
        objective=z,
    )

    qubo = Qubo(test1, penalty=1)
    with pytest.raises(ValidationError) as e1:
        qubo.transform()
    assert (
        "There are non-linear constraints. Quitting."
        in str(e1.value)
    )