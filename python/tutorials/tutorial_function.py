#!/usr/bin/python3

from pyariadne import *


def section(title):
    print(f"\n== {title} ==")


def show(label, value):
    print(f"{label}: {value}")


def tutorial_function():
    dp = DoublePrecision()

    section("Effective functions")

    x0 = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    x1 = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    f = sqrt(sqr(x0) + sqr(x1))
    g = EffectiveVectorMultivariateFunction([x0 + x1, x0 * x1])

    show("f", f)
    show("g", g)

    approximate_point = FloatDPApproximationVector([3, 4], dp)
    bounded_point = FloatDPBoundsVector([3, 4], dp)

    show("f(approximate_point)", f(approximate_point))
    show("f(bounded_point)", f(bounded_point))
    show("g(approximate_point)", g(approximate_point))
    show("g(bounded_point)", evaluate(g, bounded_point))

    show("derivative(f,0)", derivative(f, 0))

    identity2 = EffectiveVectorMultivariateFunction.identity(2)
    show("compose(f,identity2)", compose(f, identity2))
    show("join(f,g)", join(f, g))

    section("Taylor function models")

    domain = BoxDomainType([{-1: +1}, {-1: +1}])
    sweeper = ThresholdSweeperDP(dp, 1e-8)

    tx = ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(
        domain, 0, sweeper
    )
    ty = ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(
        domain, 1, sweeper
    )
    tidentity = ValidatedVectorMultivariateTaylorFunctionModelDP.identity(
        domain, sweeper
    )

    tf = 1 + 2*tx + 3*tx*ty + sqr(tx)
    tv = ValidatedVectorMultivariateTaylorFunctionModelDP([tx + ty, tx*ty])

    show("domain", tf.domain())
    show("codomain", tf.codomain())
    show("range", tf.range())
    show("tf", tf)
    show("tv", tv)

    show("derivative(tf,0)", derivative(tf, 0))
    show("antiderivative(tf,0)", antiderivative(tf, 0))

    subdomain = BoxDomainType([
        {exact(-0.5): exact(0.5)},
        {exact(-0.25): exact(0.25)},
    ])
    show("restriction(tf,subdomain)", restriction(tf, subdomain))

    outer_domain = tidentity.codomain()
    polynomial_function = sqr(x0) + x1
    outer_model = ValidatedScalarMultivariateTaylorFunctionModelDP(
        outer_domain, polynomial_function, sweeper
    )

    show(
        "compose(polynomial_function,tidentity)",
        compose(polynomial_function, tidentity),
    )
    show(
        "compose(outer_model,tidentity)",
        compose(outer_model, tidentity),
    )

    coarse = tx + FloatDPBounds(-1, 1, dp)
    show("inconsistent(tx,coarse)", inconsistent(tx, coarse))
    show("refines(tx,coarse)", refines(tx, coarse))
    show("refinement(tx,coarse)", refinement(tx, coarse))


if __name__ == "__main__":
    tutorial_function()
