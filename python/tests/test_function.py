from pyariadne import *


def test_function():
    dp = DoublePrecision()

    x = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    y = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    f = x + y
    point = FloatDPApproximationVector([1, 2], dp)

    f(point)
    derivative(f, 0)

    domain = BoxDomainType([{0: 1}, {-1: 1}])
    sweeper = ThresholdSweeperDP(dp, 1e-8)
    tx = ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(domain, 0, sweeper)
    tf = 1 + tx

    derivative(tf, 0)
    antiderivative(tf, 0)
