from pyariadne import *

def real(numerator, denominator=1):
    return Real(Rational(numerator, denominator))



def test_satisfiability():
    x = RealVariable("x")
    space = RealSpace([x])
    domain = FloatDPExactBox([
        {exact(-1): exact(1)},
    ])
    solver = SmtSolver(SmtSolverConfiguration(exact(0.125)))
    result = solver.solve(space, domain, RealExpression(x) > real(0))
    assert result.is_epsilon_sat()
    assert result.has_witness()
