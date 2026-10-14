from pyariadne import *


def test_solving():
    solver = IntervalNewtonSolver(1e-8, 12)
    str(solver)
