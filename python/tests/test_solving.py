from pyariadne import *


def test_solving():
    solver = IntervalNewtonSolver(1e-8, 12)
    str(solver)

    integrator = TaylorPicardIntegrator(1e-8)
    integrator.minimum_temporal_order()
    integrator.maximum_temporal_order()
