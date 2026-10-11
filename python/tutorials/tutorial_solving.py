#!/usr/bin/python3

from pyariadne import *


def section(title):
    print(f"\n== {title} ==")


def show(label, value):
    print(f"{label}: {value}")


def exact_box(bounds):
    return FloatDPExactBox([
        {exact(lower): exact(upper)}
        for lower, upper in bounds
    ])


def tutorial_solving():
    dp = DoublePrecision()

    section("Algebraic equations with Interval Newton")

    solver = IntervalNewtonSolver(1e-10, 16)
    show("solver", solver)
    show("maximum error", solver.maximum_error())
    show("maximum number of steps", solver.maximum_number_of_steps())

    x = EffectiveScalarMultivariateFunction.coordinate(1, 0)
    root_function = EffectiveVectorMultivariateFunction([
        (x * x + 1) * x - 1,
    ])
    root_domain = exact_box([
        (0, 1),
    ])

    root = solver.solve(root_function, root_domain)
    show("root function", root_function)
    show("search domain", root_domain)
    show("validated root", root)

    root_via_zero = solver.zero(root_function, root_domain)
    show("root via zero()", root_via_zero)

    section("Fixed points")

    contraction = EffectiveVectorMultivariateFunction([
        (x + 1) / 2,
    ])
    fixed_point_domain = exact_box([
        (0, 2),
    ])
    fixed_point = solver.fixed_point(contraction, fixed_point_domain)

    show("fixed-point map", contraction)
    show("fixed-point enclosure", fixed_point)

    section("Implicit equations")

    a = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    y = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    implicit_equation = EffectiveVectorMultivariateFunction([
        y * y - a,
    ])
    parameter_domain = exact_box([
        (0.75, 1.25),
    ])
    solution_range = exact_box([
        (0.75, 1.25),
    ])

    implicit_solution = solver.implicit(
        implicit_equation,
        parameter_domain,
        solution_range,
    )

    show("implicit equation", implicit_equation)
    show("parameter domain", parameter_domain)
    show("solution search range", solution_range)
    show("implicit solution patch", implicit_solution)
    show("implicit solution range", implicit_solution.range())

    section("Constraint solving and domain contraction")

    c0 = ValidatedScalarMultivariateFunction.coordinate(2, 0)
    c1 = ValidatedScalarMultivariateFunction.coordinate(2, 1)
    constraint_function = ValidatedVectorMultivariateFunction([
        c0 + c1,
        c0 - c1,
    ])
    constraint_domain = exact_box([
        (0, 2),
        (0, 2),
    ])
    constraint_codomain = exact_box([
        (0.875, 1.125),
        (-0.25, 0.25),
    ])

    constraint_solver = ConstraintSolver()
    feasibility, candidate = constraint_solver.feasible(
        constraint_domain,
        constraint_function,
        constraint_codomain,
    )

    show("constraint function", constraint_function)
    show("constraint domain", constraint_domain)
    show("constraint codomain", constraint_codomain)
    show("feasibility result", feasibility)
    show("candidate point", candidate)

    proved_empty, contracted_domain = constraint_solver.reduce_box(
        constraint_domain,
        constraint_function,
        constraint_codomain,
    )
    show("proved empty during reduction", proved_empty)
    show("contracted domain", contracted_domain)

    section("Validated ODE integration with Taylor-Picard")

    integrator = TaylorPicardIntegrator(1e-8)
    integrator.set_minimum_temporal_order(2)
    integrator.set_maximum_temporal_order(16)
    integrator.set_step_maximum_error(1e-8)

    show("integrator", integrator)
    show("minimum temporal order", integrator.minimum_temporal_order())
    show("maximum temporal order", integrator.maximum_temporal_order())
    show("step maximum error", integrator.step_maximum_error())

    state = EffectiveScalarMultivariateFunction.coordinate(1, 0)
    logistic = EffectiveVectorMultivariateFunction([
        state * (1 - state),
    ])
    initial_states = exact_box([
        (0.25, 0.5),
    ])
    step = Dyadic("0.125")

    flow = integrator.flow_step(logistic, initial_states, step)

    show("vector field", logistic)
    show("initial-state box", initial_states)
    show("requested step", step)
    show("validated flow patch", flow)
    show("flow domain", flow.domain())
    show("flow range", flow.range())

    section("Linear programming with the simplex solver")

    simplex = FloatDPSimplexSolver()
    objective = FloatDPVector([1, 2], dp)
    lower = FloatDPVector([0, 0], dp)
    upper = FloatDPVector([1, 1], dp)
    equality_matrix = FloatDPMatrix([[1, 1]], dp)
    equality_rhs = FloatDPVector([1], dp)

    lp_feasible = simplex.feasible(
        lower, upper, equality_matrix, equality_rhs
    )
    lp_solution = simplex.minimise(
        objective, lower, upper, equality_matrix, equality_rhs
    )

    show("LP feasible", lp_feasible)
    show("LP solution", lp_solution)

    section("Nonlinear constrained optimisation")

    q0 = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    q1 = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    five_eighths = EffectiveNumber(Rational("5/8"))
    three_eighths = EffectiveNumber(Rational("3/8"))

    objective_function = (
        sqr(q0 - five_eighths)
        + sqr(q1 - three_eighths)
    )
    band_constraint = EffectiveVectorMultivariateFunction([
        q0 + q1,
    ])
    optimisation_domain = exact_box([
        (0, 1),
        (0, 1),
    ])
    optimisation_codomain = exact_box([
        (0.75, 1.25),
    ])

    optimiser = NonlinearInteriorPointOptimiser()
    # Feasibility as a separate task was already demonstrated with
    # ConstraintSolver.  The current NonlinearInteriorPointOptimiser
    # feasibility_step path is not implemented, so this section focuses on
    # the supported minimisation workflow.
    optimum = optimiser.minimise(
        objective_function,
        optimisation_domain,
        band_constraint,
        optimisation_codomain,
    )

    show("validated local optimum", optimum)
    show("constraint at optimum", band_constraint(optimum))

    section("Choosing a solving tool")
    print("IntervalNewtonSolver: roots, fixed points and implicit functions.")
    print("ConstraintSolver: feasibility tests and interval-domain contraction.")
    print("TaylorPicardIntegrator: validated flow steps for ODEs.")
    print("FloatDPSimplexSolver: bounded linear programs with linear equalities.")
    print("NonlinearInteriorPointOptimiser: constrained nonlinear minimisation with interior-point barriers.")


if __name__ == "__main__":
    tutorial_solving()
