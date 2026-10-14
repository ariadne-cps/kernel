#!/usr/bin/python3

from pyariadne import *


def section(title):
    print(f"\n== {title} ==")


def show(label, value):
    print(f"{label}: {value}")


def tutorial_satisfiability():
    section("Bounded real epsilon-SMT")

    x = RealVariable("x")
    y = RealVariable("y")
    ex = RealExpression(x)
    ey = RealExpression(y)
    space = RealSpace([x, y])

    domain = FloatDPExactBox([
        {exact(-2), exact(2)},
        {exact(-2), exact(2)},
    ])

    formula = (
        sqr(ex) + sqr(ey) <= real(1)
    ) & (
        ex >= real(0)
    )

    solver = SmtSolver(SmtSolverConfiguration(exact(0.125)))
    result = solver.solve(space, domain, formula)

    show("formula", formula)
    show("status", result.status())
    show("epsilon-satisfiable", result.is_epsilon_sat())
    if result.has_witness():
        show("witness box", result.witness())

    statistics = result.statistics()
    show("boxes processed", statistics.boxes_processed)
    show("theory checks", statistics.theory_checks)

    section("An unsatisfiable conjunction")

    contradiction = (ex >= real(1)) & (ex <= real(0))
    unsat = solver.solve(space, domain, contradiction)

    show("formula", contradiction)
    show("status", unsat.status())
    show("unsatisfiable", unsat.is_unsat())

    section("Boolean combinations")

    disjunction = (ex < real(-1)) | (ex > real(1))
    boolean_result = solver.solve(space, domain, disjunction)

    show("formula", disjunction)
    show("status", boolean_result.status())


if __name__ == "__main__":
    tutorial_satisfiability()
