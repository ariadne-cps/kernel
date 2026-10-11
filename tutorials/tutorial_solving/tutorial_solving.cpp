#include <iostream>
#include <string>

#include "ariadne-kernel.hpp"

using namespace Ariadne;
using std::cout;

namespace {

void section(const std::string& title) {
    cout << "\n== " << title << " ==\n";
}

template<class T>
void show(const std::string& label, const T& value) {
    cout << label << ": " << value << "\n";
}

} // namespace

int main() {
    section("Algebraic equations with Interval Newton");

    // SolverBase exposes a common interface for algebraic solvers.  We use
    // IntervalNewtonSolver here as the representative implementation: the
    // tutorial is about the solving workflow, not comparing algorithms.
    IntervalNewtonSolver solver(1e-10, 16);
    show("maximum error", solver.maximum_error());
    show("maximum number of steps", solver.maximum_number_of_steps());

    auto x = EffectiveScalarMultivariateFunction::coordinate(1, 0);
    EffectiveVectorMultivariateFunction root_function({
        (x*x + 1) * x - 1
    });
    ExactBoxType root_domain({
        {0.0_x, 1.0_x}
    });

    // solve() searches for a unique transverse zero in the supplied box and
    // returns a validated enclosure of the root.
    FloatDPBoundsVector root = solver.solve(root_function, root_domain);
    show("root function", root_function);
    show("search domain", root_domain);
    show("validated root", root);

    // zero() is the explicit root-solving entry point.  It is useful when
    // writing code against SolverInterface and makes the intended equation
    // f(x)=0 clear at the call site.
    FloatDPBoundsVector root_via_zero = solver.zero(root_function, root_domain);
    show("root via zero()", root_via_zero);

    section("Fixed points");

    // fixed_point() solves f(x)=x.  For f(x)=(x+1)/2 the unique fixed point
    // in [0,2] is x=1.
    EffectiveVectorMultivariateFunction contraction({
        (x + 1) / 2
    });
    ExactBoxType fixed_point_domain({
        {0.0_x, 2.0_x}
    });
    FloatDPBoundsVector fixed_point =
        solver.fixed_point(contraction, fixed_point_domain);

    show("fixed-point map", contraction);
    show("fixed-point enclosure", fixed_point);

    section("Implicit equations");

    // The implicit solver turns an equation F(a,y)=0 into a validated
    // function y=h(a).  Here F(a,y)=y^2-a and the selected branch is near
    // y=1, so the result is a validated square-root function on a parameter
    // interval around one.
    auto z = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction implicit_equation({
        z[1] * z[1] - z[0]
    });
    ExactBoxType parameter_domain({
        {0.75_x, 1.25_x}
    });
    ExactBoxType solution_range({
        {0.75_x, 1.25_x}
    });

    ValidatedVectorMultivariateFunctionPatch implicit_solution =
        solver.implicit(implicit_equation, parameter_domain, solution_range);

    show("implicit equation", implicit_equation);
    show("parameter domain", parameter_domain);
    show("solution search range", solution_range);
    show("implicit solution patch", implicit_solution);
    show("implicit solution range", implicit_solution.range());

    section("Constraint solving and domain contraction");

    // ConstraintSolver answers a different question: does the image g(D)
    // intersect a requested codomain C?  It can also contract the domain by
    // propagating those constraints before a more expensive algorithm runs.
    ConstraintSolver constraint_solver;
    auto c = ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction constraint_function({
        c[0] + c[1],
        c[0] - c[1]
    });
    ExactBoxType constraint_domain({
        {0.0_x, 2.0_x},
        {0.0_x, 2.0_x}
    });
    ExactBoxType constraint_codomain({
        {0.875_x, 1.125_x},
        {-0.25_x, 0.25_x}
    });

    auto feasibility =
        constraint_solver.feasible(
            constraint_domain, constraint_function, constraint_codomain);

    show("constraint function", constraint_function);
    show("constraint domain", constraint_domain);
    show("constraint codomain", constraint_codomain);
    show("feasibility result", feasibility.first);
    show("candidate point", feasibility.second);

    UpperBoxType contracted_domain(constraint_domain);
    Bool proved_empty =
        constraint_solver.reduce(
            contracted_domain, constraint_function, constraint_codomain);

    show("proved empty during reduction", proved_empty);
    show("contracted domain", contracted_domain);

    section("Validated ODE integration with Taylor-Picard");

    // Only TaylorPicardIntegrator is shown here.  The important part of the
    // module API is the validated flow-step abstraction; alternative
    // integrators implement the same role with different numerical methods.
    TaylorPicardIntegrator integrator(1e-8);
    integrator.set_minimum_temporal_order(2);
    integrator.set_maximum_temporal_order(16);
    integrator.set_step_maximum_error(1e-8);

    show("minimum temporal order", integrator.minimum_temporal_order());
    show("maximum temporal order", integrator.maximum_temporal_order());
    show("step maximum error", integrator.step_maximum_error());

    auto state = EffectiveScalarMultivariateFunction::coordinate(1, 0);
    EffectiveVectorMultivariateFunction logistic({
        state * (1 - state)
    });
    ExactBoxType initial_states({
        {0.25_x, 0.5_x}
    });
    StepSizeType step = 0.125_dy;

    FlowStepModelType flow =
        integrator.flow_step(logistic, initial_states, step);

    // The returned object is a validated function patch in (initial state,
    // time).  Its domain records both the initial-state box and the time
    // interval actually covered by the step.
    show("vector field", logistic);
    show("initial-state box", initial_states);
    show("requested step", step);
    show("validated flow patch", flow);
    show("flow domain", flow.domain());
    show("flow range", flow.range());

    section("Linear programming with the simplex solver");

    // Solve
    //   minimise x0 + 2*x1
    //   subject to x0 + x1 = 1
    //              0 <= x0,x1 <= 1.
    // The optimum is attained at (1,0).
    SimplexSolver<FloatDP> simplex;
    Vector<FloatDP> objective({1.0_x, 2.0_x}, dp);
    Vector<FloatDP> lower({0.0_x, 0.0_x}, dp);
    Vector<FloatDP> upper({1.0_x, 1.0_x}, dp);
    Matrix<FloatDP> equality_matrix({{1.0_x, 1.0_x}}, dp);
    Vector<FloatDP> equality_rhs({1.0_x}, dp);

    ValidatedKleenean lp_feasible =
        simplex.feasible(lower, upper, equality_matrix, equality_rhs);
    FloatDPBoundsVector lp_solution =
        simplex.minimise(
            objective, lower, upper, equality_matrix, equality_rhs);

    show("LP feasible", lp_feasible);
    show("LP solution", lp_solution);

    section("Nonlinear constrained optimisation");

    // Minimise x^2+y^2 over [0,1]^2 subject to x+y=1.  The optimum is
    // (1/2,1/2).  NonlinearInteriorPointOptimiser returns a validated box
    // containing a feasible local optimum.
    NonlinearInteriorPointOptimiser optimiser;
    auto q = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction objective_function =
        sqr(q[0]) + sqr(q[1]);
    EffectiveVectorMultivariateFunction equality_constraint({
        q[0] + q[1]
    });
    ExactBoxType optimisation_domain({
        {0.0_x, 1.0_x},
        {0.0_x, 1.0_x}
    });
    ExactBoxType optimisation_codomain({
        {1.0_x, 1.0_x}
    });

    ValidatedKleenean nonlinear_feasible =
        optimiser.feasible(
            optimisation_domain,
            equality_constraint,
            optimisation_codomain);

    FloatDPBoundsVector optimum =
        optimiser.minimise(
            objective_function,
            optimisation_domain,
            equality_constraint,
            optimisation_codomain);

    show("nonlinear problem feasible", nonlinear_feasible);
    show("validated local optimum", optimum);
    show("constraint at optimum", equality_constraint(optimum));

    section("Choosing a solving tool");
    cout
        << "IntervalNewtonSolver: roots, fixed points and implicit functions.\n"
        << "ConstraintSolver: feasibility tests and interval-domain contraction.\n"
        << "TaylorPicardIntegrator: validated flow steps for ODEs.\n"
        << "SimplexSolver: bounded linear programs with linear equalities.\n"
        << "NonlinearInteriorPointOptimiser: constrained nonlinear feasibility and minimisation.\n";

    return 0;
}
