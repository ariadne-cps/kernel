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

    // EquationSolverBase exposes a common interface for algebraic solvers.  We use
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
    // writing code against EquationSolverInterface and makes the intended equation
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

    section("Linear programming with the simplex solver");

    // Solve
    //   minimise x0 + 2*x1
    //   subject to x0 + x1 = 1
    //              0 <= x0,x1 <= 1.
    // The optimum is attained at (1,0).
    LinearSimplexSolver<FloatDP> simplex;
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

    // Minimise the squared distance from the target (5/8,3/8) over
    // [0,1]^2, while requiring x+y to stay in the non-degenerate band
    // [3/4,5/4].  The target itself is feasible and lies in the interior of
    // both the search domain and the constraint band, which is the natural
    // use case for the barrier-based interior-point implementation.
    //
    // Feasibility as a separate task was already demonstrated above with
    // ConstraintSolver.  NonlinearInteriorPointSolver::feasible() is not
    // used here because its current feasibility_step path is not implemented.
    NonlinearInteriorPointSolver optimiser;
    auto q = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction objective_function =
        sqr(q[0] - 0.625_x) + sqr(q[1] - 0.375_x);
    EffectiveVectorMultivariateFunction band_constraint({
        q[0] + q[1]
    });
    ExactBoxType optimisation_domain({
        {0.0_x, 1.0_x},
        {0.0_x, 1.0_x}
    });
    ExactBoxType optimisation_codomain({
        {0.75_x, 1.25_x}
    });

    FloatDPBoundsVector optimum =
        optimiser.minimise(
            objective_function,
            optimisation_domain,
            band_constraint,
            optimisation_codomain);

    show("validated local optimum", optimum);
    show("constraint at optimum", band_constraint(optimum));

    section("Choosing a solving tool");
    cout
        << "IntervalNewtonSolver: roots, fixed points and implicit functions.\n"
        << "ConstraintSolver: feasibility tests and interval-domain contraction.\n"
        << "LinearSimplexSolver: bounded linear programs with linear equalities.\n"
        << "NonlinearInteriorPointSolver: constrained nonlinear minimisation with interior-point barriers.\n";

    return 0;
}
