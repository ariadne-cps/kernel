/***************************************************************************
 *            test_nonlinear_interval_interior_point_method.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "macro/testing.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "geometry/box.hpp"
#include "solving/nonlinear_interior_point_method.hpp"
#include "drawing/command_line_interface.hpp"

using namespace Ariadne;

class TestIntervalNonlinearInteriorPointSolver {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_exact_midpoint());
        ARIADNE_TEST_CALL(test_empty_dimensions());
        ARIADNE_TEST_CALL(test_affine_contraction());
        ARIADNE_TEST_CALL(test_nonlinear_contraction());
        ARIADNE_TEST_CALL(test_singular_contraction());
        ARIADNE_TEST_CALL(test_empty_local_contraction());
        ARIADNE_TEST_CALL(test_infeasible_equality());
        ARIADNE_TEST_CALL(test_uncertainty_requires_certificate());
        ARIADNE_TEST_CALL(test_domain_preconditions());
        ARIADNE_TEST_CALL(test_logging_diagnostics());
    }

    Void test_exact_midpoint() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ARIADNE_TEST_ASSERT(definitely(solver.feasible_zero(
            ExactBoxType({{-1.0_x,1.0_x}}),ValidatedVectorMultivariateFunction({x[0]}))));
        ARIADNE_TEST_ASSERT(definitely(solver.feasible_zero(
            ExactBoxType({{0.0_x,1.0_x}}),ValidatedVectorMultivariateFunction({sqr(x[0])-0.25_x}))));
    }

    Void test_empty_dimensions() {
        // With no equalities, any point in the nonempty domain is feasible.
        ValidatedVectorMultivariateFunction vacuous(0u,1u);
        ARIADNE_TEST_ASSERT(definitely(solver.feasible_zero(
            ExactBoxType({{-1.0_x,1.0_x}}),vacuous)));
        // R^0 has one point. Exact zero is certified at that point;
        // a nonzero constant is not mistaken for a feasible equality.
        ValidatedVectorMultivariateFunction zero_function=ValidatedVectorMultivariateFunction::constant(
            0u,Vector<ValidatedNumber>({ValidatedNumber(0)}));
        ValidatedVectorMultivariateFunction nonzero_function=ValidatedVectorMultivariateFunction::constant(
            0u,Vector<ValidatedNumber>({ValidatedNumber(1)}));
        ARIADNE_TEST_ASSERT(definitely(solver.feasible_zero(ExactBoxType{},zero_function)));
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(ExactBoxType{},nonzero_function)));
    }

    Void test_affine_contraction() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        // Neither midpoint coordinate solves the equality. Both primal
        // and multiplier intervals must admit a complete Newton search.
        ValidatedVectorMultivariateFunction h({x[0]-0.25_x,x[1]+0.25_x});
        ExactBoxType domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        ARIADNE_TEST_ASSERT(definitely(solver.feasible_zero(domain,h)));
    }

    Void test_nonlinear_contraction() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        // sqrt(2) is not representable as a double. A numerical point
        // or small residual cannot replace the final interval certificate.
        ValidatedVectorMultivariateFunction h({sqr(x[0])-2});
        ExactBoxType domain({{1.0_x,2.0_x}});
        ARIADNE_TEST_ASSERT(definitely(solver.feasible_zero(domain,h)));
    }

    Void test_singular_contraction() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ExactBoxType domain({{-1.0_x,1.0_x}});
        // Both roots +/-1/4 exist, but the Hessian enclosure is singular.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            domain,ValidatedVectorMultivariateFunction({sqr(x[0])-0.0625_x}))));
        // A constant equality has a singular Schur complement instead.
        ValidatedVectorMultivariateFunction constant=ValidatedVectorMultivariateFunction::constant(
            1u,Vector<ValidatedNumber>({ValidatedNumber(1)}));
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(domain,constant)));
    }

    Void test_empty_local_contraction() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ExactBoxType domain({{-1.0_x,1.0_x}});
        // The zero is in D but outside the central search box. A local
        // exclusion cannot establish global infeasibility.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            domain,ValidatedVectorMultivariateFunction({x[0]-0.75_x}))));
        // Scaling preserves the root x=1/4, but the corresponding
        // barrier multiplier lies outside the local multiplier box.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            domain,ValidatedVectorMultivariateFunction({0.25_x*(x[0]-0.25_x)}))));
    }

    Void test_infeasible_equality() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            domain,ValidatedVectorMultivariateFunction({x[0]-2}))));
        // Even an arbitrarily small numerical residual is no certificate.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            domain,ValidatedVectorMultivariateFunction({0x1p-40_x*(x[0]-2)}))));
        // The nonlinear equality has roots +/-1/4, both outside [1,2].
        // Its Newton image is disjoint from the local primal box.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            ExactBoxType({{1.0_x,2.0_x}}),
            ValidatedVectorMultivariateFunction({sqr(x[0])-0.0625_x}))));
    }

    Void test_uncertainty_requires_certificate() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        auto uncertain=ValidatedScalarMultivariateFunction::constant(
            1u,ValidatedNumber(FloatDPBounds(0.25_x,0.75_x,dp)));
        ValidatedVectorMultivariateFunction h({x[0]-uncertain});
        // The local contraction retains only part of the possible root
        // range. It cannot certify all represented constant values, and
        // an arbitrary numerical midpoint is not an equality certificate.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(
            ExactBoxType({{-1.0_x,1.0_x}}),h)));
    }

    Void test_domain_preconditions() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction h({x[0]});
        ExactBoxType wrong_domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        ExactBoxType unbounded({{-1.0_x,inf}});
        ExactBoxType singleton({{0.0_x,0.0_x}});
        // Each input violates a different public requirement before
        // constructing or evaluating any barrier reciprocal.
        ARIADNE_TEST_FAIL(solver.feasible_zero(wrong_domain,h));
        ARIADNE_TEST_FAIL(solver.feasible_zero(unbounded,h));
        ARIADNE_TEST_FAIL(solver.feasible_zero(singleton,h));
    }

    Void test_logging_diagnostics() {
        auto& configuration=Logger::instance().configuration();
        struct RestoreVerbosity {
            unsigned int value;
            ~RestoreVerbosity() {
                Logger::instance().configuration().set_verbosity(value);
            }
        } restore{configuration.verbosity()};
        configuration.set_verbosity(10u);
        test_exact_midpoint();
        test_affine_contraction();
        test_nonlinear_contraction();
        test_singular_contraction();
        test_empty_local_contraction();
        test_infeasible_equality();
        test_uncertainty_requires_certificate();
    }
  private:
    IntervalNonlinearInteriorPointSolver solver;
};

Int main(Int argc, const char* argv[]) {
    if(not kernel_command_line_interface().acquire(argc,argv)) { return -1; }
    TestIntervalNonlinearInteriorPointSolver().test();
    return ARIADNE_TEST_FAILURES;
}
