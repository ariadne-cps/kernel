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

#include <memory>
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
        ARIADNE_TEST_CALL(test_affine_equality());
        ARIADNE_TEST_CALL(test_nonlinear_equality());
        ARIADNE_TEST_CALL(test_infeasible_equality());
        ARIADNE_TEST_CALL(test_dimension_precondition());
        ARIADNE_TEST_CALL(test_logging_diagnostics());
    }

    Void test_affine_equality() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction h({x[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        // There is a zero at the midpoint. A singular interval Newton
        // enclosure cannot justify either a true or a false certificate.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(domain,h)));
    }

    Void test_nonlinear_equality() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction h({sqr(x[0])-0.25_x});
        ExactBoxType domain({{0.0_x,1.0_x}});
        // The root x=1/2 exists, but uncertain Hessian and multiplier
        // ranges prevent an invertible initial interval correction.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(domain,h)));
    }

    Void test_infeasible_equality() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction h({x[0]-2});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        // x=2 lies outside D. A failed correction alone provides no
        // certificate of feasibility or infeasibility.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible_zero(domain,h)));
    }

    Void test_dimension_precondition() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction h({x[0]});
        ExactBoxType wrong_domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        // Reject the malformed public input before evaluating h or
        // constructing any of the interval Newton matrices.
        ARIADNE_TEST_FAIL(solver.feasible_zero(wrong_domain,h));
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
        test_affine_equality();
        test_nonlinear_equality();
        test_infeasible_equality();
    }
  private:
    IntervalNonlinearInteriorPointSolver solver;
};

Int main(Int argc, const char* argv[]) {
    if(not kernel_command_line_interface().acquire(argc,argv)) { return -1; }
    TestIntervalNonlinearInteriorPointSolver().test();
    return ARIADNE_TEST_FAILURES;
}
