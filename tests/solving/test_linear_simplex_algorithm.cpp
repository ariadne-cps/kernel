/***************************************************************************
 *            test_linear_simplex_algorithm.cpp
 *
 *  Copyright  2026  Ariadne contributors
 ****************************************************************************/

/*
 * This file is part of Ariadne.
 *
 * Ariadne is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Ariadne is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "macro/testing.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "solving/linear_simplex_algorithm.hpp"
#include "drawing/command_line_interface.hpp"

using namespace Ariadne;

class TestLinearSimplexSolver
{
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_pivot_to_better_vertex());
        ARIADNE_TEST_CALL(test_validated_strict_feasibility());
    }

    Void test_validated_strict_feasibility() {
        // With x0+x1=1/2 and both variables in [0,1], the initial
        // basis gives x0=1/2, a strictly feasible basic variable.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.validated_feasible(xl,xu,A,b)));
    }

    Void test_pivot_to_better_vertex() {
        // x0+x1=1, 0<=x0,x1<=1, minimise 2*x0+x1.
        // The default first-column basis starts at (1,0), whereas the
        // unique optimum (0,1) requires exchanging the basic variable.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.0_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Vector<FloatDP> c({2.0_x,1.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.feasible(xl,xu,A,b)));
        auto optimum=solver.minimise(c,xl,xu,A,b);
        ARIADNE_TEST_EQUAL(optimum.size(),2u);
        ARIADNE_TEST_ASSERT(definitely(optimum[0]<FloatDPBounds(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(optimum[1]>FloatDPBounds(0.75_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(
            optimum[0]+optimum[1]==FloatDPBounds(1.0_x,dp)));
    }
};

Int main(Int argc, const char* argv[])
{
    if (not kernel_command_line_interface().acquire(argc,argv)) return -1;
    TestLinearSimplexSolver().test();
    return ARIADNE_TEST_FAILURES;
}
