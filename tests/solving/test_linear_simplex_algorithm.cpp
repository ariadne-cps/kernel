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
        ARIADNE_TEST_CALL(test_validated_feasibility_bound_flip());
        ARIADNE_TEST_CALL(test_hotstarted_basis_initialisation());
        ARIADNE_TEST_CALL(test_hotstarted_optimisation_from_alternate_basis());
        ARIADNE_TEST_CALL(test_hotstarted_minimise_partial_permutation());
        ARIADNE_TEST_CALL(test_consistency_check_nonsquare_system());
    }

    Void test_consistency_check_nonsquare_system() {
        // An equality system may have more variables than equations.
        // The check must form A*x-b, not A*b-x (the latter has
        // incompatible dimensions in this 1-by-2 example).
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.0_x},dp);
        Vector<FloatDPBounds> x({0.25_x,0.75_x},dp);

        ARIADNE_TEST_ASSERT(definitely(
            x[0]+x[1]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,b,x));
    }

    Void test_hotstarted_minimise_partial_permutation() {
        // A warm start may specify only the index of each basic variable.
        // Extending that permutation must preserve the feasible starting
        // basis and yield the same optimum as the complete permutation.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.0_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Vector<FloatDP> c({1.0_x,2.0_x},dp);

        Array<Slackness> short_vt(2u);
        short_vt[0]=Slackness::LOWER;
        short_vt[1]=Slackness::BASIS;
        Array<SizeType> short_p(1u);
        short_p[0]=1u;
        Matrix<FloatDPBounds> short_B=Matrix<FloatDPBounds>::identity(1u,dp);

        auto short_optimum=solver.hotstarted_minimise(c,xl,xu,A,b,short_vt,short_p,short_B);
        ARIADNE_TEST_EQUAL(short_p.size(),2u);
        ARIADNE_TEST_EQUAL(short_p[0],1u);
        ARIADNE_TEST_EQUAL(short_p[1],0u);
        ARIADNE_TEST_ASSERT(definitely(short_optimum[0]>FloatDPBounds(0.75_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(short_optimum[1]<FloatDPBounds(0.25_x,dp)));

        Array<Slackness> full_vt(2u);
        full_vt[0]=Slackness::LOWER;
        full_vt[1]=Slackness::BASIS;
        Array<SizeType> full_p(2u);
        full_p[0]=1u;
        full_p[1]=0u;
        Matrix<FloatDPBounds> full_B=Matrix<FloatDPBounds>::identity(1u,dp);

        auto full_optimum=solver.hotstarted_minimise(c,xl,xu,A,b,full_vt,full_p,full_B);
        ARIADNE_TEST_ASSERT(definitely(
            abs(short_optimum[0]-full_optimum[0])<FloatDPBounds(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(
            abs(short_optimum[1]-full_optimum[1])<FloatDPBounds(0.25_x,dp)));
    }

    Void test_hotstarted_optimisation_from_alternate_basis() {
        // Start from (0,1) with x1 basic and x0 at its lower bound.
        // Minimising x0+2*x1 must move to (1,0), even though the
        // supplied basis describes a feasible starting point.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.0_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Vector<FloatDP> c({1.0_x,2.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::LOWER;
        vt[1]=Slackness::BASIS;

        auto optimum=solver.hotstarted_minimise(c,xl,xu,A,b,vt);
        ARIADNE_TEST_EQUAL(optimum.size(),2u);
        ARIADNE_TEST_ASSERT(definitely(optimum[0]>FloatDPBounds(0.75_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(optimum[1]<FloatDPBounds(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(
            optimum[0]+optimum[1]==FloatDPBounds(1.0_x,dp)));
    }

    Void test_hotstarted_basis_initialisation() {
        // With x0+x1=1/2 and x1 fixed at its lower bound by the
        // selected basis, the reconstructed point must be (1/2,0).
        // The same basis must certify the strictly feasible basic x0.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::LOWER;

        auto x=solver.compute_x(xl,xu,A,b,vt);
        ARIADNE_TEST_EQUAL(x.size(),2u);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(solver.hotstarted_feasible(xl,xu,A,b,vt)));
        ARIADNE_TEST_EQUAL(vt[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::LOWER);
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

    Void test_validated_feasibility_bound_flip() {
        // The first basic variable violates its upper bound (x0=1.25).
        // Moving the nonbasic x1 from LOWER to UPPER leaves x0=0.25,
        // strictly inside [0,1]. No basis exchange is needed.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.25_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::LOWER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;

        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_EQUAL(vt[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::UPPER);
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(xl,xu,A,b,vt,p));
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
