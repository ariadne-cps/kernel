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

#include <sstream>
#include "macro/testing.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "solving/linear_simplex_algorithm.hpp"
#include "solving/linear_programming.hpp"
#include "drawing/command_line_interface.hpp"

using namespace Ariadne;

class TestLinearSimplexSolver
{
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_linear_programming_exception_contracts());
        ARIADNE_TEST_CALL(test_public_simplex_input_contracts());
        ARIADNE_TEST_CALL(test_slackness_output());
        ARIADNE_TEST_CALL(test_pivot_to_better_vertex());
        ARIADNE_TEST_CALL(test_verify_feasibility_basic_statuses());
        ARIADNE_TEST_CALL(test_verify_feasibility_dual_signs());
        ARIADNE_TEST_CALL(test_upper_only_nonbasic_bounds());
        ARIADNE_TEST_CALL(test_approximate_upper_entry_basis_exchange());
        ARIADNE_TEST_CALL(test_unbounded_nonbasic_initialisation());
        ARIADNE_TEST_CALL(test_lpstep_without_leaving_variable());
        ARIADNE_TEST_CALL(test_approximate_lpstep_without_leaving_variable());
        ARIADNE_TEST_CALL(test_single_variable_infeasibility());
        ARIADNE_TEST_CALL(test_lpstep_upper_bound_flip());
        ARIADNE_TEST_CALL(test_lpstep_upper_entry());
        ARIADNE_TEST_CALL(test_upper_entry_tied_upper_bounds());
        ARIADNE_TEST_CALL(test_validated_multiple_upper_violations());
        ARIADNE_TEST_CALL(test_degenerate_ratio_tie_breaking());
        ARIADNE_TEST_CALL(test_validated_upper_entry_bound_flip());
        ARIADNE_TEST_CALL(test_validated_upper_entry_basis_exchange());
        ARIADNE_TEST_CALL(test_validated_strict_feasibility());
        ARIADNE_TEST_CALL(test_validated_feasibility_stationary_at_boundary());
        ARIADNE_TEST_CALL(test_validated_feasibility_bound_flip());
        ARIADNE_TEST_CALL(test_validated_feasibility_basis_exchange());
        ARIADNE_TEST_CALL(test_hotstarted_basis_initialisation());
        ARIADNE_TEST_CALL(test_hotstarted_optimisation_from_alternate_basis());
        ARIADNE_TEST_CALL(test_hotstarted_minimise_partial_permutation());
        ARIADNE_TEST_CALL(test_invalid_partial_basis_permutations());
        ARIADNE_TEST_CALL(test_consistency_check_variable_types());
        ARIADNE_TEST_CALL(test_consistency_check_nonsquare_system());
        ARIADNE_TEST_CALL(test_hotstarted_feasibility_optional_basis());
        ARIADNE_TEST_CALL(test_compute_basis_tall_matrix_precondition());
        ARIADNE_TEST_CALL(test_compute_basis_unit_columns());
        ARIADNE_TEST_CALL(test_compute_basis_singular_after_elimination());
        ARIADNE_TEST_CALL(test_compute_basis_largest_nonunit_pivot());
        ARIADNE_TEST_CALL(test_compute_basis_rational_nonunit());
        ARIADNE_TEST_CALL(test_compute_basis_nonunit_two_by_two());
        ARIADNE_TEST_CALL(test_compute_basis_column_pivoting());
    }


    Void test_public_simplex_input_contracts() {
        // A single independent equality in two bounded variables has a
        // valid basis, so each malformed input below violates exactly one
        // documented dimension or basis invariant.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        Vector<FloatDP> c({1.0_x,2.0_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Vector<FloatDP> short_vector({0.0_x},dp);
        Vector<FloatDP> long_vector({0.0_x,0.0_x,0.0_x},dp);

        ARIADNE_TEST_FAIL(solver.feasible(xl,xu,A,long_vector));
        ARIADNE_TEST_FAIL(solver.feasible(short_vector,xu,A,b));
        ARIADNE_TEST_FAIL(solver.feasible(xl,short_vector,A,b));
        Array<Slackness> types(2u);
        types[0]=Slackness::BASIS;
        types[1]=Slackness::LOWER;
        ARIADNE_TEST_FAIL(solver.hotstarted_feasible(xl,xu,A,long_vector,types));
        ARIADNE_TEST_FAIL(solver.hotstarted_feasible(short_vector,xu,A,b,types));
        ARIADNE_TEST_FAIL(solver.hotstarted_feasible(xl,short_vector,A,b,types));

        ARIADNE_TEST_FAIL(solver.minimise(c,xl,xu,A,long_vector));
        ARIADNE_TEST_FAIL(solver.minimise(short_vector,xl,xu,A,b));
        ARIADNE_TEST_FAIL(solver.minimise(c,short_vector,xu,A,b));
        ARIADNE_TEST_FAIL(solver.minimise(c,xl,short_vector,A,b));

        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,long_vector,types));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(short_vector,xl,xu,A,b,types));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,short_vector,xu,A,b,types));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,short_vector,A,b,types));
        Array<Slackness> wrong_size(1u);
        wrong_size[0]=Slackness::BASIS;
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,b,wrong_size));
        Array<Slackness> no_basis(2u);
        no_basis[0]=Slackness::LOWER;
        no_basis[1]=Slackness::UPPER;
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,b,no_basis));

        Array<SizeType> permutation(2u);
        permutation[0]=0u;
        permutation[1]=1u;
        Matrix<FloatDPBounds> inverse=Matrix<FloatDPBounds>::identity(1u,dp);
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,long_vector,types,permutation,inverse));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(short_vector,xl,xu,A,b,types,permutation,inverse));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,short_vector,xu,A,b,types,permutation,inverse));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,short_vector,A,b,types,permutation,inverse));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,b,wrong_size,permutation,inverse));
        Array<SizeType> wrong_permutation(3u);
        wrong_permutation[0]=0u;
        wrong_permutation[1]=1u;
        wrong_permutation[2]=2u;
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,b,types,wrong_permutation,inverse));
        Matrix<FloatDPBounds> wrong_rows(2u,1u,FloatDPBounds(0.0_x,dp));
        Matrix<FloatDPBounds> wrong_columns(1u,2u,FloatDPBounds(0.0_x,dp));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,b,types,permutation,wrong_rows));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(c,xl,xu,A,b,types,permutation,wrong_columns));
    }

    Void test_linear_programming_exception_contracts() {
        const DegenerateFeasibilityProblemException degenerate;
        const UnboundedLinearProgram unbounded(StringType("unbounded programme"));
        const InfeasibleLinearProgram infeasible(StringType("infeasible programme"));
        ARIADNE_TEST_EQUAL(std::string(degenerate.what()),std::string(""));
        ARIADNE_TEST_EQUAL(std::string(unbounded.what()),std::string("unbounded programme"));
        ARIADNE_TEST_EQUAL(std::string(infeasible.what()),std::string("infeasible programme"));
        const std::runtime_error& base=infeasible;
        ARIADNE_TEST_EQUAL(std::string(base.what()),std::string("infeasible programme"));
    }

    Void test_slackness_output() {
        // Exercise all named variable states and the documented fallback.
        std::ostringstream stream;
        stream << Slackness::BASIS << Slackness::LOWER
               << Slackness::UPPER << Slackness::FIXED;
        ARIADNE_TEST_EQUAL(stream.str(),std::string("BLUE"));

        std::ostringstream unknown;
        unknown << static_cast<Slackness>(3);
        ARIADNE_TEST_EQUAL(unknown.str(),std::string("?"));
    }

    Void test_compute_basis_column_pivoting() {
        // The first column cannot pivot on row zero. The algorithm
        // must exchange the two columns and return the inverse in
        // the resulting basis order.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{0.0_x,2.0_x},{1.0_x,1.0_x}},dp);

        auto basis=solver.compute_basis(A);
        const auto& p=basis.first;
        const auto& B=basis.second;
        ARIADNE_TEST_EQUAL(p.size(),2u);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(p[1],0u);
        ARIADNE_TEST_ASSERT(definitely(B[0][0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[0][1]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[1][0]==FloatDPBounds(-0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[1][1]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,p,B));
    }

    Void test_compute_basis_tall_matrix_precondition() {
        // A two-by-one matrix cannot contain a full row basis.
        // The public entry point must reject the dimensions before
        // attempting pivot selection or indexing its permutation.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> tall(2u,1u,FloatDP(1.0_x,dp));
        ARIADNE_TEST_FAIL(solver.compute_basis(tall));
    }

    Void test_compute_basis_unit_columns() {
        // A unit pivot with no nonzero entries below it must be accepted.
        // The identity matrix needs no column exchange or elimination.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x},{0.0_x,1.0_x}},dp);

        auto basis=solver.compute_basis(A);
        const auto& p=basis.first;
        const auto& B=basis.second;
        ARIADNE_TEST_EQUAL(p.size(),2u);
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],1u);
        ARIADNE_TEST_ASSERT(definitely(B[0][0]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[0][1]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[1][0]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[1][1]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,p,B));
    }

    Void test_compute_basis_singular_after_elimination() {
        // The rows are linearly dependent, although both are nonzero.
        // The second pivot must detect singularity after elimination.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x},{1.0_x,1.0_x}},dp);
        ARIADNE_TEST_THROWS(solver.compute_basis(A),SingularLinearProgram);
    }

    Void test_compute_basis_largest_nonunit_pivot() {
        // Neither column is a unit vector. The leading entry in
        // column zero is smaller than in column one, so the pivot
        // strategy must swap the columns instead of accepting the
        // initial order. The permuted basis is invertible (det=5).
        RationalLinearSimplexSolver solver;
        Matrix<Rational> A({{Rational(1),Rational(3)},
                            {Rational(2),Rational(1)}});

        auto basis=solver.compute_basis(A);
        const auto& permutation=basis.first;
        const auto& inverse_basis=basis.second;
        ARIADNE_TEST_EQUAL(permutation.size(),2u);
        ARIADNE_TEST_EQUAL(permutation[0],1u);
        ARIADNE_TEST_EQUAL(permutation[1],0u);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,permutation,inverse_basis));
    }

    Void test_compute_basis_rational_nonunit() {
        // The rational instantiation must compute a full inverse basis
        // even when neither column is a unit column. det(A)=3.
        // This exercises fraction-preserving elimination separately
        // from the floating-point simplex instantiation.
        RationalLinearSimplexSolver solver;
        Matrix<Rational> A({{Rational(2),Rational(1)},
                            {Rational(1),Rational(2)}});

        auto basis=solver.compute_basis(A);
        const auto& p=basis.first;
        const auto& B=basis.second;
        ARIADNE_TEST_EQUAL(p.size(),2u);
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],1u);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,p,B));
    }

    Void test_compute_basis_nonunit_two_by_two() {
        // With no unit column, both pivot stages require elimination.
        // The nonzero upper off-diagonal also requires backsubstitution.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{2.0_x,1.0_x},{1.0_x,2.0_x}},dp);
        Array<SizeType> p(2u);
        Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(2u,dp);

        make_lpair(p,B)=solver.compute_basis(A);
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],1u);
        ARIADNE_TEST_ASSERT(definitely(B[0][0]>FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[1][1]>FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[0][1]<FloatDPBounds(-0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(B[1][0]<FloatDPBounds(-0.25_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,p,B));
    }

    Void test_hotstarted_feasibility_optional_basis() {
        // The full hotstart interface can construct a basis when no
        // variable types are supplied, or reconstruct the permutation
        // from caller-supplied variable types.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);

        Array<Slackness> inferred_vt(0u);
        Array<SizeType> inferred_p(0u);
        Matrix<FloatDPBounds> inferred_B=Matrix<FloatDPBounds>::identity(1u,dp);
        Vector<FloatDPBounds> inferred_x({0.0_x,0.0_x},dp);
        Vector<FloatDPBounds> inferred_y({0.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.hotstarted_feasible(
            xl,xu,A,b,inferred_vt,inferred_p,inferred_B,inferred_x,inferred_y)));
        ARIADNE_TEST_EQUAL(inferred_vt.size(),2u);
        ARIADNE_TEST_EQUAL(inferred_vt[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(inferred_vt[1],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(inferred_p.size(),2u);
        ARIADNE_TEST_EQUAL(inferred_p[0],0u);
        ARIADNE_TEST_EQUAL(inferred_p[1],1u);
        ARIADNE_TEST_ASSERT(definitely(inferred_x[0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(inferred_x[1]==FloatDPBounds(0.0_x,dp)));

        Array<Slackness> supplied_vt(2u);
        supplied_vt[0]=Slackness::BASIS;
        supplied_vt[1]=Slackness::LOWER;
        Array<SizeType> reconstructed_p(0u);
        Matrix<FloatDPBounds> reconstructed_B=Matrix<FloatDPBounds>::identity(1u,dp);
        Vector<FloatDPBounds> reconstructed_x({0.0_x,0.0_x},dp);
        Vector<FloatDPBounds> reconstructed_y({0.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.hotstarted_feasible(
            xl,xu,A,b,supplied_vt,reconstructed_p,reconstructed_B,reconstructed_x,reconstructed_y)));
        ARIADNE_TEST_EQUAL(reconstructed_p.size(),2u);
        ARIADNE_TEST_EQUAL(reconstructed_p[0],0u);
        ARIADNE_TEST_EQUAL(reconstructed_p[1],1u);
        ARIADNE_TEST_ASSERT(definitely(reconstructed_x[0]==inferred_x[0]));
        ARIADNE_TEST_ASSERT(definitely(reconstructed_x[1]==inferred_x[1]));
        ARIADNE_TEST_ASSERT(definitely(reconstructed_B[0][0]==inferred_B[0][0]));
    }

    Void test_consistency_check_variable_types() {
        // The permutation may list only basic variables or all variables;
        // in either case the result is the number of basic variables.
        LinearSimplexSolver<FloatDP> solver;
        Array<Slackness> types(3u);
        types[0]=Slackness::BASIS;
        types[1]=Slackness::LOWER;
        types[2]=Slackness::BASIS;

        Array<SizeType> basic_only(2u);
        basic_only[0]=2u;
        basic_only[1]=0u;
        ARIADNE_TEST_EQUAL(solver.consistency_check(types,basic_only),2u);

        Array<SizeType> complete(3u);
        complete[0]=2u;
        complete[1]=0u;
        complete[2]=1u;
        ARIADNE_TEST_EQUAL(solver.consistency_check(types,complete),2u);

        Array<Slackness> all_basic(2u);
        all_basic[0]=Slackness::BASIS;
        all_basic[1]=Slackness::BASIS;
        Array<SizeType> all_basic_order(2u);
        all_basic_order[0]=1u;
        all_basic_order[1]=0u;
        ARIADNE_TEST_EQUAL(solver.consistency_check(all_basic,all_basic_order),2u);

        Array<Slackness> no_basic(2u);
        no_basic[0]=Slackness::LOWER;
        no_basic[1]=Slackness::UPPER;
        Array<SizeType> no_basic_order(2u);
        no_basic_order[0]=1u;
        no_basic_order[1]=0u;
        ARIADNE_TEST_EQUAL(solver.consistency_check(no_basic,no_basic_order),0u);
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

    Void test_invalid_partial_basis_permutations() {
        // The optional partial permutation lists m basic variable indices.
        // An out-of-range index and a repeated index are distinct caller
        // errors; neither is allowed to reach basis construction.
        LinearSimplexSolver<FloatDP> solver;

        Matrix<FloatDP> A_one({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b_one({1.0_x},dp);
        Vector<FloatDP> c_one({1.0_x,2.0_x},dp);
        Vector<FloatDP> lower_one({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper_one({1.0_x,1.0_x},dp);
        Array<Slackness> types_one(2u);
        types_one[0]=Slackness::BASIS;
        types_one[1]=Slackness::LOWER;
        Array<SizeType> out_of_range(1u);
        out_of_range[0]=2u; // Valid column indices are 0 and 1.
        Matrix<FloatDPBounds> inverse_one=Matrix<FloatDPBounds>::identity(1u,dp);
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(
            c_one,lower_one,upper_one,A_one,b_one,
            types_one,out_of_range,inverse_one));

        Matrix<FloatDP> A_two({{1.0_x,0.0_x,0.0_x},
                               {0.0_x,1.0_x,0.0_x}},dp);
        Vector<FloatDP> b_two({0.5_x,0.5_x},dp);
        Vector<FloatDP> c_two({1.0_x,2.0_x,3.0_x},dp);
        Vector<FloatDP> lower_two({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDP> upper_two({1.0_x,1.0_x,1.0_x},dp);
        Array<Slackness> types_two(3u);
        types_two[0]=Slackness::BASIS;
        types_two[1]=Slackness::BASIS;
        types_two[2]=Slackness::LOWER;
        Array<SizeType> duplicate(2u);
        duplicate[0]=0u;
        duplicate[1]=0u;
        Matrix<FloatDPBounds> inverse_two=Matrix<FloatDPBounds>::identity(2u,dp);
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(
            c_two,lower_two,upper_two,A_two,b_two,
            types_two,duplicate,inverse_two));
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

    Void test_verify_feasibility_basic_statuses() {
        // With no nonbasic variables, Ax=b uniquely determines x.
        // Points on a bound are indeterminate for strict feasibility;
        // points strictly outside the box are certainly infeasible.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x}},dp);
        Vector<FloatDP> xl({0.0_x},dp);
        Vector<FloatDP> xu({1.0_x},dp);
        Array<Slackness> vt(1u);
        vt[0]=Slackness::BASIS;

        ARIADNE_TEST_ASSERT(definitely(solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({0.5_x},dp),vt)));
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({0.0_x},dp),vt)));
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({1.0_x},dp),vt)));
        ARIADNE_TEST_ASSERT(definitely(not solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({-1.0_x},dp),vt)));
        ARIADNE_TEST_ASSERT(definitely(not solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({2.0_x},dp),vt)));
    }

    Void test_verify_feasibility_dual_signs() {
        // For x0+x1=b, x0 is basic in [0,1] and nonbasic x1 in [-1,1].
        // If a basic bound is violated, the reduced-cost sign determines
        // whether the current basis certifies infeasibility or permits a pivot.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> xl({0.0_x,-1.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;

        // x1 at LOWER: b=-2 implies x0=-1 and no feasible solution.
        vt[1]=Slackness::LOWER;
        ARIADNE_TEST_ASSERT(definitely(not solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({-2.0_x},dp),vt)));
        // x0=2 at this basis, but x1 can increase to make x0 feasible.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({1.0_x},dp),vt)));

        // x1 at UPPER: b=3 implies x0=2 and no feasible solution.
        vt[1]=Slackness::UPPER;
        ARIADNE_TEST_ASSERT(definitely(not solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({3.0_x},dp),vt)));
        // x0=-1 at this basis, but x1 can decrease to make x0 feasible.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.verify_feasibility(
            xl,xu,A,Vector<FloatDP>({0.0_x},dp),vt)));
    }

    Void test_approximate_upper_entry_basis_exchange() {
        // The generic ratio test must handle a nonbasic variable entering
        // from its upper bound. For x0+x1=1.5, moving x1 down from 1
        // makes the basic variable x0 reach its upper bound at 1, causing
        // x1 to enter the basis. This also exercises a finite leaving ratio.
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDPApproximation> rhs({1.5_x},dp);
        Vector<FloatDPApproximation> lower({0.0_x,-inf},dp);
        Vector<FloatDPApproximation> upper({1.0_x,1.0_x},dp);

        Array<Slackness> variable_types(2u);
        variable_types[0]=Slackness::BASIS;
        variable_types[1]=Slackness::UPPER;
        auto basis=solver.compute_basis(A);
        auto permutation=basis.first;
        auto inverse_basis=basis.second;
        auto primal=solver.compute_x(lower,upper,A,rhs,variable_types);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,variable_types,permutation,inverse_basis,primal));

        SizeType leaving_index=solver.lpstep(
            lower,upper,A,rhs,variable_types,permutation,inverse_basis,primal,1u);
        ARIADNE_TEST_EQUAL(leaving_index,0u);
        ARIADNE_TEST_EQUAL(variable_types[0],Slackness::UPPER);
        ARIADNE_TEST_EQUAL(variable_types[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(permutation[0],1u);
        ARIADNE_TEST_EQUAL(permutation[1],0u);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,variable_types,permutation,inverse_basis,primal));
    }

    Void test_upper_only_nonbasic_bounds() {
        // A variable bounded only from above must initialise as UPPER.
        // A negative step from that bound enters x1 and makes x0 hit
        // its upper bound, rather than performing an unbounded flip.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.5_x},dp);
        Vector<FloatDP> xl({0.0_x,-inf},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.feasible(xl,xu,A,b)));
        ARIADNE_TEST_ASSERT(definitely(solver.validated_feasible(xl,xu,A,b)));

        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::UPPER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;
        Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(1u,dp);
        Vector<FloatDPBounds> x=solver.compute_x(xl,xu,A,b,vt);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));

        SizeType r=solver.lpstep(xl,xu,A,b,vt,p,B,x,1u);
        ARIADNE_TEST_EQUAL(r,0u);
        ARIADNE_TEST_EQUAL(vt[0],Slackness::UPPER);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(p[1],0u);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));
    }

    Void test_unbounded_nonbasic_initialisation() {
        // compute_vt cannot choose LOWER or UPPER for a nonbasic
        // variable whose two bounds are both infinite.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        Vector<FloatDP> xl({0.0_x,-inf},dp);
        Vector<FloatDP> xu({1.0_x,+inf},dp);

        ARIADNE_TEST_THROWS(solver.feasible(xl,xu,A,b),
                            DegenerateFeasibilityProblemException);
        ARIADNE_TEST_THROWS(solver.validated_feasible(xl,xu,A,b),
                            DegenerateFeasibilityProblemException);
    }

    Void test_approximate_lpstep_without_leaving_variable() {
        // This is the approximate-arithmetic counterpart of the FloatDP
        // regression. Here the generic compute_rt overload is used.
        // With A=(1,0), the nonbasic x1 never affects Ax=b; when
        // it is unbounded in the entering direction, no basis pivot
        // or bound flip can occur.
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x,0.0_x}},dp);
        Vector<FloatDPApproximation> b({0.5_x},dp);

        for(SizeType order=0u; order!=2u; ++order) {
            Vector<FloatDPApproximation> xl((order==0u)
                ? Vector<FloatDPApproximation>({0.0_x,0.0_x},dp)
                : Vector<FloatDPApproximation>({0.0_x,-inf},dp));
            Vector<FloatDPApproximation> xu((order==0u)
                ? Vector<FloatDPApproximation>({1.0_x,+inf},dp)
                : Vector<FloatDPApproximation>({1.0_x,0.0_x},dp));
            Array<Slackness> variable_types(2u);
            variable_types[0]=Slackness::BASIS;
            variable_types[1]=(order==0u ? Slackness::LOWER : Slackness::UPPER);

            auto basis=solver.compute_basis(A);
            auto permutation=basis.first;
            auto inverse_basis=basis.second;
            auto primal=solver.compute_x(xl,xu,A,b,variable_types);
            ARIADNE_TEST_EXECUTE(solver.consistency_check(
                xl,xu,A,b,variable_types,permutation,inverse_basis,primal));

            SizeType leaving_index=solver.lpstep(
                xl,xu,A,b,variable_types,permutation,inverse_basis,primal,1u);
            ARIADNE_TEST_EQUAL(leaving_index,2u);
            ARIADNE_TEST_EQUAL(variable_types[0],Slackness::BASIS);
            ARIADNE_TEST_EQUAL(variable_types[1],
                (order==0u ? Slackness::LOWER : Slackness::UPPER));
            ARIADNE_TEST_EQUAL(permutation[0],0u);
            ARIADNE_TEST_EQUAL(permutation[1],1u);
            ARIADNE_TEST_EXECUTE(solver.consistency_check(
                xl,xu,A,b,variable_types,permutation,inverse_basis,primal));
        }
    }

    Void test_lpstep_without_leaving_variable() {
        // x1 does not occur in Ax=b and can move indefinitely away
        // from its current bound. No basic variable can leave.
        // Cover both positive (LOWER) and negative (UPPER) entering steps.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        for(SizeType order=0u; order!=2u; ++order) {
            Vector<FloatDP> xl((order==0u)
                ? Vector<FloatDP>({0.0_x,0.0_x},dp)
                : Vector<FloatDP>({0.0_x,-inf},dp));
            Vector<FloatDP> xu((order==0u)
                ? Vector<FloatDP>({1.0_x,+inf},dp)
                : Vector<FloatDP>({1.0_x,0.0_x},dp));
            Array<Slackness> vt(2u);
            vt[0]=Slackness::BASIS;
            vt[1]=(order==0u ? Slackness::LOWER : Slackness::UPPER);
            Array<SizeType> p(2u);
            p[0]=0u;
            p[1]=1u;
            Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(1u,dp);
            Vector<FloatDPBounds> x=solver.compute_x(xl,xu,A,b,vt);

            ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.0_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));

            SizeType r=solver.lpstep(xl,xu,A,b,vt,p,B,x,1u);
            ARIADNE_TEST_EQUAL(r,2u);
            ARIADNE_TEST_EQUAL(vt[0],Slackness::BASIS);
            ARIADNE_TEST_EQUAL(vt[1],(order==0u ? Slackness::LOWER : Slackness::UPPER));
            ARIADNE_TEST_EQUAL(p[0],0u);
            ARIADNE_TEST_EQUAL(p[1],1u);
            ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.0_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));
        }
    }

    Void test_single_variable_infeasibility() {
        // Ax=b has exactly one solution. When it lies beyond either
        // bound, there is no nonbasic variable to enter the basis,
        // so both algorithms must certify infeasibility.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x}},dp);
        Vector<FloatDP> xl({0.0_x},dp);
        Vector<FloatDP> xu({1.0_x},dp);
        Vector<FloatDP> above({2.0_x},dp);
        Vector<FloatDP> below({-1.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(not solver.feasible(xl,xu,A,above)));
        ARIADNE_TEST_ASSERT(definitely(not solver.validated_feasible(xl,xu,A,above)));
        ARIADNE_TEST_ASSERT(definitely(not solver.feasible(xl,xu,A,below)));
        ARIADNE_TEST_ASSERT(definitely(not solver.validated_feasible(xl,xu,A,below)));
    }

    Void test_lpstep_upper_bound_flip() {
        // x0+x1=0.5, x0 in [-5,5] and x1 in [0,1].
        // Decreasing x1 from UPPER to LOWER does not hit a bound on x0.
        // The operation must flip the nonbasic bound without a pivot.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({0.5_x},dp);
        Vector<FloatDP> xl({-5.0_x,0.0_x},dp);
        Vector<FloatDP> xu({5.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::UPPER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;
        Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(1u,dp);
        Vector<FloatDPBounds> x=solver.compute_x(xl,xu,A,b,vt);

        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(-0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));

        SizeType r=solver.lpstep(xl,xu,A,b,vt,p,B,x,1u);
        ARIADNE_TEST_EQUAL(r,1u);
        ARIADNE_TEST_EQUAL(vt[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],1u);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));
    }

    Void test_lpstep_upper_entry() {
        // Start with a feasible basis: x0=0.5, x1=1 (UPPER).
        // Decreasing x1 reaches the upper bound of x0 at 1.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.5_x},dp);
        Vector<FloatDP> xl({0.0_x,-1.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::UPPER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;
        Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(1u,dp);
        Vector<FloatDPBounds> x=solver.compute_x(xl,xu,A,b,vt);

        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));

        SizeType r=solver.lpstep(xl,xu,A,b,vt,p,B,x,1u);
        ARIADNE_TEST_EQUAL(r,0u);
        ARIADNE_TEST_EQUAL(vt[0],Slackness::UPPER);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(p[1],0u);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));
    }

    Void test_upper_entry_tied_upper_bounds() {
        // x0+x2=x1+x2=1.5, initially x=(0.5,0.5,1).
        // Entering x2 from UPPER makes both basic variables hit
        // UPPER at the same step. Bland's rule must select x0
        // regardless of the current ordering of the basis.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x,1.0_x},
                           {0.0_x,1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.5_x,1.5_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x,1.0_x},dp);

        for(SizeType order=0u; order!=2u; ++order) {
            Array<Slackness> vt(3u);
            vt[0]=Slackness::BASIS;
            vt[1]=Slackness::BASIS;
            vt[2]=Slackness::UPPER;
            Array<SizeType> p(3u);
            Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(2u,dp);
            if(order==0u) {
                p[0]=0u;
                p[1]=1u;
            } else {
                p[0]=1u;
                p[1]=0u;
                B[0][0]=FloatDPBounds(0.0_x,dp);
                B[0][1]=FloatDPBounds(1.0_x,dp);
                B[1][0]=FloatDPBounds(1.0_x,dp);
                B[1][1]=FloatDPBounds(0.0_x,dp);
            }
            p[2]=2u;
            Vector<FloatDPBounds> x=solver.compute_x(xl,xu,A,b,vt);
            ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[2]==FloatDPBounds(1.0_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));

            SizeType r=solver.lpstep(xl,xu,A,b,vt,p,B,x,2u);
            ARIADNE_TEST_EQUAL(r,(order==0u ? 0u : 1u));
            ARIADNE_TEST_EQUAL(vt[0],Slackness::UPPER);
            ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
            ARIADNE_TEST_EQUAL(vt[2],Slackness::BASIS);
            ARIADNE_TEST_EQUAL(p[0],(order==0u ? 2u : 1u));
            ARIADNE_TEST_EQUAL(p[1],(order==0u ? 1u : 2u));
            ARIADNE_TEST_EQUAL(p[2],0u);
            ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(1.0_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(1.0_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[2]==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));
        }
    }

    Void test_validated_multiple_upper_violations() {
        // x0-x2=0.5 and x1-x2=0.75. With x2 at UPPER=1,
        // both basic variables violate their upper bounds.
        // Decreasing x2 first takes x0 to LOWER and exchanges
        // it with x2, leaving a strictly feasible new basis.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x,-1.0_x},
                           {0.0_x,1.0_x,-1.0_x}},dp);
        Vector<FloatDP> b({0.5_x,0.75_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x,-5.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x,1.0_x},dp);
        Array<Slackness> vt(3u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::BASIS;
        vt[2]=Slackness::UPPER;
        Array<SizeType> p(3u);
        p[0]=0u;
        p[1]=1u;
        p[2]=2u;

        auto initial=solver.compute_x(xl,xu,A,b,vt);
        ARIADNE_TEST_ASSERT(definitely(initial[0]==FloatDPBounds(1.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(initial[1]==FloatDPBounds(1.75_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(initial[2]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_EQUAL(vt[0],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(vt[2],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(p[0],2u);
        ARIADNE_TEST_EQUAL(p[1],1u);
        ARIADNE_TEST_EQUAL(p[2],0u);
        auto x=solver.compute_x(xl,xu,A,b,vt);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[2]==FloatDPBounds(-0.5_x,dp)));
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_ASSERT(definitely(solver.verify_feasibility(xl,xu,A,b,vt)));
    }

    Void test_degenerate_ratio_tie_breaking() {
        // At x=(1,1,0), entering x2 makes both basic variables
        // reach their lower bounds at t=1. The leaving variable
        // must be x0 (smaller original index), independently of
        // the order of the basis columns in p.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x,1.0_x},
                           {0.0_x,1.0_x,1.0_x}},dp);
        Vector<FloatDP> b({1.0_x,1.0_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({2.0_x,2.0_x,2.0_x},dp);

        for(SizeType order=0u; order!=2u; ++order) {
            Array<Slackness> vt(3u);
            vt[0]=Slackness::BASIS;
            vt[1]=Slackness::BASIS;
            vt[2]=Slackness::LOWER;
            Array<SizeType> p(3u);
            Matrix<FloatDPBounds> B=Matrix<FloatDPBounds>::identity(2u,dp);
            if(order==0u) {
                p[0]=0u;
                p[1]=1u;
            } else {
                p[0]=1u;
                p[1]=0u;
                B[0][0]=FloatDPBounds(0.0_x,dp);
                B[0][1]=FloatDPBounds(1.0_x,dp);
                B[1][0]=FloatDPBounds(1.0_x,dp);
                B[1][1]=FloatDPBounds(0.0_x,dp);
            }
            p[2]=2u;
            Vector<FloatDPBounds> x=solver.compute_x(xl,xu,A,b,vt);
            ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(1.0_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(1.0_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[2]==FloatDPBounds(0.0_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));

            SizeType r=solver.lpstep(xl,xu,A,b,vt,p,B,x,2u);
            ARIADNE_TEST_EQUAL(r,(order==0u ? 0u : 1u));
            ARIADNE_TEST_EQUAL(vt[0],Slackness::LOWER);
            ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
            ARIADNE_TEST_EQUAL(vt[2],Slackness::BASIS);
            ARIADNE_TEST_EQUAL(p[0],(order==0u ? 2u : 1u));
            ARIADNE_TEST_EQUAL(p[1],(order==0u ? 1u : 2u));
            ARIADNE_TEST_EQUAL(p[2],0u);
            ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.0_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(0.0_x,dp)));
            ARIADNE_TEST_ASSERT(definitely(x[2]==FloatDPBounds(1.0_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(xl,xu,A,b,vt,p,B,x));
        }
    }

    Void test_validated_upper_entry_bound_flip() {
        // x0-x1=1.5 starts at x=(2.5,1), with x0 above its upper
        // bound and nonbasic x1 at UPPER. Moving x1 down to LOWER
        // reaches a strictly feasible basic x0 without a basis change.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,-1.0_x}},dp);
        Vector<FloatDP> b({1.5_x},dp);
        Vector<FloatDP> xl({0.0_x,-1.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::UPPER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;

        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_EQUAL(vt[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],1u);
        auto x=solver.compute_x(xl,xu,A,b,vt);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(-1.0_x,dp)));
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_ASSERT(definitely(solver.verify_feasibility(xl,xu,A,b,vt)));
    }

    Void test_validated_upper_entry_basis_exchange() {
        // Wider bounds for x1 force x0 to hit LOWER before x1
        // reaches its lower bound. The UPPER entry must exchange
        // x1 into the basis with the resulting point (0,-1.5).
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,-1.0_x}},dp);
        Vector<FloatDP> b({1.5_x},dp);
        Vector<FloatDP> xl({0.0_x,-5.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::UPPER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;

        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_EQUAL(vt[0],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(p[1],0u);
        auto x=solver.compute_x(xl,xu,A,b,vt);
        ARIADNE_TEST_ASSERT(definitely(x[0]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(x[1]==FloatDPBounds(-1.5_x,dp)));
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_ASSERT(definitely(solver.verify_feasibility(xl,xu,A,b,vt)));
    }

    Void test_validated_feasibility_stationary_at_boundary() {
        // x0=1 is on its upper bound. The nonbasic x1 has zero
        // coefficient in A, hence its reduced cost is exactly zero.
        // No strictly improving pivot exists, although this basis
        // cannot certify strict feasibility at the boundary.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x}},dp);
        Vector<FloatDP> rhs({1.0_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x},dp);
        Array<Slackness> variable_types(2u);
        variable_types[0]=Slackness::BASIS;
        variable_types[1]=Slackness::LOWER;
        Array<SizeType> permutation(2u);
        permutation[0]=0u;
        permutation[1]=1u;

        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(
            lower,upper,A,rhs,variable_types,permutation));
        ARIADNE_TEST_EQUAL(variable_types[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(variable_types[1],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(permutation[0],0u);
        ARIADNE_TEST_EQUAL(permutation[1],1u);
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.verify_feasibility(
            lower,upper,A,rhs,variable_types)));
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

    Void test_validated_feasibility_basis_exchange() {
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,2.0_x}},dp);
        Vector<FloatDP> b({1.5_x},dp);
        Vector<FloatDP> xl({0.0_x,0.0_x},dp);
        Vector<FloatDP> xu({1.0_x,1.0_x},dp);
        Array<Slackness> vt(2u);
        vt[0]=Slackness::BASIS;
        vt[1]=Slackness::LOWER;
        Array<SizeType> p(2u);
        p[0]=0u;
        p[1]=1u;
        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(xl,xu,A,b,vt,p));
        ARIADNE_TEST_EQUAL(vt[0],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(vt[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(p[1],0u);
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(xl,xu,A,b,vt,p));
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
