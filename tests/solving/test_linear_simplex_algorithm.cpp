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
        ARIADNE_TEST_CALL(test_logging_diagnostics());
        ARIADNE_TEST_CALL(test_public_simplex_input_contracts());
        ARIADNE_TEST_CALL(test_slackness_output());
        ARIADNE_TEST_CALL(test_pivot_to_better_vertex());
        ARIADNE_TEST_CALL(test_bland_rule_with_reordered_nonbasic_variables());
        ARIADNE_TEST_CALL(test_bland_rule_rejects_inconsistent_partition());
        ARIADNE_TEST_CALL(test_rounded_simplex_feasibility_and_optimum());
        ARIADNE_TEST_CALL(test_rounded_simplex_hotstarted_paths());
        ARIADNE_TEST_CALL(test_verify_feasibility_basic_statuses());
        ARIADNE_TEST_CALL(test_verify_feasibility_input_contracts());
        ARIADNE_TEST_CALL(test_verify_feasibility_dual_signs());
        ARIADNE_TEST_CALL(test_verify_feasibility_with_narrow_bounds());
        ARIADNE_TEST_CALL(test_uncertain_feasibility_at_large_bounds());
        ARIADNE_TEST_CALL(test_uncertain_feasibility_residual());
        ARIADNE_TEST_CALL(test_upper_only_nonbasic_bounds());
        ARIADNE_TEST_CALL(test_approximate_upper_entry_basis_exchange());
        ARIADNE_TEST_CALL(test_unbounded_nonbasic_initialisation());
        ARIADNE_TEST_CALL(test_lpstep_without_leaving_variable());
        ARIADNE_TEST_CALL(test_lpstep_entering_index_contracts());
        ARIADNE_TEST_CALL(test_approximate_lpstep_without_leaving_variable());
        ARIADNE_TEST_CALL(test_single_variable_infeasibility());
        ARIADNE_TEST_CALL(test_approximate_infeasible_equality());
        ARIADNE_TEST_CALL(test_feasibility_iteration_budget());
        ARIADNE_TEST_CALL(test_optimum_on_last_permitted_iteration());
        ARIADNE_TEST_CALL(test_optimisation_iteration_budget());
        ARIADNE_TEST_CALL(test_lpstep_upper_bound_flip());
        ARIADNE_TEST_CALL(test_lpstep_upper_entry());
        ARIADNE_TEST_CALL(test_upper_entry_tied_upper_bounds());
        ARIADNE_TEST_CALL(test_validated_multiple_upper_violations());
        ARIADNE_TEST_CALL(test_degenerate_ratio_tie_breaking());
        ARIADNE_TEST_CALL(test_zero_step_lower_bound_pivot());
        ARIADNE_TEST_CALL(test_zero_step_upper_bound_pivot());
        ARIADNE_TEST_CALL(test_rounding_at_lower_basic_bound());
        ARIADNE_TEST_CALL(test_rounding_at_upper_basic_bound());
        ARIADNE_TEST_CALL(test_rigorous_warm_start_lower_correction());
        ARIADNE_TEST_CALL(test_rigorous_warm_start_upper_correction());
        ARIADNE_TEST_CALL(test_rigorous_warm_start_lower_failure());
        ARIADNE_TEST_CALL(test_rigorous_warm_start_upper_failure());
        ARIADNE_TEST_CALL(test_validated_upper_entry_bound_flip());
        ARIADNE_TEST_CALL(test_validated_upper_entry_basis_exchange());
        ARIADNE_TEST_CALL(test_validated_strict_feasibility());
        ARIADNE_TEST_CALL(test_validated_reduced_cost_certificates());
        ARIADNE_TEST_CALL(test_validated_ratio_skips_violated_lower_bound());
        ARIADNE_TEST_CALL(test_validated_ratio_skips_violated_upper_bound());
        ARIADNE_TEST_CALL(test_validated_feasibility_stationary_at_boundary());
        ARIADNE_TEST_CALL(test_validated_feasibility_without_finite_leaving_bound());
        ARIADNE_TEST_CALL(test_validated_feasibility_bound_flip());
        ARIADNE_TEST_CALL(test_validated_feasibility_basis_exchange());
        ARIADNE_TEST_CALL(test_validated_feasibility_exchange_to_upper_bound());
        ARIADNE_TEST_CALL(test_hotstarted_basis_initialisation());
        ARIADNE_TEST_CALL(test_hotstarted_optimisation_from_alternate_basis());
        ARIADNE_TEST_CALL(test_hotstarted_minimise_partial_permutation());
        ARIADNE_TEST_CALL(test_invalid_partial_basis_permutations());
        ARIADNE_TEST_CALL(test_consistency_check_variable_types());
        ARIADNE_TEST_CALL(test_consistency_check_nonsquare_system());
        ARIADNE_TEST_CALL(test_consistency_check_rejects_inconsistent_states());
        ARIADNE_TEST_CALL(test_point_reconstruction_input_contracts());
        ARIADNE_TEST_CALL(test_hotstarted_feasibility_optional_basis());
        ARIADNE_TEST_CALL(test_warm_start_validation_corrects_false_negative());
        ARIADNE_TEST_CALL(test_warm_start_validation_corrects_false_positive());
        ARIADNE_TEST_CALL(test_compute_basis_tall_matrix_precondition());
        ARIADNE_TEST_CALL(test_compute_basis_unit_columns());
        ARIADNE_TEST_CALL(test_compute_basis_singular_after_elimination());
        ARIADNE_TEST_CALL(test_compute_basis_largest_nonunit_pivot());
        ARIADNE_TEST_CALL(test_compute_basis_rational_nonunit());
        ARIADNE_TEST_CALL(test_compute_basis_nonunit_two_by_two());
        ARIADNE_TEST_CALL(test_compute_basis_column_pivoting());
    }


    Void test_logging_diagnostics() {
        // Recheck the same simplex contracts at high verbosity: pivot
        // selection, bound flips, basis exchange and optimal termination
        // must not depend on whether diagnostic streams are enabled.
        auto& configuration=Logger::instance().configuration();
        const unsigned int old_verbosity=configuration.verbosity();
        struct RestoreVerbosity {
            unsigned int value;
            ~RestoreVerbosity() {
                Logger::instance().configuration().set_verbosity(value);
            }
        } restore{old_verbosity};
        configuration.set_verbosity(10u);

        test_pivot_to_better_vertex();
        test_bland_rule_with_reordered_nonbasic_variables();
        test_rounded_simplex_feasibility_and_optimum();
        test_lpstep_upper_bound_flip();
        test_lpstep_upper_entry();
        test_approximate_upper_entry_basis_exchange();
        test_validated_upper_entry_basis_exchange();
        test_validated_feasibility_basis_exchange();
        test_lpstep_without_leaving_variable();
        test_validated_feasibility_without_finite_leaving_bound();
        test_validated_feasibility_stationary_at_boundary();
        test_single_variable_infeasibility();
        test_approximate_infeasible_equality();
        test_hotstarted_basis_initialisation();
        test_upper_only_nonbasic_bounds();
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

    Void test_warm_start_validation_corrects_false_negative() {
        check_warm_start_validation(true);
    }

    Void test_warm_start_validation_corrects_false_positive() {
        check_warm_start_validation(false);
    }

    Void check_warm_start_validation(Bool actually_feasible) {
        // An accepted small error in the warm-start inverse puts the
        // numerical x on the opposite side of a nearby upper bound from
        // the exact solution x=1. A determinate validation must prevail,
        // and the final primal/dual state must use the rebuilt inverse.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x}},dp);
        Vector<FloatDP> rhs({1.0_x},dp);
        Vector<FloatDP> lower({0.0_x},dp);
        Vector<FloatDP> upper({actually_feasible
            ? 1.0000000000004547_pr : 0.9999999999995453_pr},dp);
        Array<Slackness> types({Slackness::BASIS});
        Array<SizeType> p({0u});
        auto inverse=Matrix<FloatDPBounds>::identity(1u,dp);
        FloatDPBounds error=FloatDPBounds(1,dp)/FloatDPBounds(1099511627776.0_x,dp);
        inverse[0][0]=actually_feasible ? FloatDPBounds(1,dp)+error : FloatDPBounds(1,dp)-error;
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,p,inverse));
        Vector<FloatDPBounds> primal(1u,dp),dual(1u,dp);
        auto status=solver.hotstarted_feasible(lower,upper,A,rhs,types,p,inverse,primal,dual);
        ARIADNE_TEST_ASSERT(definitely(actually_feasible ? status : not status));
        ARIADNE_TEST_ASSERT(definitely(inverse[0][0]==FloatDPBounds(1,dp)));
        ARIADNE_TEST_ASSERT(definitely(primal[0]==FloatDPBounds(1,dp)));
        ARIADNE_TEST_ASSERT(definitely(dual[0]==FloatDPBounds(actually_feasible ? 0 : 1,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,rhs,primal));
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

    Void test_consistency_check_rejects_inconsistent_states() {
        // These public diagnostics accept caller-supplied basis and primal
        // data. Reject each inconsistent component independently.
        LinearSimplexSolver<FloatDP> solver;
        Array<Slackness> types({Slackness::BASIS,Slackness::LOWER});
        Array<SizeType> p({0u,1u});
        Array<SizeType> too_long({0u,1u,0u});
        Array<SizeType> bad_prefix({2u,1u});
        Array<SizeType> bad_suffix({0u,2u});
        Array<SizeType> partial_outside({2u});
        Array<SizeType> partial_nonbasic({1u});
        ARIADNE_TEST_FAIL(solver.consistency_check(types,too_long));
        ARIADNE_TEST_FAIL(solver.consistency_check(types,bad_prefix));
        ARIADNE_TEST_FAIL(solver.consistency_check(types,bad_suffix));
        ARIADNE_TEST_FAIL(solver.consistency_check(types,partial_outside));
        ARIADNE_TEST_FAIL(solver.consistency_check(types,partial_nonbasic));
        Array<Slackness> separated_basis({Slackness::BASIS,Slackness::LOWER,Slackness::BASIS});
        Array<SizeType> wrong_order({0u,1u,2u});
        ARIADNE_TEST_FAIL(solver.consistency_check(separated_basis,wrong_order));
        Array<Slackness> one_basis({Slackness::BASIS,Slackness::LOWER,Slackness::LOWER});
        Array<SizeType> outside_suffix({0u,1u,3u});
        ARIADNE_TEST_FAIL(solver.consistency_check(one_basis,outside_suffix));

        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> rhs({0.5_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x},dp);
        Matrix<FloatDPBounds> inverse({{1.0_x}},dp);
        Matrix<FloatDPBounds> wrong_inverse({{2.0_x}},dp);
        Vector<FloatDPBounds> primal({0.5_x,0.0_x},dp);
        Vector<FloatDPBounds> wrong_bound({0.25_x,0.25_x},dp);
        Vector<FloatDPBounds> wrong_residual({0.75_x,0.0_x},dp);
        ARIADNE_TEST_FAIL(solver.consistency_check(A,p,wrong_inverse));
        ARIADNE_TEST_FAIL(solver.consistency_check(A,rhs,wrong_residual));
        ARIADNE_TEST_FAIL(solver.consistency_check(
            lower,upper,A,rhs,types,p,wrong_inverse,primal));
        ARIADNE_TEST_FAIL(solver.consistency_check(
            lower,upper,A,rhs,types,p,inverse,wrong_bound));
        ARIADNE_TEST_FAIL(solver.consistency_check(
            lower,upper,A,rhs,types,p,inverse,wrong_residual));
        Array<Slackness> invalid_nonbasic({Slackness::BASIS,Slackness::FIXED});
        ARIADNE_TEST_FAIL(solver.consistency_check(
            lower,upper,A,rhs,invalid_nonbasic,p,inverse,primal));
    }

    Void test_point_reconstruction_input_contracts() {
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x,0.0_x}},dp);
        Vector<FloatDPApproximation> rhs({0.5_x},dp);
        Vector<FloatDPApproximation> lower({0.0_x,0.0_x},dp);
        Vector<FloatDPApproximation> upper({1.0_x,1.0_x},dp);
        Array<Slackness> incomplete_types({Slackness::BASIS});
        ARIADNE_TEST_FAIL(solver.compute_x(lower,upper,A,rhs,incomplete_types));

        // A warm-start inverse must reconstruct Ax=b to the stated
        // tolerance. This small relative error in B passes the inverse
        // check, but produces an unacceptable residual for a large rhs.
        Matrix<FloatDPApproximation> one_variable({{1.0_x}},dp);
        Vector<FloatDPApproximation> large_rhs({1073741824.0_x},dp);
        Vector<FloatDPApproximation> large_lower({0.0_x},dp);
        Vector<FloatDPApproximation> large_upper({2147483648.0_x},dp);
        Vector<FloatDPApproximation> costs({0.0_x},dp);
        Array<Slackness> types({Slackness::BASIS});
        Array<SizeType> p({0u});
        Matrix<FloatDPApproximation> inaccurate_inverse({{1.0000000000009094947017729282379150390625_x}},dp);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(one_variable,p,inaccurate_inverse));
        ARIADNE_TEST_FAIL(solver.hotstarted_minimise(costs,large_lower,large_upper,
            one_variable,large_rhs,types,p,inaccurate_inverse));

        // With no equations, a NaN nonbasic bound reaches the point
        // reconstruction check without any invalid indexing or arithmetic.
        Matrix<FloatDPApproximation> no_equations(0u,1u,FloatDPApproximation(0.0_x,dp));
        Vector<FloatDPApproximation> no_rhs(0u,dp);
        Vector<FloatDPApproximation> nan_lower(1u,FloatDPApproximation(FloatDP::nan(dp)));
        Vector<FloatDPApproximation> finite_upper({1.0_x},dp);
        Array<Slackness> nonbasic({Slackness::LOWER});
        ARIADNE_TEST_FAIL(solver.compute_x(nan_lower,finite_upper,no_equations,no_rhs,nonbasic));
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

    Void test_verify_feasibility_input_contracts() {
        // Validate each independent public input contract before basis
        // reconstruction. The underlying feasible system is A=(1,1),
        // rhs=1/2, with x0 basic and x1 at its lower bound.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> rhs({0.5_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x},dp);
        Vector<FloatDP> short_bounds({0.0_x},dp);
        Vector<FloatDP> wrong_rhs({0.5_x,0.5_x},dp);
        Array<Slackness> types(2u);
        types[0]=Slackness::BASIS;
        types[1]=Slackness::LOWER;
        Array<Slackness> wrong_types(1u);
        wrong_types[0]=Slackness::BASIS;

        ARIADNE_TEST_FAIL(solver.verify_feasibility(
            lower,upper,A,wrong_rhs,types));
        ARIADNE_TEST_FAIL(solver.verify_feasibility(
            short_bounds,upper,A,rhs,types));
        ARIADNE_TEST_FAIL(solver.verify_feasibility(
            lower,short_bounds,A,rhs,types));
        ARIADNE_TEST_FAIL(solver.verify_feasibility(
            lower,upper,A,rhs,wrong_types));

        // Singleton bounds are valid constraints, but the fixed
        // variable cannot be selected as a basic variable.
        Matrix<FloatDP> one_variable_A({{1.0_x}},dp);
        Vector<FloatDP> fixed_value({0.0_x},dp);
        Array<Slackness> fixed_basic(1u);
        fixed_basic[0]=Slackness::BASIS;
        ARIADNE_TEST_FAIL(solver.verify_feasibility(
            fixed_value,fixed_value,one_variable_A,fixed_value,fixed_basic));

        // A singleton is permitted for a nonbasic variable: x1=0
        // leaves the basic x0=1/2 strictly inside its own bounds.
        Vector<FloatDP> fixed_nonbasic_upper({1.0_x,0.0_x},dp);
        ARIADNE_TEST_ASSERT(definitely(solver.verify_feasibility(
            lower,fixed_nonbasic_upper,A,rhs,types)));
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

    Void test_verify_feasibility_with_narrow_bounds() {
        // The exact solution is (1/5,1/5), strictly inside these bounds.
        // Inverting A in interval arithmetic gives a wider enclosure
        // crossing both bounds. That uncertainty must be reported without
        // throwing or claiming a certificate of feasibility/infeasibility.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{3.0_x,2.0_x},{2.0_x,3.0_x}},dp);
        Vector<FloatDP> rhs({1.0_x,1.0_x},dp);
        Vector<FloatDP> lower({0.19999999999999998_pr,0.19999999999999998_pr},dp);
        Vector<FloatDP> upper({0.20000000000000004_pr,0.20000000000000004_pr},dp);
        Array<Slackness> types({Slackness::BASIS,Slackness::BASIS});
        Rational exact_solution=Rational(1)/5;
        ARIADNE_TEST_ASSERT(Rational(lower[0])<exact_solution);
        ARIADNE_TEST_ASSERT(exact_solution<Rational(upper[0]));
        auto primal=solver.compute_x(lower,upper,A,rhs,types);
        ARIADNE_TEST_ASSERT(possibly(primal[0]<lower[0]));
        ARIADNE_TEST_ASSERT(possibly(primal[0]>upper[0]));
        ARIADNE_TEST_ASSERT(is_indeterminate(
            solver.verify_feasibility(lower,upper,A,rhs,types)));
    }

    Void test_uncertain_feasibility_at_large_bounds() {
        // 3*x=300000 has the exact solution x=100000. Its interval
        // enclosure straddles a bound by more than the absolute tolerance.
        // Neither endpoint placement warrants an assertion failure or a
        // definite claim when the numerical bounds cannot establish it.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{3.0_x}},dp);
        Vector<FloatDP> rhs({300000.0_x},dp);
        Vector<FloatDP> at_bound({100000.0_x},dp);
        Vector<FloatDP> below({0.0_x},dp);
        Vector<FloatDP> above({200000.0_x},dp);
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible(at_bound,above,A,rhs)));
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.feasible(below,at_bound,A,rhs)));
    }

    Void test_uncertain_feasibility_residual() {
        // The near-singular but invertible basis has an exact binary
        // inverse. Applying it to the decimal rhs loses accuracy through
        // cancellation: the residual reaches the strict tolerance.
        const auto old_rounding=FloatDP::get_rounding_mode();
        struct RestoreRounding {
            decltype(FloatDP::get_rounding_mode()) value;
            ~RestoreRounding() { FloatDP::set_rounding_mode(value); }
        } restore{old_rounding};
        FloatDP::set_rounding_to_nearest();
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x,1.0_x},
            {1.0_x,1.0000000000000002220446049250313080847263336181640625_x}},dp);
        Vector<FloatDPApproximation> rhs({0.0001_pr,0.001_pr},dp);
        Vector<FloatDPApproximation> lower({-1e14_pr,-1e14_pr},dp);
        Vector<FloatDPApproximation> upper({1e14_pr,1e14_pr},dp);
        Array<Slackness> types({Slackness::BASIS,Slackness::BASIS});
        ARIADNE_TEST_ASSERT(is_indeterminate(
            solver.hotstarted_feasible(lower,upper,A,rhs,types)));
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

    Void test_lpstep_entering_index_contracts() {
        // lpstep(..., s) is a public overload: s must index a nonbasic
        // variable in p. In particular, n is a sentinel, not an index.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> rhs({0.5_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x},dp);
        Array<Slackness> statuses(2u);
        statuses[0]=Slackness::BASIS;
        statuses[1]=Slackness::LOWER;
        Array<SizeType> permutation(2u);
        permutation[0]=0u;
        permutation[1]=1u;
        Matrix<FloatDPBounds> inverse_basis=Matrix<FloatDPBounds>::identity(1u,dp);
        auto primal=solver.compute_x(lower,upper,A,rhs,statuses);

        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,statuses,permutation,inverse_basis,primal));
        ARIADNE_TEST_FAIL(solver.lpstep(
            lower,upper,A,rhs,statuses,permutation,inverse_basis,primal,2u));
        ARIADNE_TEST_FAIL(solver.lpstep(
            lower,upper,A,rhs,statuses,permutation,inverse_basis,primal,0u));
        Array<Slackness> wrong_types=statuses;
        wrong_types[1]=Slackness::BASIS;
        ARIADNE_TEST_FAIL(solver.lpstep(
            lower,upper,A,rhs,wrong_types,permutation,inverse_basis,primal,1u));
        Vector<FloatDPBounds> short_primal({0.5_x},dp);
        ARIADNE_TEST_FAIL(solver.lpstep(
            lower,upper,A,rhs,statuses,permutation,inverse_basis,short_primal,1u));
        // Rejection must leave the valid starting basis intact.
        ARIADNE_TEST_EQUAL(statuses[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(statuses[1],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(permutation[0],0u);
        ARIADNE_TEST_EQUAL(permutation[1],1u);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,statuses,permutation,inverse_basis,primal));
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

    Void test_approximate_infeasible_equality() {
        // The unique solution of x=-1 is outside [0,1]. Exercise the
        // generic approximate-arithmetic phase and its diagnostics. The
        // separate interval-style validator cannot certify approximate
        // comparisons, so its own result remains indeterminate.
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x}},dp);
        Vector<FloatDPApproximation> rhs({-1.0_x},dp);
        Vector<FloatDPApproximation> lower({0.0_x},dp);
        Vector<FloatDPApproximation> upper({1.0_x},dp);
        Array<Slackness> types({Slackness::BASIS});
        ARIADNE_TEST_ASSERT(definitely(not solver.feasible(lower,upper,A,rhs)));
        ARIADNE_TEST_ASSERT(is_indeterminate(
            solver.verify_feasibility(lower,upper,A,rhs,types)));
    }

    Void test_feasibility_iteration_budget() {
        // Restoring x0>=0 in x0-sum(xj)=-1025 requires 1025 unit
        // bound flips. The fixed iteration budget must report exhaustion
        // rather than returning a feasibility certificate for x0<0.
        LinearSimplexSolver<FloatDPApproximation> solver;
        const SizeType n=1026u;
        FloatDPApproximation zero_value(0.0_x,dp);
        Matrix<FloatDPApproximation> A(1u,n,zero_value);
        A[0][0]=1;
        for(SizeType j=1u; j!=n; ++j) { A[0][j]=-1; }
        Vector<FloatDPApproximation> rhs({-1025.0_x},dp);
        Vector<FloatDPApproximation> lower(n,zero_value);
        Vector<FloatDPApproximation> upper(n,FloatDPApproximation(1.0_x,dp));
        ARIADNE_TEST_THROWS(solver.feasible(lower,upper,A,rhs),
                            DegenerateFeasibilityProblemException);
    }

    Void test_optimum_on_last_permitted_iteration() {
        check_optimisation_budget(1024u,true);
    }

    Void test_optimisation_iteration_budget() {
        check_optimisation_budget(1025u,false);
    }

    Void check_optimisation_budget(SizeType n, Bool completes) {
        // The equation fixes x0=1/2. Each remaining independent variable
        // needs one bound flip to minimise -sum(xj), followed by one call
        // to establish optimality. Completion on call 1024 is permitted;
        // a problem requiring another call exhausts the fixed budget.
        LinearSimplexSolver<FloatDP> solver;
        FloatDP zero_value(0.0_x,dp);
        Matrix<FloatDP> A(1u,n,zero_value);
        A[0][0]=1;
        Vector<FloatDP> rhs({0.5_x},dp);
        Vector<FloatDP> lower(n,zero_value);
        Vector<FloatDP> upper(n,FloatDP(1.0_x,dp));
        Vector<FloatDP> costs(n,FloatDP(-1.0_x,dp));
        costs[0]=0;
        if(completes) {
            auto optimum=solver.minimise(costs,lower,upper,A,rhs);
            Vector<FloatDPBounds> expected(n,FloatDPBounds(1.0_x,dp));
            expected[0]=FloatDPBounds(0.5_x,dp);
            ARIADNE_TEST_ASSERT(definitely(norm(optimum-expected)==FloatDPBounds(0.0_x,dp)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(A,rhs,optimum));
        } else {
            ARIADNE_TEST_THROWS(solver.minimise(costs,lower,upper,A,rhs),
                                DegenerateFeasibilityProblemException);
        }
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

    Void test_zero_step_lower_bound_pivot() {
        check_zero_step_pivot(false);
    }

    Void test_zero_step_upper_bound_pivot() {
        check_zero_step_pivot(true);
    }

    Void check_zero_step_pivot(Bool at_upper) {
        // A basic variable already on the blocking bound leaves at t=0.
        // The basis must change without moving the feasible vertex, and
        // its new nonbasic type must describe the bound it already occupies.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,(at_upper ? -1.0_x : 1.0_x)}},dp);
        Vector<FloatDP> rhs({(at_upper ? 1.0_x : 0.0_x)},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x},dp);
        Array<Slackness> types({Slackness::BASIS,Slackness::LOWER});
        Array<SizeType> p({0u,1u});
        Matrix<FloatDPBounds> inverse=Matrix<FloatDPBounds>::identity(1u,dp);
        auto primal=solver.compute_x(lower,upper,A,rhs,types);
        const auto initial=primal;

        SizeType leaving=solver.lpstep(lower,upper,A,rhs,types,p,inverse,primal,1u);
        ARIADNE_TEST_EQUAL(leaving,0u);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(p[1],0u);
        ARIADNE_TEST_EQUAL(types[0],(at_upper ? Slackness::UPPER : Slackness::LOWER));
        ARIADNE_TEST_EQUAL(types[1],Slackness::BASIS);
        ARIADNE_TEST_ASSERT(definitely(primal[0]==initial[0]));
        ARIADNE_TEST_ASSERT(definitely(primal[1]==initial[1]));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,types,p,inverse,primal));
    }

    Void test_rounding_at_lower_basic_bound() {
        check_rounding_at_basic_bound(true);
    }

    Void test_rounding_at_upper_basic_bound() {
        check_rounding_at_basic_bound(false);
    }

    Void check_rounding_at_basic_bound(Bool enters_from_upper) {
        // The ratio test uses 3/10, whereas reconstruction uses (1/10)*3.
        // Their floating results differ by one ulp. The reconstructed
        // basic point must be clamped to its bound within the tolerance,
        // for both signs, while retaining a small equality residual.
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x,10.0_x}},dp);
        Vector<FloatDPApproximation> rhs({enters_from_upper ? -3.0_x : 3.0_x},dp);
        Vector<FloatDPApproximation> lower({enters_from_upper ? -10.0_pr : 0.0_pr,
                                            enters_from_upper ? -0.3_pr : 0.0_pr},dp);
        Vector<FloatDPApproximation> upper({enters_from_upper ? 0.0_pr : 10.0_pr,
                                            enters_from_upper ? 0.0_pr : 0.3_pr},dp);
        Array<Slackness> types({Slackness::BASIS,
            enters_from_upper ? Slackness::UPPER : Slackness::LOWER});
        auto basis=solver.compute_basis(A);
        auto p=basis.first;
        auto inverse=basis.second;
        auto primal=solver.compute_x(lower,upper,A,rhs,types);

        ARIADNE_TEST_EQUAL(solver.lpstep(lower,upper,A,rhs,types,p,inverse,primal,1u),0u);
        ARIADNE_TEST_EQUAL(p[0],1u);
        ARIADNE_TEST_EQUAL(types[0],
            (enters_from_upper ? Slackness::UPPER : Slackness::LOWER));
        ARIADNE_TEST_EQUAL(types[1],Slackness::BASIS);
        const auto blocking_bound=enters_from_upper ? lower[1] : upper[1];
        ARIADNE_TEST_ASSERT(decide(primal[1]==blocking_bound));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,types,p,inverse,primal));
    }

    Void test_rigorous_warm_start_lower_correction() {
        check_rigorous_warm_start_reconstruction(true,false);
    }

    Void test_rigorous_warm_start_upper_correction() {
        check_rigorous_warm_start_reconstruction(false,false);
    }

    Void test_rigorous_warm_start_lower_failure() {
        check_rigorous_warm_start_reconstruction(true,true);
    }

    Void test_rigorous_warm_start_upper_failure() {
        check_rigorous_warm_start_reconstruction(false,true);
    }

    Void check_rigorous_warm_start_reconstruction(Bool enters_from_upper, Bool rejects) {
        // The warm-start inverse passes the public inverse-error check.
        // Its small error makes the ratio prefer a basis exchange before
        // a nearby entering bound. Recomputing the inverse then reveals
        // an overshoot: clamp within tolerance, report degeneracy beyond it.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,10.0_x}},dp);
        Vector<FloatDP> rhs({enters_from_upper ? -3.0_x : 3.0_x},dp);
        Vector<FloatDP> lower({enters_from_upper ? -10.0_x : 0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({enters_from_upper ? 0.0_x : 10.0_x,0.0_x},dp);
        FloatDP bound(rejects ? 0.299999999995_pr : 0.2999999999999_pr,dp);
        if(enters_from_upper) { lower[1]=-bound; } else { upper[1]=bound; }
        Array<Slackness> types({Slackness::BASIS,
            enters_from_upper ? Slackness::UPPER : Slackness::LOWER});
        Array<SizeType> p({0u,1u});
        auto inverse=Matrix<FloatDPBounds>::identity(1u,dp);
        // The perturbations are exactly 2^-35 and 2^-40 respectively.
        inverse[0][0]=FloatDPBounds(1,dp)+FloatDPBounds(1,dp)/
            FloatDPBounds(rejects ? 34359738368.0_x : 1099511627776.0_x,dp);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,p,inverse));
        auto primal=solver.compute_x(lower,upper,A,rhs,types);
        if(rejects) {
            ARIADNE_TEST_THROWS(solver.lpstep(lower,upper,A,rhs,types,p,inverse,primal,1u),
                                DegenerateFeasibilityProblemException);
        } else {
            ARIADNE_TEST_EQUAL(solver.lpstep(lower,upper,A,rhs,types,p,inverse,primal,1u),0u);
            ARIADNE_TEST_EQUAL(p[0],1u);
            ARIADNE_TEST_EQUAL(types[1],Slackness::BASIS);
            FloatDP expected=enters_from_upper ? -bound : bound;
            ARIADNE_TEST_ASSERT(definitely(primal[1]==FloatDPBounds(expected)));
            ARIADNE_TEST_EXECUTE(solver.consistency_check(
                lower,upper,A,rhs,types,p,inverse,primal));
        }
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

    Void test_validated_feasibility_without_finite_leaving_bound() {
        // Starting at (2,0), x0+x1=2 has a feasible point (1,1).
        // Increasing the LOWER nonbasic x1 improves the violated upper
        // bound on x0, but both bounds in this direction are infinite.
        // No variable can leave the basis in a finite simplex step.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> rhs({2.0_x},dp);
        Vector<FloatDP> lower({-inf,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,+inf},dp);
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
        // A valid improving ray rules out an infeasibility certificate
        // for the current basis; feasibility is not yet established.
        ARIADNE_TEST_ASSERT(is_indeterminate(solver.verify_feasibility(
            lower,upper,A,rhs,variable_types)));
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

    Void test_validated_reduced_cost_certificates() {
        // With x0+x1=b, a lower x1 cannot fix x0<0, and an upper
        // x1 cannot fix x0>1. The reduced costs exclude those directions.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x},dp);
        Array<Slackness> types({Slackness::BASIS,Slackness::LOWER});
        Array<SizeType> p({0u,1u});
        Vector<FloatDP> below({-0.5_x},dp);
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(lower,upper,A,below,types,p));
        ARIADNE_TEST_ASSERT(definitely(not solver.verify_feasibility(lower,upper,A,below,types)));

        types[1]=Slackness::UPPER;
        Vector<FloatDP> above({3.0_x},dp);
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(lower,upper,A,above,types,p));
        ARIADNE_TEST_ASSERT(definitely(not solver.verify_feasibility(lower,upper,A,above,types)));

        // An unrelated upper variable has zero reduced cost. It cannot
        // improve the violated equality, but supplies no strict certificate.
        Matrix<FloatDP> independent({{1.0_x,0.0_x}},dp);
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(lower,upper,independent,below,types,p));
        ARIADNE_TEST_ASSERT(is_indeterminate(
            solver.verify_feasibility(lower,upper,independent,below,types)));
        ARIADNE_TEST_EQUAL(types[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(types[1],Slackness::UPPER);
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],1u);
    }

    Void test_validated_ratio_skips_violated_lower_bound() {
        check_validated_ratio_with_violated_bound(false);
    }

    Void test_validated_ratio_skips_violated_upper_bound() {
        check_validated_ratio_with_violated_bound(true);
    }

    Void check_validated_ratio_with_violated_bound(Bool towards_upper) {
        // The first basic variable already violates the bound towards
        // which it moves. That bound cannot block a positive step. The
        // other basic variable supplies the admissible blocking ratio.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,0.0_x,towards_upper ? -1.0_x : 1.0_x},
                           {0.0_x,1.0_x,towards_upper ? -2.0_x : 2.0_x}},dp);
        Vector<FloatDP> rhs({towards_upper ? 1.5_x : -0.5_x,
                              towards_upper ? -0.5_x : 1.5_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,1.0_x,2.0_x},dp);
        Array<Slackness> types({Slackness::BASIS,Slackness::BASIS,Slackness::LOWER});
        Array<SizeType> p({0u,1u,2u});
        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(lower,upper,A,rhs,types,p));
        ARIADNE_TEST_EQUAL(p[0],0u);
        ARIADNE_TEST_EQUAL(p[1],2u);
        ARIADNE_TEST_EQUAL(p[2],1u);
        ARIADNE_TEST_EQUAL(types[1],towards_upper ? Slackness::UPPER : Slackness::LOWER);
        auto primal=solver.compute_x(lower,upper,A,rhs,types);
        ARIADNE_TEST_ASSERT(definitely(primal[0]==FloatDPBounds(towards_upper ? 2.25_x : -1.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(primal[1]==FloatDPBounds(towards_upper ? 1.0_x : 0.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(primal[2]==FloatDPBounds(0.75_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,rhs,primal));
    }

    Void test_validated_feasibility_exchange_to_upper_bound() {
        // The initial basis fixes x1 at its LOWER bound, giving
        // x0=-1.5 below its allowed interval [0,1]. Increasing x1
        // increases x0 until x0 reaches UPPER; x1 must enter the basis.
        // This tests the positive d[r]*t branch of the validated pivot.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,-1.0_x}},dp);
        Vector<FloatDP> rhs({-1.5_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({1.0_x,5.0_x},dp);
        Array<Slackness> variable_types(2u);
        variable_types[0]=Slackness::BASIS;
        variable_types[1]=Slackness::LOWER;
        Array<SizeType> permutation(2u);
        permutation[0]=0u;
        permutation[1]=1u;

        auto initial=solver.compute_x(lower,upper,A,rhs,variable_types);
        ARIADNE_TEST_ASSERT(definitely(initial[0]==FloatDPBounds(-1.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(initial[1]==FloatDPBounds(0.0_x,dp)));

        ARIADNE_TEST_ASSERT(!solver.validated_feasibility_step(
            lower,upper,A,rhs,variable_types,permutation));
        ARIADNE_TEST_EQUAL(variable_types[0],Slackness::UPPER);
        ARIADNE_TEST_EQUAL(variable_types[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(permutation[0],1u);
        ARIADNE_TEST_EQUAL(permutation[1],0u);

        auto feasible_point=solver.compute_x(lower,upper,A,rhs,variable_types);
        ARIADNE_TEST_ASSERT(definitely(
            feasible_point[0]==FloatDPBounds(1.0_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(
            feasible_point[1]==FloatDPBounds(2.5_x,dp)));
        ARIADNE_TEST_ASSERT(solver.validated_feasibility_step(
            lower,upper,A,rhs,variable_types,permutation));
        ARIADNE_TEST_ASSERT(definitely(solver.verify_feasibility(
            lower,upper,A,rhs,variable_types)));
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

    Void test_rounded_simplex_hotstarted_paths() {
        // A strict initial basic point (1/2,0) needs no feasibility
        // pivot. A subsequent hotstarted minimisation of 2*x0+x1
        // must instead pivot to the optimum (0,1/2). In particular,
        // exercise the RoundedFloatDP hotstart instantiation rather
        // than relying on the FloatDP tests of the same interfaces.
        LinearSimplexSolver<RoundedFloatDP> solver;
        Matrix<RoundedFloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<RoundedFloatDP> rhs({0.5_x},dp);
        Vector<RoundedFloatDP> lower({0.0_x,0.0_x},dp);
        Vector<RoundedFloatDP> upper({1.0_x,1.0_x},dp);
        Vector<RoundedFloatDP> costs({2.0_x,1.0_x},dp);

        Array<Slackness> feasibility_basis(2u);
        feasibility_basis[0]=Slackness::BASIS;
        feasibility_basis[1]=Slackness::LOWER;
        ARIADNE_TEST_ASSERT(definitely(solver.hotstarted_feasible(
            lower,upper,A,rhs,feasibility_basis)));
        ARIADNE_TEST_EQUAL(feasibility_basis[0],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(feasibility_basis[1],Slackness::LOWER);

        Array<Slackness> optimisation_basis(2u);
        optimisation_basis[0]=Slackness::BASIS;
        optimisation_basis[1]=Slackness::LOWER;
        auto optimum=solver.hotstarted_minimise(
            costs,lower,upper,A,rhs,optimisation_basis);
        ARIADNE_TEST_EQUAL(optimum.size(),2u);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,rhs,optimum));
        ARIADNE_TEST_ASSERT(decide(optimum[0]<RoundedFloatDP(0.125_x,dp)));
        ARIADNE_TEST_ASSERT(decide(optimum[1]>RoundedFloatDP(0.375_x,dp)));
        ARIADNE_TEST_ASSERT(decide(optimum[1]<RoundedFloatDP(0.625_x,dp)));
    }

    Void test_rounded_simplex_feasibility_and_optimum() {
        // Cover the RoundedFloatDP instantiation on a complete LP,
        // independently of the FloatDP and Approximation instantiations.
        // On x0+x1=1 in [0,1]^2, minimising 2*x0+x1 must select
        // x0=0, x1=1 rather than the initial vertex (1,0).
        LinearSimplexSolver<RoundedFloatDP> solver;
        Matrix<RoundedFloatDP> A({{1.0_x,1.0_x}},dp);
        Vector<RoundedFloatDP> rhs({1.0_x},dp);
        Vector<RoundedFloatDP> lower({0.0_x,0.0_x},dp);
        Vector<RoundedFloatDP> upper({1.0_x,1.0_x},dp);
        Vector<RoundedFloatDP> costs({2.0_x,1.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.feasible(lower,upper,A,rhs)));
        auto optimum=solver.minimise(costs,lower,upper,A,rhs);
        ARIADNE_TEST_EQUAL(optimum.size(),2u);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(A,rhs,optimum));

        // RoundedFloatDP comparisons are ApproximateKleenean, not
        // validated comparisons. Decide them rather than using
        // definitely(), which does not certify an approximate result.
        // The unique optimum is (0,1), not the initial vertex (1,0).
        ARIADNE_TEST_ASSERT(decide(optimum[0]<RoundedFloatDP(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(decide(optimum[1]>RoundedFloatDP(0.75_x,dp)));
        ARIADNE_TEST_ASSERT(decide(optimum[0]<optimum[1]));
    }

    Void test_bland_rule_with_reordered_nonbasic_variables() {
        // Bland's rule chooses the eligible variable with the smallest
        // ORIGINAL index, even when it appears later in the permutation.
        // Both x1 and x2 reduce 2*x0+x1+x2, but p=[0,2,1] puts x2 first.
        LinearSimplexSolver<FloatDP> solver;
        Matrix<FloatDP> A({{1.0_x,1.0_x,1.0_x}},dp);
        Vector<FloatDP> rhs({1.0_x},dp);
        Vector<FloatDP> costs({2.0_x,1.0_x,1.0_x},dp);
        Vector<FloatDP> lower({0.0_x,0.0_x,0.0_x},dp);
        Vector<FloatDP> upper({2.0_x,2.0_x,2.0_x},dp);
        Array<Slackness> types(3u);
        types[0]=Slackness::BASIS;
        types[1]=Slackness::LOWER;
        types[2]=Slackness::LOWER;
        Array<SizeType> permutation(3u);
        permutation[0]=0u;
        permutation[1]=2u;
        permutation[2]=1u;
        Matrix<FloatDPBounds> inverse_basis=Matrix<FloatDPBounds>::identity(1u,dp);
        auto primal=solver.compute_x(lower,upper,A,rhs,types);
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,types,permutation,inverse_basis,primal));

        ARIADNE_TEST_ASSERT(!solver.lpstep(
            costs,lower,upper,A,rhs,types,permutation,inverse_basis,primal));
        // x1 must enter, although it is in the LAST nonbasic position.
        ARIADNE_TEST_EQUAL(permutation[0],1u);
        ARIADNE_TEST_EQUAL(permutation[1],2u);
        ARIADNE_TEST_EQUAL(permutation[2],0u);
        ARIADNE_TEST_EQUAL(types[0],Slackness::LOWER);
        ARIADNE_TEST_EQUAL(types[1],Slackness::BASIS);
        ARIADNE_TEST_EQUAL(types[2],Slackness::LOWER);
        ARIADNE_TEST_ASSERT(definitely(primal[0]<FloatDPBounds(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(primal[1]>FloatDPBounds(0.75_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(primal[2]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_EXECUTE(solver.consistency_check(
            lower,upper,A,rhs,types,permutation,inverse_basis,primal));
        ARIADNE_TEST_ASSERT(solver.lpstep(
            costs,lower,upper,A,rhs,types,permutation,inverse_basis,primal));
    }

    Void test_bland_rule_rejects_inconsistent_partition() {
        // The leading column describes a basis variable. Labelling all
        // variables LOWER violates that partition even though dimensions,
        // the permutation and inverse matrix are otherwise correct.
        LinearSimplexSolver<FloatDPApproximation> solver;
        Matrix<FloatDPApproximation> A({{1.0_x,1.0_x}},dp);
        Vector<FloatDPApproximation> rhs({0.5_x},dp);
        Vector<FloatDPApproximation> costs({0.0_x,-1.0_x},dp);
        Vector<FloatDPApproximation> lower({0.0_x,0.0_x},dp);
        Vector<FloatDPApproximation> upper({1.0_x,1.0_x},dp);
        Vector<FloatDPApproximation> primal({0.5_x,0.0_x},dp);
        Array<Slackness> types({Slackness::LOWER,Slackness::LOWER});
        Array<SizeType> p({0u,1u});
        auto inverse=Matrix<FloatDPApproximation>::identity(1u,dp);
        ARIADNE_TEST_FAIL(solver.lpstep(costs,lower,upper,A,rhs,types,p,inverse,primal));
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
