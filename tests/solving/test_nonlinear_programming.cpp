/***************************************************************************
 *            test_nonlinear_programming.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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

#include <iostream>
#include <fstream>
#include "macro/testing.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "function/formula.hpp"
#include "solving/nonlinear_programming.hpp"
#include "solving/nonlinear_interior_point_method.hpp"
#include "solving/nonlinear_penalty_function_method.hpp"
#include "geometry/box.hpp"
#include "drawing/command_line_interface.hpp"

using namespace std;
using namespace Ariadne;

class TestNonlinearProgrammingSolver
{
  private:
    std::unique_ptr<NonlinearProgrammingSolverInterface> optimiser;
    DoublePrecision pr;
  public:
    TestNonlinearProgrammingSolver(const NonlinearProgrammingSolverInterface& opt)
        : optimiser(opt.clone()) { }

    Void test() {
        ARIADNE_TEST_CALL(test_feasibility_check());
        ARIADNE_TEST_CALL(test_unconstrained_optimisation());
        ARIADNE_TEST_CALL(test_constrained_optimisation());
        ARIADNE_TEST_CALL(test_equality_constrained_optimisation());
        ARIADNE_TEST_CALL(test_linear_feasibility());
        ARIADNE_TEST_CALL(test_nonlinear_feasibility());
        ARIADNE_TEST_CALL(test_nonlinear_equality_feasibility());
    }

    Void test_unconstrained_optimisation() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveScalarMultivariateFunction x0s = sqr(x[0]);
        EffectiveScalarMultivariateFunction f(x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0]));
        ARIADNE_TEST_PRINT(f);
        EffectiveVectorMultivariateFunction g(0u,2u);
        ARIADNE_TEST_PRINT(g);
        ExactBoxType D=ExactBoxType{{-1,2},{-3,5}};
        ExactBoxType C=ExactBoxType{};
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));

        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g(x_optimal),C);
        //ExactDouble required_accuracy=1e-8_pr;
        //ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
    }

    Void test_equality_constrained_optimisation() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveScalarMultivariateFunction f=(sqr(x[0])+sqr(x[1]));
        ARIADNE_TEST_PRINT(f);
        Real a(1.5_x); Real b(0.25_x);
        EffectiveVectorMultivariateFunction g={a+x[0]+2*x[1]+b*x[0]*x[1]};
        ARIADNE_TEST_PRINT(g);
        ExactBoxType C={{0.0_x,0.0_x}};
        ExactBoxType D=ExactBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x}};
        ARIADNE_TEST_PRINT(Ariadne::make_tuple(f,D,g,C));

        ExactDouble required_accuracy=1e-7_pr;
        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_LESS(norm(g(x_optimal)),required_accuracy);
    }

    Void test_constrained_optimisation() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(3);
        EffectiveScalarMultivariateFunction x0s = sqr(x[0]);
        EffectiveScalarMultivariateFunction f = x0s*(12+x0s*(Decimal(6.3)+x0s))+6*x[1]*(x[1]-x[0])+x[2];
        ARIADNE_TEST_PRINT(f);
        //EffectiveVectorMultivariateFunction g( (x[0]-1, x[0]+x[1]*x[1], x[1]*x[1]) );
        ExactBoxType D = ExactBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x},{-3.0_x,5.0_x}};
        ARIADNE_TEST_PRINT(D);
        EffectiveVectorMultivariateFunction g = {2*x[1]+x[0], x[0]+x[1]*x[1]-Real(0.875_x)};
        ARIADNE_TEST_PRINT(g);
        ExactBoxType C = ExactBoxType{{0.0_x,inf},{0.0_x,inf}};
        ARIADNE_TEST_PRINT(C);

        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,g,C);
        ARIADNE_TEST_BINARY_PREDICATE(element,x_optimal,D);
        ARIADNE_TEST_BINARY_PREDICATE(element,g(x_optimal),C);
        //ExactDouble required_accuracy=1e-6_pr;
        //ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
    }

    Void test_mixed_constrained_optimisation() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(3);
        EffectiveScalarMultivariateFunction f(+(sqr(x[0])+sqr(x[1])+x[1]*x[2]));
        ARIADNE_TEST_PRINT(f);
        ExactBoxType D = ExactBoxType{{-1.0_x,2.0_x},{-3.0_x,5.0_x},{1.25_x,2.25_x}};
        ARIADNE_TEST_PRINT(D);
        EffectiveScalarMultivariateFunction g = x[0]*x[1]-x[0]*Real(1.25_x);
        EffectiveVectorMultivariateFunction h = {Real(1.5_x)+x[0]+2*x[1]+Real(0.25_x)*x[0]*x[1]};
        EffectiveVectorMultivariateFunction gh=join(g,h);
        ARIADNE_TEST_PRINT(gh);
        ExactBoxType C = ExactBoxType {{-1.0_x,-0.5_x},{0.0_x,0.0_x}};
        ARIADNE_TEST_PRINT(C);

        FloatDPBoundsVector x_optimal=optimiser->minimise(f,D,gh,C);
        ExactDouble required_accuracy=1e-8_pr;
        ARIADNE_TEST_LESS(norm(h(x_optimal)),required_accuracy);
    }

    Void test_linear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveVectorMultivariateFunction g=EffectiveVectorMultivariateFunction(1u, 2*x[0]+x[1]);
        ARIADNE_TEST_PRINT(g);
        ExactBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        ExactBoxType C = ExactBoxType{{-2.0_x,1.0_x}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=ExactBoxType{{1.0_x,1.5_x}};
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=ExactBoxType{{1.0_x,1.5_x},{0.5_x,1.0_x}};
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

    Void test_nonlinear_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveVectorMultivariateFunction g = {2*x[0]+x[1]+x[0]*x[1]/8};
        ARIADNE_TEST_PRINT(g);
        ExactBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        ExactBoxType C = ExactBoxType{{-2.0_x,1.0_x}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        C=ExactBoxType{{1.0_x,1.5_x}};
        ARIADNE_TEST_ASSERT(optimiser->feasible(D,g,C));
        D=ExactBoxType{{1.0_x,1.5_x},{0.5_x,1.0_x}};
        ARIADNE_TEST_ASSERT(!optimiser->feasible(D,g,C));
    }

    Void test_nonlinear_equality_feasibility() {
        // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        EffectiveVectorMultivariateFunction h = { 2*x[0]-x[1]+x[0]*x[1]/8 };
        ARIADNE_TEST_PRINT(h);
        ExactBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        ExactBoxType C = ExactBoxType{{0.0_x,0.0_x}};

        ARIADNE_TEST_ASSERT(optimiser->feasible(D,h,C));
    }

    Void test_feasibility_check() {
        EffectiveVectorMultivariateFunction x=EffectiveVectorMultivariateFunction::identity(2);
        ARIADNE_TEST_CONSTRUCT( EffectiveVectorMultivariateFunction, g, ({sqr(x[0])+2*sqr(x[1])-1}) );
        ARIADNE_TEST_CONSTRUCT( ExactBoxType, D, ({{-1.0_x, 1.0_x},{-1.0_x,1.0_x}}) );
        ARIADNE_TEST_CONSTRUCT( ExactBoxType, C, ({{0.0_x,0.0_x}}) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X1, ({{0.296875_x,0.406250_x},{0.593750_x,0.703125_x}},pr) );
        ARIADNE_TEST_ASSERT( definitely(optimiser->contains_feasible_point(D,g,C,X1)) );

        // The following test fails since it is difficult to find the feasible
        // point in the box.
        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X2, ({{0.296875_x,0.406250_x},{0.656250_x,0.656250_x}},pr) );
        ARIADNE_TEST_ASSERT( optimiser->contains_feasible_point(D,g,C,X2) );

        ARIADNE_TEST_CONSTRUCT( FloatDPBoundsVector, X3, ({{0.296875_x,0.406250_x},{0.656250_x,0.687500_x}},pr) );
        ARIADNE_TEST_ASSERT( definitely(optimiser->contains_feasible_point(D,g,C,X3)) );

        ARIADNE_TEST_CONSTRUCT(FloatDPVector, x2, ({0.343750_x,0.656250_x},pr) );
        ARIADNE_TEST_ASSERT( optimiser->validate_feasibility(D,g,C,x2) );
    }

};

class NonlinearAdapterProbe : public NonlinearProgrammingSolverBase {
  public:
    using NonlinearProgrammingSolverBase::minimise;
    using NonlinearProgrammingSolverBase::feasible;

    mutable SizeType minimise_calls=0u;
    mutable SizeType feasible_calls=0u;
    mutable Bool correct_constraints=false;

    NonlinearProgrammingSolverInterface* clone() const override {
        return new NonlinearAdapterProbe(*this);
    }

    Void record(ValidatedVectorMultivariateFunction g, ExactBoxType C) const {
        auto values=g(FloatDPVector({0.5_x},dp));
        correct_constraints = C.size()==2u && values.size()==2u
            && C[0].lower_bound()==-inf && C[0].upper_bound()==0
            && C[1].lower_bound()==0 && C[1].upper_bound()==0
            && definitely(values[0]<FloatDPBounds(0,dp))
            && definitely(values[1]>FloatDPBounds(1,dp));
    }

    Vector<ValidatedNumericType> minimise(
        ValidatedScalarMultivariateFunction, ExactBoxType,
        ValidatedVectorMultivariateFunction g, ExactBoxType C) const override {
        ++minimise_calls;
        record(g,C);
        Vector<ValidatedNumericType> result(1u,dp);
        result[0]=ValidatedNumericType(0.5_x,dp);
        return result;
    }

    ValidatedKleenean feasible(
        ExactBoxType, ValidatedVectorMultivariateFunction g, ExactBoxType C) const override {
        ++feasible_calls;
        record(g,C);
        return true;
    }
};

class TestNonlinearAdapterProbe {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_standard_form());
    }

    Void test_standard_form() {
        NonlinearAdapterProbe solver;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(x[0]);
        ValidatedVectorMultivariateFunction g({x[0]-1});
        ValidatedVectorMultivariateFunction h({x[0]+1});
        ExactBoxType D({{0.0_x,1.0_x}});

        auto minimiser=solver.minimise(f,D,g,h);
        ARIADNE_TEST_EQUAL(minimiser.size(),1u);
        ARIADNE_TEST_EQUAL(solver.minimise_calls,1u);
        ARIADNE_TEST_ASSERT(solver.correct_constraints);
        ARIADNE_TEST_ASSERT(definitely(solver.feasible(D,g,h)));
        ARIADNE_TEST_EQUAL(solver.feasible_calls,1u);
        ARIADNE_TEST_ASSERT(solver.correct_constraints);
    }
};

class TestNonlinearInteriorPointSolverBasics {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_disjoint_constraint_minimisation());
        ARIADNE_TEST_CALL(test_relaxation_parameter());
        ARIADNE_TEST_CALL(test_initial_feasibility_state());
        ARIADNE_TEST_CALL(test_feasible_interior_minimisation());
        ARIADNE_TEST_CALL(test_nonlinear_equality_minimisation_backtracking());
        ARIADNE_TEST_CALL(test_stationary_minimisation_step());
    }

    Void test_feasible_interior_minimisation() {
        // The interior stationary point x=1/2 solves this constrained
        // minimisation problem. Exercise the complete public driver,
        // including its twelve penalty-reduction iterations.
        NonlinearInteriorPointSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction objective=sqr(coordinate[0]-0.5_x);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.25_x,0.75_x}});

        auto optimum=solver.minimise(objective,domain,g,codomain);
        ARIADNE_TEST_EQUAL(optimum.size(),1u);
        ARIADNE_TEST_ASSERT(definitely(optimum[0]==FloatDPBounds(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(contains(domain,optimum)));
        ARIADNE_TEST_ASSERT(definitely(contains(codomain,g(optimum))));
    }

    Void test_nonlinear_equality_minimisation_backtracking() {
        // At x=1/8, linearising h(x)=x^2-1/4 predicts x=17/16,
        // outside [0,1]. Backtracking must shorten the step, even
        // though the nonlinear equality has a feasible root at 1/2.
        NonlinearInteriorPointSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction objective=sqr(coordinate[0]-0.5_x);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ValidatedVectorMultivariateFunction h({sqr(coordinate[0])-0.25_x});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.0_x,1.0_x}});
        FloatDPApproximationVector x({0.125_x},dp);
        FloatDPApproximationVector w({0.125_x},dp);
        FloatDPApproximationVector kappa({0.125_x},dp);
        FloatDPApproximationVector lambda({0.125_x},dp);
        FloatDPApproximation mu(0.25_x,dp);

        solver.minimisation_step(objective,domain,g,codomain,h,
                                  x,w,kappa,lambda,mu);
        ARIADNE_TEST_ASSERT(contains(domain,cast_exact(x)));
        ARIADNE_TEST_ASSERT(contains(codomain,cast_exact(w)));
        ARIADNE_TEST_ASSERT(decide(x[0]>0.5_approx));
        ARIADNE_TEST_ASSERT(decide(x[0]<1.0_approx));
        ARIADNE_TEST_EQUAL(kappa.size(),1u);
        ARIADNE_TEST_EQUAL(lambda.size(),1u);
    }

    Void test_stationary_minimisation_step() {
        NonlinearInteriorPointSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction objective=sqr(coordinates[0]-0.5_x);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ValidatedVectorMultivariateFunction h({coordinates[0]-0.5_x});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType C({{0.25_x,0.75_x}});
        FloatDPApproximationVector x({0.5_x},dp);
        FloatDPApproximationVector w({0.5_x},dp);
        FloatDPApproximationVector kappa({0.0_x},dp);
        FloatDPApproximationVector lambda({0.0_x},dp);
        FloatDPApproximation mu(0.25_x,dp);

        solver.minimisation_step(objective,D,g,C,h,x,w,kappa,lambda,mu);

        ARIADNE_TEST_ASSERT(decide(x[0]==0.5_approx));
        ARIADNE_TEST_ASSERT(decide(w[0]==0.5_approx));
        ARIADNE_TEST_ASSERT(decide(kappa[0]==0.0_approx));
        ARIADNE_TEST_ASSERT(decide(lambda[0]==0.0_approx));
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(x)));
        ARIADNE_TEST_ASSERT(contains(C,cast_exact(w)));
    }

    Void test_initial_feasibility_state() {
        NonlinearInteriorPointSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.25_x,0.75_x}});
        FloatDPApproximationVector primal,dual;

        solver.setup_feasibility(domain,g,codomain,primal,dual);
        ARIADNE_TEST_EQUAL(primal.size(),4u);
        ARIADNE_TEST_EQUAL(dual.size(),1u);
        ARIADNE_TEST_ASSERT(decide(dual[0]==0.5_approx));
        for(SizeType i=0u; i!=primal.size(); ++i) {
            ARIADNE_TEST_ASSERT(decide(primal[i]==0.25_approx));
        }

        ValidatedVectorMultivariateFunction no_constraints(0u,1u);
        solver.setup_feasibility(domain,no_constraints,ExactBoxType{},primal,dual);
        ARIADNE_TEST_EQUAL(primal.size(),2u);
        ARIADNE_TEST_EQUAL(dual.size(),1u);
        ARIADNE_TEST_ASSERT(decide(dual[0]==0.5_approx));
        for(SizeType i=0u; i!=primal.size(); ++i) {
            ARIADNE_TEST_ASSERT(decide(primal[i]==0.5_approx));
        }
    }

    Void test_disjoint_constraint_minimisation() {
        NonlinearInteriorPointSolver solver;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction objective=sqr(x[0]);
        ValidatedVectorMultivariateFunction constraints({x[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{2.0_x,3.0_x}});
        ARIADNE_TEST_THROWS(
            solver.minimise(objective,domain,constraints,codomain),
            InfeasibleProblemException);
    }

    Void test_relaxation_parameter() {
        NonlinearInteriorPointSolver solver;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction constraints({x[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        FloatDPApproximationVector point({0.5_x},dp);
        FloatDPApproximationVector positive_multiplier({2.0_x},dp);
        FloatDPApproximationVector negative_multiplier({-2.0_x},dp);
        ExactBoxType equality({{0.5_x,0.5_x}});
        ExactBoxType upper_bounded({{-inf,0.75_x}});
        ExactBoxType lower_bounded({{0.25_x,inf}});
        ExactBoxType bounded({{0.25_x,0.75_x}});

        ARIADNE_TEST_ASSERT(decide(solver.compute_mu(
            domain,constraints,equality,point,positive_multiplier)==0.0_approx));
        ARIADNE_TEST_ASSERT(is_inf(upper_bounded[0].lower_bound()));
        ARIADNE_TEST_ASSERT(is_inf(lower_bounded[0].upper_bound()));
        ARIADNE_TEST_ASSERT(decide(solver.compute_mu(
            domain,constraints,upper_bounded,point,positive_multiplier)==-0.5_approx));
        ARIADNE_TEST_ASSERT(decide(solver.compute_mu(
            domain,constraints,upper_bounded,point,negative_multiplier)==0.5_approx));
        ARIADNE_TEST_ASSERT(decide(solver.compute_mu(
            domain,constraints,lower_bounded,point,positive_multiplier)==0.5_approx));
        ARIADNE_TEST_ASSERT(decide(solver.compute_mu(
            domain,constraints,bounded,point,negative_multiplier)==0.5_approx));
        ARIADNE_TEST_ASSERT(decide(solver.compute_mu(
            domain,constraints,bounded,point,positive_multiplier)==0.5_approx));
    }
};

class TestNonlinearInfeasibleInteriorPointSolver
{
  private:
    NonlinearInfeasibleInteriorPointSolver optimiser;
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_common_point_feasibility_and_tolerance());
        ARIADNE_TEST_CALL(test_base_validation());
        ARIADNE_TEST_CALL(test_feasible_box_without_equalities());
        ARIADNE_TEST_CALL(test_feasible_box_rejections());
        ARIADNE_TEST_CALL(test_feasible_box_singular_equality());
        ARIADNE_TEST_CALL(test_base_infeasibility());
        ARIADNE_TEST_CALL(test_unconstrained_infeasibility_certificate());
        ARIADNE_TEST_CALL(test_equality_solution_rechecks_inequalities());
        ARIADNE_TEST_CALL(test_infeasible_minimisation_rejects_false_optima());
        ARIADNE_TEST_CALL(test_candidate_linear());
        ARIADNE_TEST_CALL(test_candidate_nonlinear());
        ARIADNE_TEST_CALL(test_candidate_equality());
        ARIADNE_TEST_CALL(test_candidate_infeasible());
        ARIADNE_TEST_CALL(test_candidate_seven_dimensional());
        ARIADNE_TEST_CALL(test_candidate_narrow_codomain());
        ARIADNE_TEST_CALL(test_candidate_anisotropic_domain());
        ARIADNE_TEST_CALL(test_candidate_multiple_constraints());
        ARIADNE_TEST_CALL(test_candidate_nonlinear_seven_dimensional());
        ARIADNE_TEST_CALL(test_candidate_boundary());
        ARIADNE_TEST_CALL(test_candidate_rank_deficient());
        ARIADNE_TEST_CALL(test_candidate_overdetermined());
        ARIADNE_TEST_CALL(test_candidate_deterministic());
    }

    Void test_equality_solution_rechecks_inequalities() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]-0.5_x,x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType permissive({{0.0_x,0.0_x},{0.0_x,0.75_x}});
        ExactBoxType restrictive({{0.0_x,0.0_x},{0.0_x,0.375_x}});
        FloatDPVector initial({0.25_x},dp);

        // The equality has its only root at x=0.5. Both inequality
        // codomains accept the initial guess, but only one accepts
        // the validated equality solution.
        ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(D,g,permissive,initial));
        ARIADNE_TEST_ASSERT(!optimiser.validate_feasibility(D,g,restrictive,initial));
    }

    Void test_unconstrained_infeasibility_certificate() {
        ValidatedVectorMultivariateFunction g(0u,1u);
        ExactBoxType C=ExactBoxType{};
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType outside({{2.0_x,3.0_x}});
        ExactBoxType empty_domain=intersection(D,outside);
        FloatDPVector multipliers(0u,dp);

        ARIADNE_TEST_ASSERT(!D.is_empty());
        ARIADNE_TEST_ASSERT(empty_domain.is_empty());
        ARIADNE_TEST_ASSERT(!optimiser.is_infeasibility_certificate(D,g,C,multipliers));
        ARIADNE_TEST_ASSERT(optimiser.is_infeasibility_certificate(
            empty_domain,g,C,multipliers));
    }

    Void test_base_infeasibility() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType far({{2.0_x,3.0_x}});
        ExactBoxType overlap({{0.0_x,1.0_x}});
        FloatDPVector centre({0.5_x},dp);
        FloatDPVector multiplier({1.0_x},dp);
        ARIADNE_TEST_ASSERT(optimiser.validate_infeasibility(D,g,far,centre,multiplier));
        ARIADNE_TEST_ASSERT(!optimiser.validate_infeasibility(D,g,overlap,centre,multiplier));
    }

    Void test_feasible_box_rejections() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType wide({{-1.0_x,2.0_x}});
        ExactBoxType disjoint({{2.0_x,3.0_x}});
        ExactBoxType overlap({{0.5_x,2.0_x}});
        FloatDPBoundsVector inside({{0.25_x,0.75_x}},dp);
        FloatDPBoundsVector outside({{2.0_x,3.0_x}},dp);
        FloatDPBoundsVector boundary({{0.75_x,1.25_x}},dp);

        auto domain_disjoint=optimiser.contains_feasible_point(D,g,wide,outside);
        ARIADNE_TEST_ASSERT(definitely(not domain_disjoint));

        auto domain_partial=optimiser.contains_feasible_point(D,g,wide,boundary);
        ARIADNE_TEST_ASSERT(!definitely(domain_partial));
        ARIADNE_TEST_ASSERT(!definitely(not domain_partial));

        // X overlaps D, but every point satisfying C lies outside D.
        // A partial overlap cannot establish a feasible point in X.
        ExactBoxType beyond_domain({{1.125_x,1.375_x}});
        auto separated_inside_domain=optimiser.contains_feasible_point(D,g,beyond_domain,boundary);
        ARIADNE_TEST_ASSERT(!definitely(separated_inside_domain));

        auto constraint_disjoint=optimiser.contains_feasible_point(D,g,disjoint,inside);
        ARIADNE_TEST_ASSERT(definitely(not constraint_disjoint));

        auto constraint_overlap=optimiser.contains_feasible_point(D,g,overlap,inside);
        ARIADNE_TEST_ASSERT(!definitely(constraint_overlap));
        ARIADNE_TEST_ASSERT(!definitely(not constraint_overlap));
    }

    Void test_feasible_box_without_equalities() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType C({{-1.0_x,2.0_x}});
        FloatDPBoundsVector box({{0.25_x,0.75_x}},dp);
        ARIADNE_TEST_ASSERT(definitely(optimiser.contains_feasible_point(D,g,C,box)));
    }

    Void test_feasible_box_singular_equality() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]-x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType C({{0.0_x,0.0_x}});
        FloatDPBoundsVector box({{0.25_x,0.75_x}},dp);
        auto status=optimiser.contains_feasible_point(D,g,C,box);
        ARIADNE_TEST_ASSERT(!definitely(status));
        ARIADNE_TEST_ASSERT(!definitely(not status));
    }

    Void test_base_validation() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType C({{0.25_x,0.75_x}});
        FloatDPVector point({0.5_x},dp);
        FloatDPVector outside({0.875_x},dp);
        FloatDPVector multiplier({1.0_x},dp);
        ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(D,g,C,point,multiplier));
        ARIADNE_TEST_ASSERT(!optimiser.validate_feasibility(D,g,C,outside,multiplier));
    }

    Void test_common_point_feasibility_and_tolerance() {
        // A point can be outside C but still satisfy the caller's
        // explicitly requested tolerance. Domain membership remains exact.
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({x[0]});
        ExactBoxType D({{0.0_x,1.0_x}});
        ExactBoxType C({{0.25_x,0.75_x}});
        FloatDPVector inside({0.5_x},dp);
        FloatDPVector outside({0.875_x},dp);
        FloatDPVector beyond_domain({1.5_x},dp);
        FloatDPApproximationVector approx_inside({0.5_x},dp);
        FloatDPApproximationVector approx_outside({0.875_x},dp);
        FloatDPApproximationVector approx_beyond_domain({1.5_x},dp);
        FloatDPApproximation zero_tolerance(0,dp);
        FloatDPApproximation quarter_tolerance(0.25_x,dp);

        ARIADNE_TEST_ASSERT(optimiser.is_feasible_point(D,g,C,inside));
        ARIADNE_TEST_ASSERT(!optimiser.is_feasible_point(D,g,C,outside));
        ARIADNE_TEST_ASSERT(!optimiser.is_feasible_point(D,g,C,beyond_domain));
        ARIADNE_TEST_ASSERT(optimiser.almost_feasible_point(
            D,g,C,approx_inside,zero_tolerance));
        ARIADNE_TEST_ASSERT(!optimiser.almost_feasible_point(
            D,g,C,approx_outside,zero_tolerance));
        ARIADNE_TEST_ASSERT(optimiser.almost_feasible_point(
            D,g,C,approx_outside,quarter_tolerance));
        ARIADNE_TEST_ASSERT(!optimiser.almost_feasible_point(
            D,g,C,approx_beyond_domain,quarter_tolerance));
    }

    Void test_infeasible_minimisation_rejects_false_optima() {
        // g(x)=x takes values in [-1,1]. Neither codomain intersects
        // g(D), so minimisation must not return a certified optimum.
        // Both an infeasibility certificate and an indeterminate
        // outcome are sound results for an iterative solver.
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction objective=sqr(coordinates[0]);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType positive({{1.25_x,1.75_x}});
        ExactBoxType negative({{-1.75_x,-1.25_x}});

        ARIADNE_TEST_ASSERT(definitely(disjoint(Ariadne::apply(g,domain),positive)));
        ARIADNE_TEST_ASSERT(definitely(disjoint(Ariadne::apply(g,domain),negative)));

        for(SizeType i=0u; i!=2u; ++i) {
            const ExactBoxType& codomain=(i==0u ? positive : negative);
            bool rejected=false;
            try {
                (void)optimiser.minimise(objective,domain,g,codomain);
            } catch(const InfeasibleProblemException&) {
                rejected=true;
            } catch(const IndeterminateFeasibilityException&) {
                rejected=true;
            }
            ARIADNE_TEST_ASSERT(rejected);
        }
    }

    Void test_candidate_linear() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({x[0]+x[1]});
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{0.19_pr,0.21_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(definitely(candidate_result.first));
        ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
            D,g,C,cast_exact(candidate_result.second)));
    }

    Void test_candidate_nonlinear() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({sqr(x[0])+x[1]});
        ExactBoxType D({{0.0_x,2.0_x},{0.0_x,2.0_x}});
        ExactBoxType C({{0.49_pr,0.51_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_equality() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({2*x[0]-x[1]+x[0]*x[1]/8});
        ExactBoxType D({{0.0_x,2.0_x},{0.0_x,2.0_x}});
        ExactBoxType C({{0.0_x,0.0_x}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(definitely(candidate_result.first));
        ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
            D,g,C,cast_exact(candidate_result.second)));
    }

    Void test_candidate_infeasible() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({x[0]+x[1]});
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{3.0_x,4.0_x}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(definitely(not candidate_result.first));
    }

    Void test_candidate_seven_dimensional() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(7);
        ValidatedScalarMultivariateFunction sum=x[0];
        for(SizeType i=1u; i!=7u; ++i) {
            sum=sum+x[i];
        }
        ValidatedVectorMultivariateFunction g({sum});
        ExactBoxType D({
            {0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x},
            {0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x},
            {0.0_x,1.0_x}
        });
        ExactBoxType C({{1.175_pr,1.425_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }
    Void test_candidate_narrow_codomain() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({x[0]+x[1]});
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{0.199999_pr,0.200001_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_anisotropic_domain() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({1000*x[0]+x[1]});
        ExactBoxType D({{0.0_x,0.001_pr},{0.0_x,1000.0_x}});
        ExactBoxType C({{0.49_pr,0.51_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        ARIADNE_TEST_ASSERT(not definitely(not candidate_result.first));
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_multiple_constraints() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(3);
        ValidatedVectorMultivariateFunction g({
            x[0]+x[1]+x[2],
            x[0]-x[1]
        });
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{1.19_pr,1.21_pr},{-0.01_pr,0.01_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_nonlinear_seven_dimensional() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(7);
        ValidatedScalarMultivariateFunction expression=
            sqr(x[0])+x[1]+x[2]+x[3]+x[4]+x[5]+x[6];
        ValidatedVectorMultivariateFunction g({expression});
        ExactBoxType D({
            {0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x},
            {0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x},
            {0.0_x,1.0_x}
        });
        ExactBoxType C({{1.29_pr,1.31_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_boundary() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({x[0]+x[1]});
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{0.0_x,0.01_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_rank_deficient() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({
            x[0]+x[1],
            2*x[0]+2*x[1]
        });
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{0.49_pr,0.51_pr},{0.98_pr,1.02_pr}});

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_overdetermined() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({
            x[0]+x[1],
            x[0]-x[1],
            x[0]+2*x[1]
        });
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({
            {0.99_pr,1.01_pr},
            {-0.01_pr,0.01_pr},
            {1.49_pr,1.51_pr}
        });

        auto candidate_result=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(contains(D,cast_exact(candidate_result.second)));
        ARIADNE_TEST_ASSERT(candidate_result.second.size()==D.dimension());
        if(definitely(candidate_result.first)) {
            ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
                D,g,C,cast_exact(candidate_result.second)));
        }
    }

    Void test_candidate_deterministic() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(3);
        ValidatedVectorMultivariateFunction g({
            sqr(x[0])+x[1]+x[2],
            x[0]-x[1]
        });
        ExactBoxType D({{0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType C({{0.99_pr,1.01_pr},{-0.01_pr,0.01_pr}});

        auto first=optimiser.feasible_candidate(D,g,C);
        auto second=optimiser.feasible_candidate(D,g,C);
        ARIADNE_TEST_ASSERT(first.first==second.first);
        ARIADNE_TEST_EQUAL(first.second.size(),second.second.size());
        for(SizeType i=0u; i!=first.second.size(); ++i) {
            ARIADNE_TEST_EQUAL(first.second[i].raw(),second.second[i].raw());
        }
    }

};

class TestNonlinearPenaltyFunctionSolver
{
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_point_feasibility_with_inequalities());
        ARIADNE_TEST_CALL(test_point_feasibility_with_equalities());
        ARIADNE_TEST_CALL(test_stationary_feasibility());
        ARIADNE_TEST_CALL(test_approximate_penalty_step());
        ARIADNE_TEST_CALL(test_one_sided_codomain_feasibility());
        ARIADNE_TEST_CALL(test_approximate_penalty_backtracking());
    }

    Void test_point_feasibility_with_inequalities() {
        NonlinearPenaltyFunctionSolver solver;
        std::unique_ptr<NonlinearPenaltyFunctionSolver> cloned(solver.clone());
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-0.5_x,0.5_x}});
        FloatDPVector inside({0.0_x},dp);
        FloatDPVector outside({0.75_x},dp);
        FloatDPVector multipliers({0.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(solver.check_feasibility(
            domain,g,codomain,inside,multipliers)));
        ARIADNE_TEST_ASSERT(definitely(not cloned->check_feasibility(
            domain,g,codomain,outside,multipliers)));
    }

    Void test_point_feasibility_with_equalities() {
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType interior_domain({{-1.0_x,1.0_x}});
        ExactBoxType boundary_domain({{0.0_x,1.0_x}});
        ExactBoxType equality({{0.0_x,0.0_x}});
        FloatDPVector root({0.0_x},dp);
        FloatDPVector multipliers({0.0_x},dp);

        // The interval Newton inclusion certifies the interior root.
        ARIADNE_TEST_ASSERT(definitely(solver.check_feasibility(
            interior_domain,g,equality,root,multipliers)));

        // A boundary root is feasible, but the symmetric inclusion box
        // crosses the domain boundary: it must not be declared infeasible.
        auto boundary_status=solver.check_feasibility(
            boundary_domain,g,equality,root,multipliers);
        ARIADNE_TEST_ASSERT(not definitely(not boundary_status));
    }

    Void test_stationary_feasibility() {
        // The midpoint solves g(x)=0 strictly inside both boxes.
        // The iterative penalty feasibility routine should preserve
        // this exact stationary solution.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});

        ARIADNE_TEST_ASSERT(definitely(solver.feasible(domain,g,codomain)));
    }

    Void test_one_sided_codomain_feasibility() {
        // The initial slack point is chosen one unit inside the only
        // finite codomain bound, avoiding an infinite midpoint.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType upper_only({{-inf,1.0_x}});
        ExactBoxType lower_only({{-1.0_x,+inf}});

        ARIADNE_TEST_ASSERT(definitely(solver.feasible(domain,g,upper_only)));
        ARIADNE_TEST_ASSERT(definitely(solver.feasible(domain,g,lower_only)));
    }

    Void test_approximate_penalty_backtracking() {
        // For g(x)=4x, x=1/4 and w=3/4, the full Newton update
        // leaves the boxes [-1,1]. Backtracking must reduce the
        // step to keep both updated components admissible.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({4*coordinates[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        FloatDPApproximationVector point({0.25_x},dp);
        FloatDPApproximationVector slack({0.75_x},dp);
        FloatDPApproximation weight(1.0_x,dp);

        solver.feasibility_step(domain,g,codomain,point,slack,weight);
        ARIADNE_TEST_ASSERT(decide(weight==FloatDPApproximation(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(contains(domain,cast_exact(point)));
        ARIADNE_TEST_ASSERT(contains(codomain,cast_exact(slack)));
        ARIADNE_TEST_ASSERT(decide(point[0]<FloatDPApproximation(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(decide(slack[0]<FloatDPApproximation(0.0_x,dp)));
    }

    Void test_approximate_penalty_step() {
        // A regular, one-dimensional problem has an invertible
        // penalised Hessian. Its first step brings w to g(x) while
        // preserving x and halves the penalty weight.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        FloatDPApproximationVector point({0.25_x},dp);
        FloatDPApproximationVector slack({0.0_x},dp);
        FloatDPApproximation weight(1.0_x,dp);

        solver.feasibility_step(domain,g,codomain,point,slack,weight);
        ARIADNE_TEST_ASSERT(decide(weight==FloatDPApproximation(0.5_x,dp)));
        ARIADNE_TEST_ASSERT(contains(domain,cast_exact(point)));
        ARIADNE_TEST_ASSERT(contains(codomain,cast_exact(slack)));
        ARIADNE_TEST_ASSERT(decide(point[0]==FloatDPApproximation(0.25_x,dp)));
        ARIADNE_TEST_ASSERT(decide(slack[0]==FloatDPApproximation(0.25_x,dp)));
    }
};

Int main(Int argc, const char* argv[]) {
    if (not kernel_command_line_interface().acquire(argc,argv)) return -1;

    TestNonlinearAdapterProbe().test();

    NonlinearInfeasibleInteriorPointSolver nlio;
    TestNonlinearProgrammingSolver(nlio).test();
    TestNonlinearInfeasibleInteriorPointSolver().test();
    TestNonlinearInteriorPointSolverBasics().test();
    TestNonlinearPenaltyFunctionSolver().test();
    return ARIADNE_TEST_FAILURES;
    NonlinearInteriorPointSolver nlo;
    TestNonlinearProgrammingSolver(nlo).test();

    ApproximateNonlinearInteriorPointSolver appo;
    TestNonlinearProgrammingSolver(appo).test_nonlinear_equality_feasibility();

    IntervalNonlinearInteriorPointSolver ivlo;
    TestNonlinearProgrammingSolver(ivlo).test_nonlinear_equality_feasibility();
    return ARIADNE_TEST_FAILURES;
}

