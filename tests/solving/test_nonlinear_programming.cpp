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
#include "solving/nonlinear_programming.tpl.hpp"
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
        ARIADNE_TEST_CALL(test_minimisation_dimension_preconditions());
        ARIADNE_TEST_CALL(test_nonlinear_equality_minimisation_backtracking());
        ARIADNE_TEST_CALL(test_stationary_minimisation_step());
#ifndef NDEBUG
        ARIADNE_TEST_CALL(test_minimisation_step_debug_preconditions());
#endif
        ARIADNE_TEST_CALL(test_feasibility_dimension_preconditions());
    }

#ifndef NDEBUG
    Void test_minimisation_step_debug_preconditions() {
        // Verify each existing public debug precondition independently.
        // Release builds omit these checks, so this test is Debug-only.
        NonlinearInteriorPointSolver solver;
        auto c1=ValidatedScalarMultivariateFunction::coordinates(1);
        auto c2=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedScalarMultivariateFunction objective=sqr(c1[0]-0.5_x);
        ValidatedScalarMultivariateFunction wrong_objective=sqr(c2[0]-0.5_x);
        ValidatedVectorMultivariateFunction g({c1[0]});
        ValidatedVectorMultivariateFunction h({c1[0]-0.5_x});
        ValidatedVectorMultivariateFunction wrong_arguments({c2[0]});
        ValidatedVectorMultivariateFunction wrong_results({c1[0],c1[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.25_x,0.75_x}});
        FloatDPApproximationVector point({0.5_x},dp);
        FloatDPApproximationVector slack({0.5_x},dp);
        FloatDPApproximationVector kappa({0.0_x},dp);
        FloatDPApproximationVector lambda({0.0_x},dp);
        FloatDPApproximationVector wrong_size({0.0_x,0.0_x},dp);
        FloatDPApproximationVector outside({1.25_x},dp);
        FloatDPApproximation zero_mu(0.0_x,dp);
        FloatDPApproximation positive_mu(0.25_x,dp);

        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,g,codomain,h,point,wrong_size,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            wrong_objective,domain,g,codomain,h,point,slack,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,wrong_arguments,codomain,h,point,slack,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,g,codomain,wrong_arguments,point,slack,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,wrong_results,codomain,h,point,slack,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,g,codomain,wrong_results,point,slack,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,g,codomain,h,outside,slack,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,g,codomain,h,point,outside,kappa,lambda,positive_mu));
        ARIADNE_TEST_FAIL(solver.minimisation_step(
            objective,domain,g,codomain,h,point,slack,kappa,lambda,zero_mu));
    }
#endif

    Void test_feasibility_dimension_preconditions() {
        // The public feasibility API must reject mismatched dimensions
        // before numerical initialisation or an iterative solver step.
        NonlinearInteriorPointSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.0_x,1.0_x}});
        ExactBoxType wrong_domain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType wrong_codomain({{0.0_x,1.0_x},{0.0_x,1.0_x}});

        ARIADNE_TEST_FAIL(solver.feasible(wrong_domain,g,codomain));
        ARIADNE_TEST_FAIL(solver.feasible(domain,g,wrong_codomain));
    }

    Void test_minimisation_dimension_preconditions() {
        // Independent public input contracts: objective arguments,
        // constraint arguments and the number of constraint bounds.
        // Reject each mismatch before evaluating the nonlinear system.
        NonlinearInteriorPointSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        auto coordinates2=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedScalarMultivariateFunction f=sqr(coordinates[0]-0.5_x);
        ValidatedScalarMultivariateFunction wrong_f=sqr(coordinates2[0]-0.5_x);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ValidatedVectorMultivariateFunction wrong_g({coordinates2[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.25_x,0.75_x}});
        ExactBoxType wrong_codomain({{0.25_x,0.75_x},{0.25_x,0.75_x}});

        ARIADNE_TEST_THROWS(solver.minimise(wrong_f,domain,g,codomain),
                            std::runtime_error);
        ARIADNE_TEST_THROWS(solver.minimise(f,domain,wrong_g,codomain),
                            std::runtime_error);
        ARIADNE_TEST_THROWS(solver.minimise(f,domain,g,wrong_codomain),
                            std::runtime_error);
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
        ARIADNE_TEST_CALL(test_logging_solver_diagnostics());
        ARIADNE_TEST_CALL(test_overdetermined_equality_validation());
        ARIADNE_TEST_CALL(test_logging_diagnostics());
        ARIADNE_TEST_CALL(test_validation_dimension_preconditions());
        ARIADNE_TEST_CALL(test_feasible_box_without_equalities());
        ARIADNE_TEST_CALL(test_feasible_box_rejections());
        ARIADNE_TEST_CALL(test_feasible_box_dimension_preconditions());
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

    Void test_logging_solver_diagnostics() {
        // Validate that enabling diagnostic output leaves observable
        // solver results unchanged for both feasible and infeasible
        // paths. Restore the setting even if a test throws.
        auto& configuration=Logger::instance().configuration();
        const unsigned int old_verbosity=configuration.verbosity();
        struct RestoreVerbosity {
            unsigned int value;
            ~RestoreVerbosity() {
                Logger::instance().configuration().set_verbosity(value);
            }
        } restore{old_verbosity};
        configuration.set_verbosity(10u);

        test_candidate_linear();
        test_candidate_infeasible();
        test_candidate_conflicting_equalities();
        test_candidate_boundary();
        test_candidate_rank_deficient();
        test_candidate_narrow_codomain();
        // These more demanding valid cases may return an indeterminate
        // candidate or reach boundary backtracking. Recheck their
        // established semantic contracts with diagnostic output on.
        test_candidate_anisotropic_domain();
        test_candidate_multiple_constraints();
        // These two existing regression cases can return indeterminate
        // after all candidate iterations. Run them with logging active
        // to exercise diagnostics on that legitimate return path.
        test_candidate_nonlinear();
        test_candidate_seven_dimensional();
        test_candidate_nonlinear_seven_dimensional();
        test_candidate_overdetermined();
        test_infeasible_minimisation_rejects_false_optima();
        test_stationary_minimisation();
    }

    Void test_stationary_minimisation() {
        // The initial midpoint is a stationary constrained minimiser.
        // Check the complete public driver, not just a Newton step;
        // the solver must return a certified feasible candidate.
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction objective=sqr(coordinates[0]-0.5_x);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.25_x,0.75_x}});

        auto optimum=optimiser.minimise(objective,domain,g,codomain);
        ARIADNE_TEST_EQUAL(optimum.size(),1u);
        ARIADNE_TEST_ASSERT(definitely(contains(domain,optimum)));
        ARIADNE_TEST_ASSERT(definitely(contains(codomain,g(optimum))));
        ARIADNE_TEST_ASSERT(definitely(optimum[0]==FloatDPBounds(0.5_x,dp)));

        // Each malformed input violates exactly one of the public
        // dimension preconditions and must be rejected before setup
        // or any numerical iteration.
        ExactBoxType wrong_domain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ValidatedVectorMultivariateFunction wrong_arguments(0u,2u);
        ExactBoxType wrong_codomain({{0.25_x,0.75_x},{0.25_x,0.75_x}});

        ARIADNE_TEST_FAIL(optimiser.minimise(
            objective,wrong_domain,g,codomain));
        ARIADNE_TEST_FAIL(optimiser.minimise(
            objective,domain,wrong_arguments,codomain));
        ARIADNE_TEST_FAIL(optimiser.minimise(
            objective,domain,g,wrong_codomain));
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

    Void test_feasible_box_dimension_preconditions() {
        // These are public API contracts: an invalid dimension must be
        // rejected before indexing g(X) or constructing the Newton step.
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType codomain({{0.0_x,1.0_x}});
        ExactBoxType wrong_domain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType wrong_codomain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        FloatDPBoundsVector box({{0.25_x,0.75_x}},dp);
        FloatDPBoundsVector wrong_box({{0.25_x,0.75_x},{0.25_x,0.75_x}},dp);

        ARIADNE_TEST_FAIL(optimiser.contains_feasible_point(
            wrong_domain,g,codomain,box));
        ARIADNE_TEST_FAIL(optimiser.contains_feasible_point(
            domain,g,wrong_codomain,box));
        ARIADNE_TEST_FAIL(optimiser.contains_feasible_point(
            domain,g,codomain,wrong_box));
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

    Void test_logging_diagnostics() {
        // Exercise public feasibility and infeasibility diagnostics with
        // logging enabled, checking that the reported status remains sound.
        auto& configuration=Logger::instance().configuration();
        unsigned int const old_verbosity=configuration.verbosity();
        configuration.set_verbosity(10u);

        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType inequality({{0.25_x,0.75_x}});
        ExactBoxType equality({{0.5_x,0.5_x}});
        ExactBoxType infeasible({{2.0_x,3.0_x}});
        FloatDPVector feasible_point({0.5_x},dp);
        FloatDPVector infeasible_point({0.875_x},dp);
        FloatDPVector initial({0.25_x},dp);
        FloatDPVector multiplier({1.0_x},dp);
        FloatDPBoundsVector interior({{0.25_x,0.75_x}},dp);

        ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
            domain,g,inequality,feasible_point));
        ARIADNE_TEST_ASSERT(!optimiser.validate_feasibility(
            domain,g,inequality,infeasible_point));
        ARIADNE_TEST_ASSERT(optimiser.validate_feasibility(
            domain,g,equality,initial));
        ARIADNE_TEST_ASSERT(definitely(optimiser.contains_feasible_point(
            domain,g,inequality,interior)));
        ARIADNE_TEST_ASSERT(definitely(optimiser.contains_feasible_point(
            domain,g,equality,interior)));
        ARIADNE_TEST_ASSERT(optimiser.validate_infeasibility(
            domain,g,infeasible,initial,multiplier));
        ARIADNE_TEST_ASSERT(!optimiser.validate_infeasibility(
            domain,g,inequality,initial,multiplier));
        ARIADNE_TEST_ASSERT(optimiser.is_infeasibility_certificate(
            domain,g,infeasible,multiplier));

        configuration.set_verbosity(old_verbosity);
    }

    Void test_validation_dimension_preconditions() {
        // Each public dimension contract rejects a different malformed
        // caller input before evaluating the objective or its derivatives.
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType constraints({{0.0_x,1.0_x}});
        ExactBoxType wrong_domain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        ExactBoxType wrong_constraints({{0.0_x,1.0_x},{0.0_x,1.0_x}});
        FloatDPVector point({0.5_x},dp);
        FloatDPVector wrong_point({0.5_x,0.5_x},dp);
        FloatDPVector multiplier({1.0_x},dp);
        FloatDPVector wrong_multiplier({1.0_x,1.0_x},dp);

        ARIADNE_TEST_FAIL(optimiser.validate_feasibility(
            wrong_domain,g,constraints,point));
        ARIADNE_TEST_FAIL(optimiser.validate_feasibility(
            domain,g,wrong_constraints,point));
        ARIADNE_TEST_FAIL(optimiser.validate_feasibility(
            domain,g,constraints,wrong_point));

        ARIADNE_TEST_FAIL(optimiser.validate_infeasibility(
            wrong_domain,g,constraints,point,multiplier));
        ARIADNE_TEST_FAIL(optimiser.validate_infeasibility(
            domain,g,wrong_constraints,point,multiplier));
        ARIADNE_TEST_FAIL(optimiser.validate_infeasibility(
            domain,g,constraints,wrong_point,multiplier));
        ARIADNE_TEST_FAIL(optimiser.validate_infeasibility(
            domain,g,constraints,point,wrong_multiplier));
    }

    Void test_overdetermined_equality_validation() {
        // Two distinct equality constraints on one variable cannot
        // both hold. The Newton system has two rows but rank at most one.
        // Failure to certify must return false, not throw on a singular
        // correction matrix or report a spurious feasible point.
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinates[0],coordinates[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType equalities({{0.0_x,0.0_x},{1.0_x,1.0_x}});
        FloatDPVector initial({0.0_x},dp);

        ARIADNE_TEST_ASSERT(!optimiser.validate_feasibility(
            domain,g,equalities,initial));
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

        // Independent public contracts of feasible_candidate():
        // the function input dimension must match D, and its output
        // dimension must match C. Reuse the already validated problem
        // to avoid introducing additional numerical trajectories.
        ExactBoxType wrong_domain({{0.0_x,1.0_x}});
        ExactBoxType wrong_codomain({{0.19_pr,0.21_pr},{0.19_pr,0.21_pr}});
        ARIADNE_TEST_FAIL(optimiser.feasible_candidate(wrong_domain,g,C));
        ARIADNE_TEST_FAIL(optimiser.feasible_candidate(D,g,wrong_codomain));
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

    Void test_candidate_conflicting_equalities() {
        // Both constraints evaluate to the same coordinate, yet their
        // singleton target intervals are disjoint. No feasible point
        // exists. The iterative method may establish infeasibility or
        // return indeterminate, but must never claim feasibility.
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0],coordinate[0]});
        ExactBoxType domain({{0.0_x,1.0_x}});
        ExactBoxType contradictory({{0.0_x,0.0_x},{1.0_x,1.0_x}});

        auto candidate_result=optimiser.feasible_candidate(domain,g,contradictory);
        ARIADNE_TEST_ASSERT(!definitely(candidate_result.first));
        ARIADNE_TEST_EQUAL(candidate_result.second.size(),domain.dimension());
        ARIADNE_TEST_ASSERT(contains(domain,cast_exact(candidate_result.second)));
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
        ARIADNE_TEST_CALL(test_logging_diagnostics());
        ARIADNE_TEST_CALL(test_point_feasibility_with_inequalities());
        ARIADNE_TEST_CALL(test_point_feasibility_rounding_overlap());
        ARIADNE_TEST_CALL(test_penalty_minimisation_unconstrained());
        ARIADNE_TEST_CALL(test_penalty_minimisation_lower_domain_projection());
        ARIADNE_TEST_CALL(test_penalty_minimisation_active_constraint());
        ARIADNE_TEST_CALL(test_penalty_minimisation_detects_infeasibility());
        ARIADNE_TEST_CALL(test_penalty_minimisation_input_contracts());
        ARIADNE_TEST_CALL(test_penalty_minimisation_upper_codomain());
        ARIADNE_TEST_CALL(test_penalty_minimisation_boundary_requires_certificate());
        ARIADNE_TEST_CALL(test_penalty_minimisation_lower_boundary_requires_certificate());
        ARIADNE_TEST_CALL(test_point_feasibility_input_contracts());
        ARIADNE_TEST_CALL(test_point_feasibility_with_equalities());
        ARIADNE_TEST_CALL(test_point_feasibility_multiplier_infeasibility_certificate());
        ARIADNE_TEST_CALL(test_stationary_feasibility());
        ARIADNE_TEST_CALL(test_approximate_penalty_step());
        ARIADNE_TEST_CALL(test_validated_penalty_feasibility_contraction());
        ARIADNE_TEST_CALL(test_validated_penalty_step_input_contracts());
        ARIADNE_TEST_CALL(test_validated_penalty_negative_derivative());
        ARIADNE_TEST_CALL(test_validated_penalty_derivative_containing_zero());
        ARIADNE_TEST_CALL(test_validated_penalty_multidimensional_contraction());
        ARIADNE_TEST_CALL(test_validated_penalty_inconsistent_enclosures());
        ARIADNE_TEST_CALL(test_one_sided_codomain_feasibility());
        ARIADNE_TEST_CALL(test_unbounded_codomain_feasibility());
        ARIADNE_TEST_CALL(test_primal_dual_penalty_equality_step());
        ARIADNE_TEST_CALL(test_primal_dual_penalty_step_input_contracts());
        ARIADNE_TEST_CALL(test_primal_dual_penalty_backtracking());
        ARIADNE_TEST_CALL(test_approximate_penalty_backtracking());
    }

    Void test_logging_diagnostics() {
        // Recheck already-established semantic contracts with logging
        // enabled. Diagnostic output must not change solver outcomes.
        // Restore the global configuration even if an assertion throws.
        auto& configuration=Logger::instance().configuration();
        const unsigned int old_verbosity=configuration.verbosity();
        struct RestoreVerbosity {
            unsigned int value;
            ~RestoreVerbosity() {
                Logger::instance().configuration().set_verbosity(value);
            }
        } restore{old_verbosity};
        configuration.set_verbosity(10u);

        // Cover the public driver, all three step overloads, validated
        // contraction and the true/false/indeterminate certificate paths.
        test_penalty_minimisation_unconstrained();
        test_stationary_feasibility();
        test_approximate_penalty_step();
        test_primal_dual_penalty_equality_step();
        test_primal_dual_penalty_backtracking();
        test_validated_penalty_feasibility_contraction();
        test_point_feasibility_with_inequalities();
        test_point_feasibility_with_equalities();
        test_point_feasibility_multiplier_infeasibility_certificate();
        test_point_feasibility_rounding_overlap();
    }

    Void test_point_feasibility_input_contracts() {
        // check_feasibility is a public entry point. Each invalid call
        // violates exactly one dimension contract and must be rejected
        // before function evaluation or multiplier arithmetic.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-0.5_x,0.5_x}});
        ExactBoxType two_dimensional_domain({{-1.0_x,1.0_x},
                                             {-1.0_x,1.0_x}});
        ExactBoxType two_dimensional_codomain({{-0.5_x,0.5_x},
                                               {-0.5_x,0.5_x}});
        FloatDPVector point({0.0_x},dp);
        FloatDPVector multipliers({0.0_x},dp);
        FloatDPVector wrong_size({0.0_x,0.0_x},dp);

        ARIADNE_TEST_FAIL(solver.check_feasibility(
            two_dimensional_domain,g,codomain,point,multipliers));
        ARIADNE_TEST_FAIL(solver.check_feasibility(
            domain,g,two_dimensional_codomain,point,multipliers));
        ARIADNE_TEST_FAIL(solver.check_feasibility(
            domain,g,codomain,wrong_size,multipliers));
        ARIADNE_TEST_FAIL(solver.check_feasibility(
            domain,g,codomain,point,wrong_size));
    }

    Void test_penalty_minimisation_lower_domain_projection() {
        // The unconstrained minimiser of (x+2)^2 is x=-2.
        // On [-1,1], its unique constrained minimiser is x=-1.
        // The initial Gauss-Newton step crosses the lower domain bound,
        // so the projected line search must clamp the candidate to -1.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]+Real(2.0_x));
        ValidatedVectorMultivariateFunction g(0u,1u);
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain(0u);

        FloatDPBoundsVector optimum=solver.minimise(f,domain,g,codomain);
        ARIADNE_TEST_EQUAL(optimum.size(),1u);
        ARIADNE_TEST_ASSERT(same(optimum[0],FloatDPBounds(-1.0_x,dp)));
    }

    Void test_penalty_minimisation_unconstrained() {
        // The unique minimiser of (x-1/4)^2 on [-1,1] is x=1/4.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]-Real(0.25_x));
        ValidatedVectorMultivariateFunction g(0u,1u);
        ExactBoxType domain({{-1.0_x,1.0_x}});
        // An empty exact box must be constructed by dimension, not from
        // an empty initializer list (which has no precision information).
        ExactBoxType codomain(0u);
        FloatDPBoundsVector optimum=solver.minimise(f,domain,g,codomain);
        ARIADNE_TEST_EQUAL(optimum.size(),1u);
        auto midpoint_optimum=midpoint(optimum);
        ARIADNE_TEST_ASSERT(decide(abs(midpoint_optimum[0]-
            FloatDPApproximation(0.25_x,dp))<FloatDPApproximation(0.015625_x,dp)));
    }

    Void test_penalty_minimisation_active_constraint() {
        // Start at x=0, which violates g(x)>=1/2; the constrained
        // optimum of (x-3/4)^2 lies strictly inside [1/2,1].
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]-Real(0.75_x));
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{0.5_x,1.0_x}});
        FloatDPBoundsVector optimum=solver.minimise(f,domain,g,codomain);
        ARIADNE_TEST_ASSERT(definitely(contains(codomain,g(optimum))));
        auto midpoint_optimum=midpoint(optimum);
        ARIADNE_TEST_ASSERT(decide(abs(midpoint_optimum[0]-
            FloatDPApproximation(0.75_x,dp))<FloatDPApproximation(0.015625_x,dp)));
    }

    Void test_penalty_minimisation_input_contracts() {
        // These are distinct violations of the documented public
        // minimisation contract. Reject them before evaluating f or g.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinates=ValidatedScalarMultivariateFunction::coordinates(1);
        auto coordinates2=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedScalarMultivariateFunction f=sqr(coordinates[0]);
        ValidatedScalarMultivariateFunction wrong_f=sqr(coordinates2[0]);
        ValidatedVectorMultivariateFunction g({coordinates[0]});
        ValidatedVectorMultivariateFunction wrong_g({coordinates2[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        ExactBoxType wrong_codomain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        ExactBoxType unbounded_below({{-inf,1.0_x}});
        ExactBoxType unbounded_above({{-1.0_x,+inf}});

        ARIADNE_TEST_FAIL(solver.minimise(wrong_f,domain,g,codomain));
        ARIADNE_TEST_FAIL(solver.minimise(f,domain,wrong_g,codomain));
        ARIADNE_TEST_FAIL(solver.minimise(f,domain,g,wrong_codomain));
        ARIADNE_TEST_FAIL(solver.minimise(f,unbounded_below,g,codomain));
        ARIADNE_TEST_FAIL(solver.minimise(f,unbounded_above,g,codomain));
    }

    Void test_penalty_minimisation_upper_codomain() {
        // g(x)=x initially exceeds the upper constraint -1/2 at x=0.
        // The unconstrained optimum -3/4 lies strictly within the
        // feasible interval [-1,-1/2] and must remain optimal.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]+Real(0.75_x));
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,-0.5_x}});
        FloatDPBoundsVector optimum=solver.minimise(f,domain,g,codomain);
        ARIADNE_TEST_ASSERT(definitely(contains(codomain,g(optimum))));
        auto midpoint_optimum=midpoint(optimum);
        ARIADNE_TEST_ASSERT(decide(abs(midpoint_optimum[0]+
            FloatDPApproximation(0.75_x,dp))<FloatDPApproximation(0.015625_x,dp)));
    }

    Void test_penalty_minimisation_boundary_requires_certificate() {
        // The constrained minimum of (x-1)^2 with x<=0 is at x=0.
        // At every finite penalty weight rho, the penalised minimiser
        // lies at x=2/(rho+2)>0. Such an infeasible numerical iterate
        // cannot be silently returned as a validated solution.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]-Real(1.0_x));
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,0.0_x}});

        ARIADNE_TEST_THROWS(
            solver.minimise(f,domain,g,codomain),
            IndeterminateFeasibilityException);
    }

    Void test_penalty_minimisation_lower_boundary_requires_certificate() {
        // The constrained minimum of (x+1)^2 with x>=1/2 is at x=1/2.
        // The exact minimiser of the quadratic-penalty objective for
        // finite rho is 1/2-3/(rho+2), which violates the lower bound.
        // A numerical iterate must not be returned as a validated
        // feasible optimum merely because its violation is small.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]+Real(1.0_x));
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{0.5_x,1.0_x}});

        ARIADNE_TEST_THROWS(
            solver.minimise(f,domain,g,codomain),
            IndeterminateFeasibilityException);
    }

    Void test_penalty_minimisation_detects_infeasibility() {
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedScalarMultivariateFunction f=sqr(coordinate[0]);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType impossible({{2.0_x,3.0_x}});
        ARIADNE_TEST_FAIL(solver.minimise(f,domain,g,impossible));
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

    Void test_point_feasibility_rounding_overlap() {
        // sqrt(2) is not exactly representable in binary floating-point.
        // Its validated point evaluation can straddle the lower bound
        // of C. This neither certifies feasibility of x=2 nor proves
        // the problem infeasible, since g(D) intersects C for x>2.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({sqrt(coordinate[0])});
        ExactBoxType domain({{1.0_x,3.0_x}});
        ExactBoxType codomain({{1.4142135623730951_pr,1.5_x}});
        FloatDPVector point({2.0_x},dp);
        FloatDPVector multiplier({0.0_x},dp);

        const auto status=solver.check_feasibility(
            domain,g,codomain,point,multiplier);
        ARIADNE_TEST_ASSERT(not definitely(status));
        ARIADNE_TEST_ASSERT(not definitely(not status));
    }

    Void test_point_feasibility_multiplier_infeasibility_certificate() {
        // g(D)=[-1,1] is disjoint from C={5/2}, so the feasibility
        // problem has no solution. At the supplied centering point
        // x=5/2 (which need not belong to D), g(x) exactly satisfies
        // the equality. The immediate pointwise check therefore cannot
        // prove infeasibility; the multiplier y=1 separates g(D) and C.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType equality({{2.5_x,2.5_x}});
        FloatDPVector center({2.5_x},dp);
        FloatDPVector certificate({1.0_x},dp);
        FloatDPVector no_certificate({0.0_x},dp);

        ARIADNE_TEST_ASSERT(definitely(not solver.check_feasibility(
            domain,g,equality,center,certificate)));

        // The same infeasible problem has no separating information
        // with y=0. Do not mistake absence of a certificate for
        // evidence of feasibility.
        auto unproved=solver.check_feasibility(
            domain,g,equality,center,no_certificate);
        ARIADNE_TEST_ASSERT(not definitely(unproved));
        ARIADNE_TEST_ASSERT(not definitely(not unproved));
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

    Void test_unbounded_codomain_feasibility() {
        // Both bounds on the slack variable are infinite. Its initial
        // value must be finite (zero), not the midpoint of [-inf,+inf].
        // g(x)=x at the domain midpoint is feasible independently of
        // how far the unconstrained slack variable may range.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType unconstrained_codomain({{-inf,+inf}});

        ARIADNE_TEST_ASSERT(definitely(
            solver.feasible(domain,g,unconstrained_codomain)));
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

    Void test_primal_dual_penalty_step_input_contracts() {
        // This public overload accepts a primal point, dual multipliers
        // and constraint slack. Reject each inconsistent dimension before
        // evaluating a function or assembling a Newton system.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        ExactBoxType wrong_domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        ExactBoxType wrong_codomain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        FloatDPApproximationVector point({0.25_x},dp);
        FloatDPApproximationVector multipliers({0.0_x},dp);
        FloatDPApproximationVector slack({0.0_x},dp);
        FloatDPApproximationVector wrong_size({0.0_x,0.0_x},dp);

        ARIADNE_TEST_FAIL(solver.feasibility_step(
            wrong_domain,g,codomain,point,multipliers,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,wrong_codomain,point,multipliers,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,codomain,wrong_size,multipliers,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,codomain,point,wrong_size,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,codomain,point,multipliers,wrong_size));
    }

    Void test_primal_dual_penalty_equality_step() {
        // The three-vector penalty step must preserve a singleton slack
        // constraint, while reducing the residual g(x)-w for g(x)=x.
        // This exercises the equality row replacement in its Newton system.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType equality({{0.0_x,0.0_x}});
        FloatDPApproximationVector point({0.25_x},dp);
        FloatDPApproximationVector multipliers({0.0_x},dp);
        FloatDPApproximationVector slack({0.0_x},dp);

        solver.feasibility_step(domain,g,equality,point,multipliers,slack);
        ARIADNE_TEST_ASSERT(contains(domain,cast_exact(point)));
        ARIADNE_TEST_ASSERT(contains(equality,cast_exact(slack)));
        ARIADNE_TEST_ASSERT(decide(slack[0]==FloatDPApproximation(0.0_x,dp)));
        ARIADNE_TEST_ASSERT(decide(abs(point[0])<FloatDPApproximation(0.05_a,dp)));
    }

    Void test_primal_dual_penalty_backtracking() {
        // For g(x)=x^2, the full Newton step from x=0.4, w=0.2,
        // y=2 exceeds both [-1,1] bounds. A 0.75-scaled step lies
        // strictly inside each box, so line search must reduce alpha.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({sqr(coordinate[0])});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        FloatDPApproximationVector point({FloatDPApproximation(0.4_a,dp)});
        FloatDPApproximationVector multipliers({2.0_x},dp);
        FloatDPApproximationVector slack({FloatDPApproximation(0.2_a,dp)});

        solver.feasibility_step(domain,g,codomain,point,multipliers,slack);
        ARIADNE_TEST_ASSERT(contains(domain,cast_exact(point)));
        ARIADNE_TEST_ASSERT(contains(codomain,cast_exact(slack)));
        // Expected reduced step is approximately (-0.843,-0.824).
        ARIADNE_TEST_ASSERT(decide(point[0]<FloatDPApproximation(-0.7_a,dp)));
        ARIADNE_TEST_ASSERT(decide(point[0]>FloatDPApproximation(-0.95_a,dp)));
        ARIADNE_TEST_ASSERT(decide(slack[0]<FloatDPApproximation(-0.7_a,dp)));
        ARIADNE_TEST_ASSERT(decide(slack[0]>FloatDPApproximation(-0.95_a,dp)));
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

    Void test_validated_penalty_step_input_contracts() {
        // Every public vector and function dimension must be compatible
        // before evaluating g or refining any interval enclosures.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        ExactBoxType wrong_domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        ExactBoxType wrong_codomain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        FloatDPBoundsVector point(1u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector slack(1u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector wrong_point(2u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector wrong_slack(2u,FloatDPBounds(-1,+1,dp));

        ARIADNE_TEST_FAIL(solver.feasibility_step(
            wrong_domain,g,codomain,point,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,wrong_codomain,point,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,codomain,wrong_point,slack));
        ARIADNE_TEST_FAIL(solver.feasibility_step(
            domain,g,codomain,point,wrong_slack));
    }

    Void test_validated_penalty_negative_derivative() {
        // g(x)=-x is strictly decreasing. Contracting g(x)=1/4
        // must select x=-1/4, exercising the negative-Jacobian case.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({-coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType equality({{0.25_x,0.25_x}});
        FloatDPBoundsVector point(1u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector slack(1u,FloatDPBounds(-1,+1,dp));

        solver.feasibility_step(domain,g,equality,point,slack);
        ARIADNE_TEST_ASSERT(definitely(point[0]==FloatDPBounds(-0.25_x,dp)));
        ARIADNE_TEST_ASSERT(definitely(slack[0]==FloatDPBounds(0.25_x,dp)));
    }

    Void test_validated_penalty_derivative_containing_zero() {
        // For g(x)=x^2 on [-1,1], the derivative enclosure [-2,2]
        // includes zero. Inverse-derivative contraction is therefore
        // invalid; the step must preserve a sound enclosure of x=0.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({sqr(coordinate[0])});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType equality({{0.0_x,0.0_x}});
        FloatDPBoundsVector point(1u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector slack(1u,FloatDPBounds(-1,+1,dp));

        solver.feasibility_step(domain,g,equality,point,slack);
        // Equality of overlapping validated intervals is indeterminate.
        // Compare their representations to verify that the enclosure
        // was preserved exactly when the derivative contains zero.
        ARIADNE_TEST_ASSERT(same(point[0],FloatDPBounds(-1,+1,dp)));
        ARIADNE_TEST_ASSERT(definitely(slack[0]==FloatDPBounds(0,dp)));
    }

    Void test_validated_penalty_multidimensional_contraction() {
        // g(x0,x1)=x0+x1=0 has infinitely many solutions in
        // [-1,1]^2. The validated step must contract the slack to
        // zero without applying a one-dimensional inverse derivative.
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction g({coordinate[0]+coordinate[1]});
        ExactBoxType domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});
        ExactBoxType equality({{0.0_x,0.0_x}});
        FloatDPBoundsVector point(2u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector slack(1u,FloatDPBounds(-1,+1,dp));

        solver.feasibility_step(domain,g,equality,point,slack);
        ARIADNE_TEST_ASSERT(definitely(slack[0]==FloatDPBounds(0.0_x,dp)));
        ARIADNE_TEST_EQUAL(point.size(),2u);
    }

    Void test_validated_penalty_feasibility_contraction() {
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType equality({{0.0_x,0.0_x}});
        FloatDPBoundsVector point(1u,FloatDPBounds(-1,+1,dp));
        FloatDPBoundsVector slack(1u,FloatDPBounds(-1,+1,dp));

        solver.feasibility_step(domain,g,equality,point,slack);
        // The mean-value contractor for g(x)=x narrows both unknowns
        // to the exactly feasible singleton zero.
        ARIADNE_TEST_ASSERT(definitely(point[0]==FloatDPBounds(0,dp)));
        ARIADNE_TEST_ASSERT(definitely(slack[0]==FloatDPBounds(0,dp)));
    }

    Void test_validated_penalty_inconsistent_enclosures() {
        NonlinearPenaltyFunctionSolver solver;
        auto coordinate=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction g({coordinate[0]});
        ExactBoxType domain({{-1.0_x,1.0_x}});
        ExactBoxType codomain({{-1.0_x,1.0_x}});
        FloatDPBoundsVector point(1u,FloatDPBounds(-1.0_x,-0.5_x,dp));
        FloatDPBoundsVector slack(1u,FloatDPBounds(0.5_x,1.0_x,dp));

        // No x in [-1,-0.5] can satisfy g(x)=w with w in [0.5,1].
        ARIADNE_TEST_THROWS(
            solver.feasibility_step(domain,g,codomain,point,slack),
            InfeasibleProblemException);
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


Void test_nonlinear_header_contracts()
{
    const IndeterminateFeasibilityException indeterminate_exception;
    const DegenerateNonlinearFeasibilityProblemException degenerate;
    ARIADNE_TEST_EQUAL(std::string(indeterminate_exception.what()),std::string("IndeterminateFeasibilityException"));
    ARIADNE_TEST_EQUAL(std::string(degenerate.what()),std::string("DegenerateNonlinearFeasibilityProblemException"));

    // Clone through the public solver interface to check dynamic dispatch,
    // concrete type preservation and independent ownership.
    NonlinearInteriorPointSolver interior;
    ApproximateNonlinearInteriorPointSolver approximate;
    IntervalNonlinearInteriorPointSolver interval;
    std::unique_ptr<NonlinearProgrammingSolverInterface> c0(interior.clone());
    std::unique_ptr<NonlinearProgrammingSolverInterface> c1(
        static_cast<const NonlinearProgrammingSolverInterface&>(approximate).clone());
    std::unique_ptr<NonlinearProgrammingSolverInterface> c2(
        static_cast<const NonlinearProgrammingSolverInterface&>(interval).clone());
    ARIADNE_TEST_ASSERT(dynamic_cast<NonlinearInteriorPointSolver*>(c0.get())!=nullptr);
    ARIADNE_TEST_ASSERT(dynamic_cast<ApproximateNonlinearInteriorPointSolver*>(c1.get())!=nullptr);
    ARIADNE_TEST_ASSERT(dynamic_cast<IntervalNonlinearInteriorPointSolver*>(c2.get())!=nullptr);
}

Void test_nonlinear_template_feasibility_products()
{
    using X=FloatDPApproximation;
    Matrix<X> A(1u,1u,X(0.0_x,dp));
    A[0][0]=X(2.0_x,dp);
    Vector<X> v(4u,dp);
    v[0]=X(1.0_x,dp);
    v[1]=X(2.0_x,dp);
    v[2]=X(3.0_x,dp);
    v[3]=X(4.0_x,dp);

    // [A -A I -I; 1 1 1 1] * [1 2 3 4] = [-3,10].
    auto result=Detail::feasibility_mul(A,v);
    ARIADNE_TEST_EQUAL(result.size(),2u);
    ARIADNE_TEST_EQUAL(result[0].raw(),X(-3.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(result[1].raw(),X(10.0_x,dp).raw());

    Vector<X> w(2u,dp);
    w[0]=X(3.0_x,dp);
    w[1]=X(5.0_x,dp);
    // Transpose product = [11,-1,8,2].
    auto transposed=Detail::feasibility_trmul(A,w);
    ARIADNE_TEST_EQUAL(transposed.size(),4u);
    ARIADNE_TEST_EQUAL(transposed[0].raw(),X(11.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(transposed[1].raw(),X(-1.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(transposed[2].raw(),X(8.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(transposed[3].raw(),X(2.0_x,dp).raw());

    Matrix<X> H(1u,1u,X(7.0_x,dp));
    // With D=[1,2,3,4], the symmetric augmented product is
    // [26,-3; -3,10].
    auto augmented=Detail::feasibility_adat(H,A,v);
    ARIADNE_TEST_EQUAL(augmented.row_size(),2u);
    ARIADNE_TEST_EQUAL(augmented.column_size(),2u);
    ARIADNE_TEST_EQUAL(augmented[0][0].raw(),X(26.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(augmented[0][1].raw(),X(-3.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(augmented[1][0].raw(),X(-3.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(augmented[1][1].raw(),X(10.0_x,dp).raw());

    // Elementwise helpers implement the same arithmetic on each component.
    Vector<X> x(2u,dp),d(2u,dp);
    x[0]=X(2.0_x,dp); x[1]=X(4.0_x,dp);
    d[0]=X(4.0_x,dp); d[1]=X(8.0_x,dp);
    auto shifted=Detail::esub(x,X(1.0_x,dp));
    auto multiplied=Detail::emul(x,d);
    auto divided=Detail::ediv(x,d);
    ARIADNE_TEST_EQUAL(shifted[0].raw(),X(1.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(shifted[1].raw(),X(3.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(multiplied[0].raw(),X(8.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(multiplied[1].raw(),X(32.0_x,dp).raw());
    ARIADNE_TEST_EQUAL(divided[0].raw(),X(0.5_x,dp).raw());
    ARIADNE_TEST_EQUAL(divided[1].raw(),X(0.5_x,dp).raw());
}

// Verify the exact interval hull of scalar samples, including a singleton.
Void test_nonlinear_header_interval_hull()
{
    RawFloatDPVector samples(3u,FloatDP(0.0_x,dp));
    samples[0]=FloatDP(4.0_x,dp);
    samples[1]=FloatDP(-3.0_x,dp);
    samples[2]=FloatDP(2.0_x,dp);
    ARIADNE_TEST_ASSERT(definitely(
        Detail::eivl(samples)==ExactIntervalType(-3.0_x,4.0_x)));

    RawFloatDPVector singleton(1u,FloatDP(2.0_x,dp));
    ARIADNE_TEST_ASSERT(definitely(
        Detail::eivl(singleton)==ExactIntervalType(2.0_x,2.0_x)));
}

// Symmetric A^T D A must accumulate, and reject incorrect dimensions.
Void test_nonlinear_header_atda()
{
    using X=FloatDPApproximation;
    const X zero(0.0_x,dp);
    Matrix<X> A(1u,1u,zero);
    A[0][0]=X(2.0_x,dp);
    Vector<X> diagonal(1u,dp);
    diagonal[0]=X(3.0_x,dp);
    Matrix<X> result(1u,1u,zero);
    result[0][0]=X(5.0_x,dp);
    Detail::atda(result,A,diagonal);
    ARIADNE_TEST_EQUAL(result[0][0].raw(),X(17.0_x,dp).raw());

    Matrix<X> nonsquare(1u,2u,zero);
    Matrix<X> wrong_columns(2u,2u,zero);
    Vector<X> wrong_diagonal(2u,dp);
    wrong_diagonal[0]=X(1.0_x,dp);
    wrong_diagonal[1]=X(2.0_x,dp);
    ARIADNE_TEST_FAIL(Detail::atda(nonsquare,A,diagonal));
    ARIADNE_TEST_FAIL(Detail::atda(wrong_columns,A,diagonal));
    ARIADNE_TEST_FAIL(Detail::atda(result,A,wrong_diagonal));
}

// Reject invalid dimensions before constructing augmented feasibility
// vectors and normal matrices; each failure protects a different index.
Void test_nonlinear_header_feasibility_preconditions()
{
    using X=FloatDPApproximation;
    const X zero(0.0_x,dp);
    RawFloatDPVector empty_samples(0u,FloatDP(0.0_x,dp));
    ARIADNE_TEST_FAIL(Detail::eivl(empty_samples));

    Matrix<X> A(1u,1u,zero);
    A[0][0]=X(2.0_x,dp);
    Vector<X> short_primal(3u,dp);
    Vector<X> long_primal(5u,dp);
    // One equation and one variable require four primal coordinates.
    ARIADNE_TEST_FAIL(Detail::feasibility_mul(A,short_primal));
    ARIADNE_TEST_FAIL(Detail::feasibility_mul(A,long_primal));

    Vector<X> short_dual(1u,dp);
    Vector<X> long_dual(3u,dp);
    // The transpose takes one multiplier plus the final scalar.
    ARIADNE_TEST_FAIL(Detail::feasibility_trmul(A,short_dual));
    ARIADNE_TEST_FAIL(Detail::feasibility_trmul(A,long_dual));

    Matrix<X> H(1u,1u,zero);
    Matrix<X> wrong_rows(2u,1u,zero);
    Matrix<X> wrong_columns(1u,2u,zero);
    Vector<X> valid_diagonal(4u,dp);
    Vector<X> short_diagonal(3u,dp);
    // The augmented normal matrix requires a square m-by-m H and
    // exactly 2*(m+n) diagonal elements. Check independently.
    ARIADNE_TEST_FAIL(Detail::feasibility_adat(wrong_rows,A,valid_diagonal));
    ARIADNE_TEST_FAIL(Detail::feasibility_adat(wrong_columns,A,valid_diagonal));
    ARIADNE_TEST_FAIL(Detail::feasibility_adat(H,A,short_diagonal));
}

Int main(Int argc, const char* argv[]) {
    if (not kernel_command_line_interface().acquire(argc,argv)) return -1;

    ARIADNE_TEST_CALL(test_nonlinear_header_contracts());
    ARIADNE_TEST_CALL(test_nonlinear_template_feasibility_products());
    ARIADNE_TEST_CALL(test_nonlinear_header_interval_hull());
    ARIADNE_TEST_CALL(test_nonlinear_header_atda());
    ARIADNE_TEST_CALL(test_nonlinear_header_feasibility_preconditions());
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

