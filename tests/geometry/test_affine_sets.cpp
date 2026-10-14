/***************************************************************************
 *            test_affine_sets.cpp
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

#include <iostream>
#include "function/function.hpp"
#include "geometry/box.hpp"
#include "geometry/affine_set.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_paving.hpp"
#include "drawing/figure.hpp"
#include "drawing/geometry2d.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;
using namespace std;

static const Colour colour(0.5,1.0,1.0);
static const Colour expected_colour(1.0,0.25,0.25);

namespace Ariadne {

template<class X> AffineConstraint<X> operator<=(const Dyadic& l, const Affine<X>& a) {
    return AffineConstraint<X>({l,a.value().precision()},a,{+infty,a.value().precision()}); }
template<class X> AffineConstraint<X> operator<=(const Affine<X>& a, const Dyadic& u) {
    return AffineConstraint<X>(-infty,a,u); }
template<class X> AffineConstraint<X> operator<=(const AffineConstraint<X>& ac, const Dyadic& u) {
    ARIADNE_ASSERT(decide(ac.upper_bound()==infty));
    return AffineConstraint<X>(ac.lower_bound(),ac.function(),{u,ac.lower_bound().precision()}); }
template<class X> AffineConstraint<X> operator==(const Affine<X>& a, const Dyadic& u) {
    return AffineConstraint<X>(u,a,u); }
}

struct FloatDPVector2d : FloatDPVector, Vector2d {
    FloatDPVector2d(double vx, double vy) : FloatDPVector({ExactDouble(vx),ExactDouble(vy)},dp), Vector2d(vx,vy) { }
};

class TestAffineSet
{
  private:
    Figure figure;
    ExactBoxType D;
    Vector< Affine<FloatDPBounds> > x;
    ValidatedAffineConstrainedImageSet set;
  public:
    TestAffineSet() :
        figure(ApproximateBoxType({{-1,+1},{-1,+1}}),Projection2d(2,0,1)),
        set(ExactBoxType::unit_box(2),Matrix<FloatDP>(2,2,dp),Vector<FloatDP>(2,dp)) { }

    Void test_basic_coverage() {
        RealBox real_box(2,RealInterval(-1,1));
        ValidatedAffineConstrainedImageSet from_real(real_box);
        ARIADNE_TEST_EQUALS(from_real.dimension(),2u);
        ARIADNE_TEST_ASSERT(definitely(from_real.is_bounded()));
        ARIADNE_TEST_EQUALS(from_real.number_of_constraints(),0u);

        ExactBoxType exact_box({{-1.0_x,+1.0_x},{-2.0_x,+2.0_x}});
        ValidatedAffineConstrainedImageSet from_exact(exact_box);
        ARIADNE_TEST_EQUALS(from_exact.domain(),exact_box);
        ARIADNE_TEST_ASSERT(definitely(from_exact.inside(ExactBoxType{{-2.0_x,+2.0_x},{-3.0_x,+3.0_x}})));

        Vector<ValidatedAffineModelDP> models=ValidatedAffineModelDP::scalings(ExactBoxType::unit_box(2),dp);
        ValidatedAffineConstrainedImageSet from_models(models);
        List<ValidatedAffineModelConstraintDP> constraints;
        ValidatedAffineConstrainedImageSet from_models_and_constraints(models,constraints);
        ARIADNE_TEST_EQUALS(from_models.number_of_parameters(),2u);
        ARIADNE_TEST_EQUALS(from_models_and_constraints.number_of_constraints(),0u);

        auto y=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction h({y[0]+y[1],y[0]-y[1]});
        auto transformed=image(from_exact,h);
        ARIADNE_TEST_EQUALS(transformed.dimension(),2u);
    }

    Void test_pure_constraint() {
        figure.clear();

        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2,dp);

        set=ValidatedAffineConstrainedImageSet(D, {x[0]+2,x[0]+x[1]} );

        set.new_parameter_constraint(x[0]+0.5_x*x[1]+0.75_x <= 0.0_x);
        set.new_parameter_constraint(-0.5_x*x[0]+1.5_x*x[1]+0.5_x <= 0.0_x);

        ARIADNE_TEST_PRINT(set);
        figure.draw(set);

        set=ValidatedAffineConstrainedImageSet(D, {x[0]+x[1]+2,x[0]+x[1]+2});
        set.new_parameter_constraint(-0.5_x*x[0]-1.5_x*x[1]+0.5_x <= 0.0_x);

        ARIADNE_TEST_PRINT(set);
        figure.draw(set);
    }

    Void test_nonunit_domain() {
        D=ExactBoxType({{0.0_x,2.0_x},{-3.0_x,1.0_x}});
        x=Affine<FloatDPBounds>::variables(2,dp);

        set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]},{x[0]+x[1]<=1.0_x});
        ARIADNE_TEST_EQUALS(set.domain(),D);
        ARIADNE_TEST_EQUALS(set.number_of_constraints(),1u);
        ARIADNE_TEST_SAME(set.bounding_box(),UpperBoxType({{0.0_x,2.0_x},{-3.0_x,1.0_x}}));

        ValidatedAffineConstrainedImageSet unconstrained(D,{x[0],x[1]});
        ARIADNE_TEST_SAME(unconstrained.bounding_box(),UpperBoxType({{0.0_x,2.0_x},{-3.0_x,1.0_x}}));
    }

    Void test_constrained_image() {
        D = ExactBoxType::unit_box(3);
        x = Affine<FloatDPBounds>::variables(3,dp);

        set=ValidatedAffineConstrainedImageSet(D, {2*x[0]+3*x[1]+x[2],x[0]+x[1]-2} );

        set.new_parameter_constraint(1.0_x*x[0]+0.5_x*x[1]+0.25_x*x[2]<=0.75_x);
        set.new_parameter_constraint(-2.0_x*x[0]-3.0_x*x[1]-1.0_x*x[2]<=3.00_x);

        ARIADNE_TEST_PRINT(set)
        Polytope2d expected_set({{-3.0,-3.0}, {-3.0,-3.66667}, {2.0,-2.0}, {4.0,-1.0},
                                     {3.0,-0.5}, {1.0,-1.0}, {-2.0,-2.0}, {-3.0,-2.5}});
        ARIADNE_TEST_PRINT(expected_set)

        figure.clear(); figure.set_bounding_box(ExactBoxType{{-7.0_x,+7.0_x},{-5.0_x,+1.0_x}}); figure.set_fill_opacity(0.05);
        figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        figure.write("test_affine_set-constrained_image");
    }


    Void test_separated() {
        D=ExactBoxType::unit_box(3);
        x=Affine<FloatDPBounds>::variables(3,dp);
        set=ValidatedAffineConstrainedImageSet(D, {-0.9375_x+0.0625_x*x[0]+0.5_x*x[1],0.87890625_x-0.1171875_x*x[0]+x[1]+0.00390625_x*x[2]});
        set.new_parameter_constraint(-1.1875_x+0.0625_x*x[0]+x[1]<=0.0_x);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell1,({{-1.0_x,-0.9375_x},{0.8125_x, 0.875_x}})); // subset
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell2,({{-0.5625_x,-0.50_x},{1.4375_x,1.5_x}})); // overlaps
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell3,({{-0.875_x,-0.75_x},{0.625_x,0.7578125_x}})); // touches at (-0.875_x,0.7578125_x);
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell4,({{-1.1875_x,-1.125_x},{0.0625_x,0.125_x}})); // almost touches
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell5,({{-0.9375_x,-0.875_x},{0.4375_x,0.5_x}})); // disjoint
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell6,({{-1.5_x,-1.375_x},{0.5_x,0.625_x}})); // disjoint; regression test

        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell1)));
        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell2)));
        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell3)));
        //ARIADNE_TEST_ASSERT(!definitely(affine_set.overlaps(cell3)));
        ARIADNE_TEST_ASSERT(possibly(set.separated(cell4)));
        ARIADNE_TEST_ASSERT(definitely(set.separated(cell5)));
        ARIADNE_TEST_ASSERT(definitely(set.separated(cell6)));

        figure.clear();
        figure.set_bounding_box(widen(set.bounding_box(),0.125_x));
        figure << set
               << fill_colour(0,0,1) << cell1
               << fill_colour(0,1,0) << cell2
               << fill_colour(1,1,0) << cell3
               << fill_colour(1,0,0) << cell4
               << fill_colour(1,0,0) << cell5
               << fill_colour(1,0,0) << cell6
               ;
        figure.write("test_affine_set-separated");
    }

    Void test_robust_outer_approximation() {
        D=ExactBoxType::unit_box(2);
        x=Affine<FloatDPBounds>::variables(2,dp);
        set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
        set.new_parameter_constraint(x[0]+x[1]<=0.5_x);

        GridTreePaving paving(Grid(2));
        set.robust_adjoin_outer_approximation_to(paving,0);
        ARIADNE_TEST_ASSERT(paving.size()>0);

        // Repeating the ordinary approximation must recognise cells already
        // contained in the paving and leave the result unchanged.
        SizeType size=paving.size();
        set.adjoin_outer_approximation_to(paving,0);
        ARIADNE_TEST_EQUALS(paving.size(),size);

        ValidatedAffineConstrainedImageSet empty_set(D,{x[0],x[1]});
        empty_set.new_parameter_constraint(x[0]<=-2.0_x);
        GridTreePaving empty_paving(Grid(2));
        empty_set.robust_adjoin_outer_approximation_to(empty_paving,0);
        ARIADNE_TEST_EQUALS(empty_paving.size(),0u);

        ValidatedAffineConstrainedImageSet infeasible_set(D,{x[0],x[1]});
        infeasible_set.new_parameter_constraint(x[0]<=-2.0_x);
        infeasible_set.new_parameter_constraint(-x[0]<=-2.0_x);
        GridTreePaving infeasible_paving(Grid(2));
        infeasible_set.robust_adjoin_outer_approximation_to(infeasible_paving,0);
        ARIADNE_TEST_EQUALS(infeasible_paving.size(),0u);

        ExactBoxType small_box({{0.25_x,0.75_x},{0.25_x,0.75_x}});
        ValidatedAffineConstrainedImageSet small_set(small_box);
        GridTreePaving covered_paving(Grid(2));
        small_set.robust_adjoin_outer_approximation_to(covered_paving,0);
        ARIADNE_TEST_ASSERT(covered_paving.size()>0);
        SizeType covered_size=covered_paving.size();
        small_set.robust_adjoin_outer_approximation_to(covered_paving,0);
        ARIADNE_TEST_EQUALS(covered_paving.size(),covered_size);

    }

    Void test_outer_approximation() {

        {
            D=ExactBoxType::unit_box(2);
            x=Affine<FloatDPBounds>::variables(2,dp);
            set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
            set.new_parameter_constraint(x[0]+x[1]<=-0.5_x);
            ARIADNE_TEST_PRINT(set);

            GridTreePaving paving = set.outer_approximation(Grid(2),4);
            paving.recombine();
            std::cout<<std::setprecision(17);
            ARIADNE_TEST_PRINT(paving);

            figure.clear();
            figure.set_bounding_box(widen(set.bounding_box(),0.125_x));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-4");
        }

        {
            // The following set has difficulties
            set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
            set.new_parameter_constraint(2.0_x*x[0]+1.0_x*x[1]<=0.5_x);
            set.new_parameter_constraint(-0.5_x*x[0]+1.0_x*x[1]<=0.75_x);
            set.new_parameter_constraint(-0.5_x*x[0]-1.0_x*x[1]<=0.875_x);

            figure.clear();
            figure.set_bounding_box(widen(ExactBoxType{{-1.0_x,+1.0_x},{-1.0_x,+1.0_x}},+0.125_x));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(set.outer_approximation(Grid(2),3));
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);

            figure.write("test_affine_set-outer_approximation-2");
            figure.clear();
        }

        {
            D=ExactBoxType::unit_box(3);
            x=Affine<FloatDPBounds>::variables(3,dp);

            set=ValidatedAffineConstrainedImageSet( D, {2.0_pr*x[0]+3.0_pr*x[1]+1.0_pr*x[2]+0.05_pr, x[0]+x[1]+2.051_pr} );

            set.new_parameter_constraint(+1.01_pr*x[0]+0.51_pr*x[1]+0.251_pr*x[2]<=0.751_pr);
            set.new_parameter_constraint(-0.51_pr*x[0]+1.51_pr*x[1]+0.01_pr*x[2]<=0.51_pr);
            set.new_parameter_constraint(-2.01_pr*x[0]-3.01_pr*x[1]-1.02_pr*x[2]<=3.01_pr);

            Grid grid(2);
            GridTreePaving paving(grid);
            set.adjoin_outer_approximation_to(paving,3);
            paving.recombine();

            figure.clear();
            figure.set_bounding_box(ExactBoxType{{-7.0_x,+7.0_x},{-1.0_x,+5.0_x}});
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.set_fill_opacity(0.5);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-1");
            figure.clear();
        }

        {
            set=ValidatedAffineConstrainedImageSet(ExactBoxType::unit_box(3),{-0.9375_x+0.0625_x*x[0]+0.5_x*x[1],0.87890625_x-0.1171875_x*x[0]+x[1]+0.00390625_x*x[2]});
            set.new_parameter_constraint(0.0625_x*x[0]+x[1]<=1.1875_x);
            ARIADNE_TEST_PRINT(set);

            GridTreePaving paving = set.outer_approximation(Grid(2),4);
            paving.recombine();
            std::cout<<std::setprecision(17);
            ARIADNE_TEST_PRINT(paving);

            figure.clear();
            figure.set_bounding_box(widen(set.bounding_box(),+0.125_x));
            figure.set_fill_opacity(0.5);
            figure.set_fill_colour(1.0,0.0,0.0);
            figure.draw(paving);
            figure.set_fill_colour(0.0,0.5,1.0);
            figure.draw(set);
            figure.write("test_affine_set-outer_approximation-3");
        }


    }

    Void test_empty() {
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2,dp);

        std::cerr<<get_characteristics(0.125_x*x[0]+1.0_x)<<"\n";
        set=ValidatedAffineConstrainedImageSet(D, {0.125_x*x[0]+1.0_x, 0.125_x*x[1]+1.0_x} );
        set.new_parameter_constraint(-1.0_x*x[1] <= -2.0_x);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.is_empty());

        figure.clear();
        figure.draw(set);
        figure.write("test_affine_set-empty");
    }

    Void test_invalid() {
        {
            D = ExactBoxType::unit_box(1);
            x = Affine<FloatDPBounds>::variables(2,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D, {0.125_x*x[0]+1.0_x, 0.125_x*x[1]+1.0_x} ));
        }
        {
            D = ExactBoxType::unit_box(2);
            x = Affine<FloatDPBounds>::variables(2,dp);
            auto a = Affine<FloatDPBounds>::variables(1,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D, {0.125_x*x[0]+1.0_x, 0.125_x*x[1]+1.0_x}, {a[0]<=-0.25_x} ));

            set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
            ARIADNE_TEST_FAIL(set.new_parameter_constraint(a[0]<=-0.25_x));
            ARIADNE_TEST_FAIL(set.separated(ExactBoxType::unit_box(1)));
            ARIADNE_TEST_FAIL(set.inside(ExactBoxType::unit_box(1)));

            GridTreePaving bad_paving(Grid(1));
            ARIADNE_TEST_FAIL(set.adjoin_outer_approximation_to(bad_paving,0));
            ARIADNE_TEST_FAIL(set.robust_adjoin_outer_approximation_to(bad_paving,0));

            Vector<ValidatedAffineModelDP> bad_models=ValidatedAffineModelDP::scalings(ExactBoxType::unit_box(1),dp);
            List<ValidatedAffineModelConstraintDP> no_constraints;
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,bad_models,no_constraints));

            ValidatedAffineModelConstraintDP bad_constraint(-1.0_x,bad_models[0],+1.0_x);
            List<ValidatedAffineModelConstraintDP> bad_constraints;
            bad_constraints.append(bad_constraint);
            Vector<ValidatedAffineModelDP> good_models=ValidatedAffineModelDP::scalings(D,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,good_models,bad_constraints));
            ARIADNE_TEST_FAIL(set.new_constraint(bad_constraint));

            ValidatedAffineConstrainedImageSet empty_set(ExactBoxType(0u));
            ValidatedAffineModelDP zero_model(0u,dp);
            ValidatedAffineModelConstraintDP zero_model_constraint(-1.0_x,zero_model,+1.0_x);
            Affine<FloatDPBounds> zero_affine(0u,dp);
            ARIADNE_TEST_FAIL(empty_set.new_constraint(zero_model_constraint));
            ARIADNE_TEST_FAIL(empty_set.new_parameter_constraint(zero_affine<=0.0_x));
            ARIADNE_TEST_FAIL(empty_set.number_of_parameters());

            Vector<Affine<FloatDPBounds>> empty_affine_models;
            Vector<ValidatedAffineModelDP> empty_models;
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,empty_affine_models));
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,empty_models,no_constraints));
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet{empty_models});

            Matrix<FloatDP> bad_matrix(1u,1u,dp);
            Vector<FloatDP> bad_offset(2u,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,bad_matrix,bad_offset));
            Matrix<FloatDP> empty_matrix(0u,0u,dp);
            Vector<FloatDP> empty_offset(0u,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(ExactBoxType(0u),empty_matrix,empty_offset));
        }
    }

    Void test_large_zonotope_boundary() {
        SizeType const n=501u;
        D=ExactBoxType::unit_box(n);

        Affine<FloatDPBounds> fx(n,dp);
        Affine<FloatDPBounds> fy(n,dp);
        for(SizeType i=0; i!=n; ++i) {
            Nat j=static_cast<Nat>(i);
            fx[j]=FloatDPBounds(1,dp);
            fy[j]=FloatDPBounds(static_cast<Int>(i+1u),dp);
        }

        set=ValidatedAffineConstrainedImageSet(D,{fx,fy});
        List<Point2d> vertices=set.boundary(0,1);
        ARIADNE_TEST_ASSERT(vertices.size()>1000u);
    }

    Void test_logging_branches() {
        unsigned int const old_verbosity=Logger::instance().configuration().verbosity();
        Logger::instance().configuration().set_verbosity(10u);

        D=ExactBoxType::unit_box(2);
        x=Affine<FloatDPBounds>::variables(2,dp);
        set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});

        ARIADNE_TEST_ASSERT(!definitely(set.separated(ExactBoxType{{-0.5_x,+0.5_x},{-0.5_x,+0.5_x}})));

        GridTreePaving paving(Grid(2));
        set.robust_adjoin_outer_approximation_to(paving,0);
        ARIADNE_TEST_ASSERT(paving.size()>0);

        List<Point2d> vertices=set.boundary(0,1);
        ARIADNE_TEST_ASSERT(!vertices.empty());
        figure.draw(set);
        figure.write("test_affine_set-logging");

        Logger::instance().configuration().set_verbosity(old_verbosity);
    }

    Void test_uniform_error() {
        FloatDP e(0.25_x,dp);
        FloatDPBounds E = FloatDPBounds(-e,+e);
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2,dp);
        ValidatedAffineConstrainedImageSet set1=ValidatedAffineConstrainedImageSet(D, {0.5_x*x[0]+0.25_x*x[1]+E/2, 0.25_x*x[0]-0.25_x*x[1]+E} );
        ARIADNE_TEST_PRINT(set1);
        D = ExactBoxType::unit_box(4);
        x = Affine<FloatDPBounds>::variables(4,dp);
        ValidatedAffineConstrainedImageSet set2=ValidatedAffineConstrainedImageSet(D, {0.5_x*x[0]+0.25_x*x[1]+e*x[2]/2, 0.25_x*x[0]-0.25_x*x[1]+e*x[3]} );
        ARIADNE_TEST_PRINT(set2);

        figure.clear();
        figure.set_bounding_box(ExactBoxType{{-1.0_x,+1.0_x},{-1.0_x,+1.0_x}});
        figure.set_fill_opacity(0.5);
        figure.set_fill_colour(0,1,0);
        figure.draw(set2);
        figure.set_fill_colour(0,0,1);
        figure.draw(set1);
        figure.write("test_affine_set-uniform_error");

        Grid g(2);
        GridTreePaving paving1=outer_approximation(set1,g,4);
        GridTreePaving paving2=outer_approximation(set2,g,4);

        figure.clear();
        figure.set_fill_colour(1,0,0);
        figure.draw(paving1);
        figure.set_fill_colour(0,0,1);
        figure.draw(set1);
        figure.write("test_affine_set-uniform_error-paving");

        ARIADNE_TEST_ASSERT(paving1==paving2);

    }


    Void test_draw() {
        FloatDPVector2d offsets{0.0,0.0}; // Offsets
        FloatDPVector& o=offsets;
        Polytope2d expected_set({{0.0,0.0}});
        Vector< Affine<FloatDPBounds> > a;
        ExactBoxType dom;
        figure.clear();
        figure.set_bounding_box(ExactBoxType{{-1.0_x,+15.0_x},{-1.0_x,+15.0_x}});
        figure.set_fill_colour(0.5,1.0,1.0);
        figure.set_fill_opacity(0.5);

        {
            // Draw overlapping sets to check colours
            figure << fill_colour(expected_colour) << Polytope2d({{0.,-0.5},{2.,-0.5},{2.,0.},{0.,0.}})
                   << fill_colour(colour) << Polytope2d({{1.,-0.5},{3.,-0.5},{3.,0.},{1.,0.}});
        }

        {
            // Test draw of nondegenerate two-dimensional image set
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{1.0,13.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_x*a[0],o[1]+0.25_x*a[1]});
            expected_set=Polytope2d({{-0.5,-0.25}, {+0.5,-0.25}, {+0.5,+0.25},{-0.5,+0.25}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of nondegenerate two-dimensional constrained image set
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d(4.0,13.0);
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_x*a[0],o[1]+0.25_x*a[1]},{a[0]+a[1]<=0.5_x});
            expected_set=Polytope2d({{-0.5,-0.25}, {+0.5,-0.25}, {+0.5,-0.125}, {-0.25,+0.25},{-0.5,+0.25}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of two-dimensional image set with errors
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{7.0,13.0};
            FloatDPBounds e=FloatDPBounds(-1.0_x,+1.0_x,dp);

            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.25_x*a[0]+0.25_x*a[1]+0.5_x*e,o[1]+0.5_x*a[0]-0.5_x*a[1]});
            expected_set=Polytope2d({{-1.0,0.0}, {-0.5,-1.0}, {+0.5,-1.0}, {+1.0,0.0}, {+0.5,+1.0}, {-0.5,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of nondegenerate one-dimensional image set
            a=Affine<FloatDPBounds>::variables(1,dp);
            dom=ExactBoxType::unit_box(1);
            offsets=FloatDPVector2d{1.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_x*a[0],o[1]+0.25_x*a[0]},{0.0_x*a[0]<=1.0_x,a[0]==0.75_x});
            expected_set=Polytope2d({{0.375, 0.1875}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional image set
            a=Affine<FloatDPBounds>::variables(1,dp);
            dom=ExactBoxType::unit_box(1);
            offsets=FloatDPVector2d{4.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.5_x*a[0],o[1]+0.25_x*a[0]});
            expected_set=Polytope2d({{-0.5,-0.25}, {0.5,0.25}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional constraint set
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{7.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{0.0_x*a[0]<=1.0_x,a[0]+a[1]==0.5_x});
            expected_set=Polytope2d({{-0.5,+1.0}, {+1.0,-0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of one-dimensional constraint set with one proper constraint on boundary
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{10.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{1.0_x<=a[0], a[1]<=0.5_x});
            expected_set=Polytope2d({{+1.0,-1.0}, {+1.0,+0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with constraint a<=0 and a>=0
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{13.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+0.5_x*a[1]<=0.75_x,0.75_x<=a[0]+0.5_x*a[1]});
            expected_set=Polytope2d({{+0.25,+1.0}, {+1.0,-0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with constraint 0<=a<=0
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{13.0,1.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{0.75_x<=a[0]+0.5_x*a[1]<=0.75_x});
            expected_set=Polytope2d({{+0.25,+1.0}, {+1.0,-0.5}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with degenerate constraints at corner of box
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{1.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=3,1.5_x*(a[0]+a[1])<=3,2*a[0]+a[1]<=3});
            expected_set=Polytope2d({{-1.0,-1.0}, {+1.0,-1.0}, {+1.0,+1.0}, {-1.0,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set with degenerate constraints near corner of box
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{4.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2*a[1]<=2,1.5_x*(a[0]+a[1])<=2,2*a[0]+a[1]<=2});
            expected_set=Polytope2d({{-1.0,-1.0}, {+1.0,-1.0}, {+1.0,0.0}, {+0.667,+0.667}, {0.0,+1.0}, {-1.0,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of set repeated constraints
            a=Affine<FloatDPBounds>::variables(2,dp);
            dom=ExactBoxType::unit_box(2);
            offsets=FloatDPVector2d{7.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2.0_x*a[1]<=2.0_x,a[0]/3.0_x+a[1]*2.0_x/3.0_x<=0.333333333333333333_pr});
            expected_set=Polytope2d({{-1.0,-1.0}, {+1.0,-1.0}, {+1.0,0.5}, {0.0,+1.0}, {-1.0,+1.0}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of projection of three-dimensional box with single equality constraint
            a=Affine<FloatDPBounds>::variables(3,dp);
            dom=ExactBoxType::unit_box(3);
            offsets=FloatDPVector2d{10.0,4.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+a[0],o[1]+a[1]},{a[0]+2.0_x*a[1]+a[2]==1.5_x});
            expected_set=Polytope2d({{+1.0,-0.25}, {+1.0,+0.75}, {+0.5,+1.0}, {-1.0,+1.0}, {-1.0,+0.75}}) + offsets;
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        {
            // Test draw of two-dimensional set with nondegenerate inequality and equality constraints
            a=Affine<FloatDPBounds>::variables(3,dp);
            dom=ExactBoxType::unit_box(3);
            offsets=FloatDPVector2d{1.0,7.0};
            set=ValidatedAffineConstrainedImageSet(dom, {o[0]+0.3_pr*a[0]+0.20_pr*a[1]+0.05_pr*a[2],o[1]-0.10_pr*a[0]+0.1_pr*a[1]+0.05_pr*a[2]},{a[1]-a[2]<=-0.25_x,a[0]+a[1]+a[2]==0.5_x});
            expected_set=Polytope2d({{+0.0,0.0}}) + offsets; // Unknown
            figure << fill_colour(expected_colour) << expected_set << fill_colour(colour) << set;
        }

        figure.write("test_affine_set-draw");
    }

    Void test() {
        figure.set_bounding_box(ExactBoxType{{-4.0_x,+4.0_x},{-4.0_x,+4.0_x}});
        ARIADNE_TEST_CALL(test_basic_coverage());
        ARIADNE_TEST_CALL(test_empty());
        ARIADNE_TEST_CALL(test_invalid());
        ARIADNE_TEST_CALL(test_pure_constraint());
        ARIADNE_TEST_CALL(test_constrained_image());
        ARIADNE_TEST_CALL(test_nonunit_domain());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_outer_approximation());
        ARIADNE_TEST_CALL(test_robust_outer_approximation());
        ARIADNE_TEST_CALL(test_large_zonotope_boundary());
        ARIADNE_TEST_CALL(test_logging_branches());
        ARIADNE_TEST_CALL(test_uniform_error());
        ARIADNE_TEST_CALL(test_draw());
    }

};



Int main()
{
    TestAffineSet().test();
    std::cerr<<"INCOMPLETE ";
    return ARIADNE_TEST_FAILURES;
}
