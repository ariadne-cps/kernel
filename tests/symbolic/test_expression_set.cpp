/***************************************************************************
 *            test_expression_set.cpp
 *
 *  Copyright  2009-21  Luca Geretti
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
#include <sstream>
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/valuation.hpp"
#include "symbolic/space.hpp"
#include "symbolic/formula.hpp"
#include "geometry/function_set.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "io/command_line_interface.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;

using RealVariablesPoint = VariablesPoint<Real>;

class TestExpressionSet {
  private:
    RealVariable x, y, z;

  public:
    TestExpressionSet() : x("x"), y("y"), z("z") { }

    Void test() {
        ARIADNE_TEST_CALL(test_labelled_set());
        ARIADNE_TEST_CALL(test_variables_point());
        ARIADNE_TEST_CALL(test_template_streaming());
        ARIADNE_TEST_CALL(test_variables_box());
        ARIADNE_TEST_CALL(test_real_expression_constraint_set());
        ARIADNE_TEST_CALL(test_real_expression_bounded_constraint_set());
    }

    Void test_labelled_set() const {
        LabelledSet<RealBox> set(RealSpace({x,y,z}),RealBox({{0.1_dec,0.5_dec},{-2.3_dec,1},{0,3}}));
        ARIADNE_TEST_EQUALS(set.dimension(),3);
        RealVariable u("u");
        ARIADNE_TEST_EXECUTE(set.euclidean_set(RealSpace({y,z})));
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({y,u})));
    }

    Void test_variables_point() const {
        Map<RealVariable,Real> bindings;
        bindings.insert(x,0.1_dec);
        bindings.insert(y,-1);
        RealVariablesPoint p(bindings);
        ARIADNE_TEST_EQUALS(p.variables().size(),2);
        ARIADNE_TEST_EQUALS(p.euclidean_set(RealSpace({x,y})).dimension(),2);
        ARIADNE_TEST_FAIL(p.euclidean_set(RealSpace({x,z})));

        auto labelled_set = LabelledSet<RealPoint>(p);
        ARIADNE_TEST_EQUALS(labelled_set.dimension(),2);
    }

    Void test_template_streaming() const {
        BoundedConstraintSet bounded_set(RealBox({{-1,1}}));
        ExactBoxType exact_box({ExactIntervalType(-1,1)});
        Symbolic<Separated,BoundedConstraintSet,ExactBoxType> separated_symbolic(
            Separated(),bounded_set,exact_box);
        std::ostringstream separated_stream;
        separated_stream << separated_symbolic;
        ARIADNE_TEST_ASSERT(not separated_stream.str().empty());
    }

    Void test_variables_box() const {
        Real real_upper(2);
        auto real_upper_interval = real_upper >= x;
        ARIADNE_TEST_ASSERT(real_upper_interval.variable()==x);
        ARIADNE_TEST_EQUAL(real_upper_interval.upper_bound(),real_upper);

        Real real_value(1);
        auto variable_equal_real = x == real_value;
        auto real_equal_variable = real_value == x;
        ARIADNE_TEST_ASSERT(variable_equal_real.variable()==x);
        ARIADNE_TEST_EQUAL(variable_equal_real.lower_bound(),real_value);
        ARIADNE_TEST_EQUAL(variable_equal_real.upper_bound(),real_value);
        ARIADNE_TEST_ASSERT(real_equal_variable.variable()==x);
        ARIADNE_TEST_EQUAL(real_equal_variable.lower_bound(),real_value);
        ARIADNE_TEST_EQUAL(real_equal_variable.upper_bound(),real_value);

        FloatDP float_lower(0,dp);
        auto float_lower_interval = float_lower <= x;
        ARIADNE_TEST_ASSERT(float_lower_interval.variable()==x);
        ARIADNE_TEST_EQUAL(float_lower_interval.lower_bound(),float_lower);

        Dyadic dyadic_lower(1);
        auto dyadic_lower_interval = dyadic_lower <= x;
        ARIADNE_TEST_ASSERT(dyadic_lower_interval.variable()==x);
        ARIADNE_TEST_EQUAL(dyadic_lower_interval.lower_bound(),dyadic_lower);

        double double_lower=0.25;
        auto double_lower_interval = double_lower <= x;
        ARIADNE_TEST_ASSERT(double_lower_interval.variable()==x);
        ARIADNE_TEST_EQUAL(double_lower_interval.lower_bound(),double_lower);

        RealVariablesBox bx1 = {0<=x<=1.2_dec};
        ARIADNE_TEST_ASSERT(not bx1.is_empty());
        ARIADNE_TEST_PRINT(bx1);

        RealVariablesBox bx2(RealSpace({x,y}),RealBox({{0.1_dec,0.5_dec},{-2.3_dec,1}}));
        ARIADNE_TEST_EQUALS(bx2.variables().size(),2);
        ARIADNE_TEST_PRINT(bx2);

        auto labelled_set = LabelledSet<RealBox>(bx2);
        ARIADNE_TEST_EQUALS(labelled_set.dimension(),2);

        ARIADNE_TEST_EQUALS(variables(bx2.bounds()).size(),2);

        ARIADNE_TEST_PRINT(over_approximation(bx2));
        ARIADNE_TEST_PRINT(approximation(bx2));
        ARIADNE_TEST_PRINT(under_approximation(bx2));

        ARIADNE_TEST_EXECUTE(make_set(RealSpace({y,x}),bx2));
        ARIADNE_TEST_EXECUTE(make_box(RealSpace({x}),bx2));
        ARIADNE_TEST_FAIL(make_box(RealSpace({x,z}),bx2));
    }

    Void test_real_expression_constraint_set() const {
        RealExpressionConstraintSet set({y >= 0, z + y <= 2, sqr(z) + sqr(y) >= 0});
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_EQUALS(set.variables().size(), 2);
        ARIADNE_TEST_EQUALS(set.constraints().size(), 3);
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({y})));
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({x})));
        ARIADNE_TEST_EXECUTE(set.euclidean_set(RealSpace({y, z})));
        ARIADNE_TEST_EXECUTE(make_set(RealSpace({z, y}),set));
    }

    Void test_real_expression_bounded_constraint_set() const {
        RealExpressionConstraintSet unbounded_set({y>=0,z+y<=2,sqr(z)+sqr(y)>=0});
        ARIADNE_TEST_PRINT(unbounded_set);
        RealVariable u("u");

        List<RealVariableInterval> bounds_list={0.1_dec<=y<=0.5_dec,-2.3_dec<=z<=1};
        List<ContinuousPredicate> constraints_list={y>=0,z+y<=2};
        ARIADNE_TEST_EQUALS(variables(bounds_list).size(),2u);

        RealExpressionBoundedConstraintSet from_initializer({0.1_dec<=y<=0.5_dec,-2.3_dec<=z<=1});
        ARIADNE_TEST_EQUALS(from_initializer.variables().size(),2u);

        RealExpressionBoundedConstraintSet from_initializer_constraints(
            {0.1_dec<=y<=0.5_dec,-2.3_dec<=z<=1},{y>=0,z+y<=2});
        ARIADNE_TEST_EQUALS(from_initializer_constraints.constraints().size(),2u);

        RealExpressionBoundedConstraintSet from_list(bounds_list);
        ARIADNE_TEST_EQUALS(from_list.bounds().size(),2u);

        RealExpressionBoundedConstraintSet from_list_constraints(bounds_list,constraints_list);
        ARIADNE_TEST_EQUALS(from_list_constraints.constraints().size(),2u);

        Map<RealVariable,RealInterval> bounds_map=RealVariablesBox(bounds_list).bounds();
        RealExpressionBoundedConstraintSet from_map(bounds_map);
        ARIADNE_TEST_EQUALS(from_map.variables().size(),2u);

        RealVariablesBox bounds_box(bounds_list);
        RealExpressionBoundedConstraintSet from_box(bounds_box);
        ARIADNE_TEST_EQUALS(from_box.variables().size(),2u);

        List<RealVariableInterval> duplicate_bounds={0.1_dec<=y<=0.5_dec,-1<=y<=1};
        ARIADNE_TEST_FAIL((void)RealExpressionBoundedConstraintSet(duplicate_bounds));
        List<ContinuousPredicate> duplicate_compatible_constraints={y>=0};
        ARIADNE_TEST_FAIL((void)RealExpressionBoundedConstraintSet(duplicate_bounds,duplicate_compatible_constraints));
        ARIADNE_TEST_FAIL(RealExpressionBoundedConstraintSet(RealVariablesBox({0.1_dec<=x<=0.5_dec,-2.3_dec<=u<=1}),unbounded_set));
        ARIADNE_TEST_EXECUTE(intersection(RealVariablesBox({0.1_dec<=z<=0.5_dec,-2.3_dec<=y<=1}),unbounded_set));
        ARIADNE_TEST_EXECUTE(make_set(RealSpace({y,z}),RealVariablesBox({0.1_dec<=z<=0.5_dec,-2.3_dec<=y<=1}),unbounded_set));

        RealExpressionBoundedConstraintSet set(RealVariablesBox({0.1_dec<=y<=0.5_dec,-2.3_dec<=z<=1}),unbounded_set);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_FAIL(approximate_euclidean_set(set,RealSpace({x})));
        ARIADNE_TEST_EXECUTE(approximate_euclidean_set(set,RealSpace({y,z})));
        ARIADNE_TEST_FAIL(approximate_euclidean_set(set,RealSpace({x,y,z})));
        ARIADNE_TEST_FAIL(make_set(RealSpace({y,z,u}),set));

        RealExpressionBoundedConstraintSet set2(RealVariablesBox({0.3_dec<=x<=1,0.1_dec<=y<=0.5_dec,-2.3_dec<=z<=1}),unbounded_set);
        ARIADNE_TEST_PRINT(set2);
        ARIADNE_TEST_FAIL(approximate_euclidean_set(set2,RealSpace({x})));
        ARIADNE_TEST_EXECUTE(approximate_euclidean_set(set2,RealSpace({y,z})));
        ARIADNE_TEST_EXECUTE(approximate_euclidean_set(set2,RealSpace({x,y,z})));
        ARIADNE_TEST_FAIL(make_set(RealSpace({y,z,u}),set2));
    }
};

Int main(Int argc, const char* argv[])
{
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;

    TestExpressionSet().test();
    return ARIADNE_TEST_FAILURES;
}
