/***************************************************************************
 *            test_symbolic_templates.cpp
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

#include <sstream>
#include "numeric/numeric.hpp"
#include "algebra/algebra.hpp"
#include "function/formula.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/valuation.hpp"
#include "symbolic/templates.tpl.hpp"
#include "symbolic/expression.tpl.hpp"
#include "macro/testing.hpp"

namespace Ariadne {
extern template class Expression<Real>;
extern template class Expression<Boolean>;
extern template class Expression<String>;
}

using namespace Ariadne;

class TestSymbolicTemplates {
    RealVariable x{"x"}, y{"y"}, z{"z"};
  public:
    Void test_symbolic_template_operations() {
        RealExpression ex(x), ey(y), zero(0);
        ARIADNE_TEST_ASSERT(identical(simplify(zero/ex),zero));
        ARIADNE_TEST_ASSERT(is_affine_in(pos(ex),{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(pos(sqr(ex)),{x}));

        // Construct the supported Nul node explicitly: nul(ex) eagerly returns a constant.
        RealExpression null_node=make_expression<Real>(Nul(),sqr(ex));
        ARIADNE_TEST_ASSERT(is_affine_in(null_node,{x}));
        ARIADNE_TEST_ASSERT(is_polynomial_in(null_node,{x}));
        RealValuation values;
        values.set(x,Real(3));
        ARIADNE_TEST_EQUAL(evaluate(null_node,values),Real(0));

        std::ostringstream product_stream, quotient_stream;
        product_stream << OperatorExpressionWriter()(ex*(ey-RealExpression(z)));
        quotient_stream << OperatorExpressionWriter()(ex/(ey+RealExpression(z)));
        ARIADNE_TEST_EQUALS(product_stream.str(),"x*(y-z)");
        ARIADNE_TEST_EQUALS(quotient_stream.str(),"x/(y+z)");

        BooleanVariable p("p");
        BooleanExpression ep(p);
        std::ostringstream logical_stream;
        logical_stream << OperatorExpressionWriter()(!ep);
        ARIADNE_TEST_EQUALS(logical_stream.str(),"!(p)");
        StringVariable label("label");
        StringValuation labels;
        labels.set(label,String("ok"));
        ARIADNE_TEST_ASSERT(evaluate(!(label==String("ok")),labels)==Boolean(false));
        ARIADNE_TEST_ASSERT(evaluate(!(label==String("other")),labels)==Boolean(true));

        RealVectorExpression vector({x,y});
        RealVectorExpression constants(RealVector({1,2}));
        // Same operator code, different node types: structural identity is not commutativity.
        RealVectorExpression left_scaled=RealExpression(2)*vector;
        RealVectorExpression right_scaled=vector*RealExpression(2);
        ARIADNE_TEST_ASSERT(not identical(left_scaled,right_scaled));
        ARIADNE_TEST_ASSERT(not identical(right_scaled,left_scaled));
        values.set(y,Real(4));
        ARIADNE_TEST_EQUAL(evaluate(left_scaled,values)[0],Real(6));
        ARIADNE_TEST_EQUAL(evaluate(right_scaled,values)[1],Real(8));
        ARIADNE_TEST_ASSERT(not identical(vector,constants));
        ARIADNE_TEST_ASSERT(not identical(constants,vector));
        ARIADNE_TEST_ASSERT(is_constant_in(neg(vector),{z}));
        ARIADNE_TEST_ASSERT(not is_constant_in(neg(vector),{x}));
        ARIADNE_TEST_ASSERT(is_constant_in(vector+constants,{z}));
        ARIADNE_TEST_ASSERT(not is_constant_in(vector+constants,{x}));
        ARIADNE_TEST_ASSERT(not is_constant_in(constants+vector,{x}));
        ARIADNE_TEST_ASSERT(is_constant_in(vector-constants,{z}));
        ARIADNE_TEST_ASSERT(is_constant_in(RealExpression(2)*vector,{z}));
        ARIADNE_TEST_ASSERT(not is_constant_in(ex*vector,{x}));
        ARIADNE_TEST_ASSERT(not is_constant_in(vector/ex,{x}));
        ARIADNE_TEST_ASSERT(is_affine_in(vector[0],{x,y}));
        RealVectorExpression nonlinear({sqr(ex),ey});
        ARIADNE_TEST_ASSERT(not is_affine_in(nonlinear[0],{x,y}));
    }

};

Int main() {
    TestSymbolicTemplates().test_symbolic_template_operations();
    return ARIADNE_TEST_FAILURES;
}
