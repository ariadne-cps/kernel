/***************************************************************************
 *            test_validated_constrained_image_set.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "geometry/function_set.hpp"
#include "geometry/affine_set.hpp"
#include "geometry/grid_paving.hpp"
#include "function/taylor_function.hpp"
#include "drawing/null_graphics.hpp"
#include "macro/testing.hpp"
#include <memory>
#include <sstream>
using namespace Ariadne;

Void test_default_image_sets() {
    ConstrainedImageSet effective;
    ValidatedConstrainedImageSet validated;
    ARIADNE_TEST_EQUALS(effective.dimension(),0u);
    ARIADNE_TEST_EQUALS(validated.dimension(),0u);
    ARIADNE_TEST_EQUALS(effective.number_of_parameters(),0u);
    ARIADNE_TEST_EQUALS(validated.number_of_parameters(),0u);
    ARIADNE_TEST_EQUALS(effective.number_of_constraints(),0u);
    ARIADNE_TEST_EQUALS(validated.number_of_constraints(),0u);
    ARIADNE_TEST_EQUALS(effective.bounding_box().dimension(),0u);
    ARIADNE_TEST_EQUALS(validated.bounding_box().dimension(),0u);
}

Void test_validated_image_operations() {
    auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction mapping({x+y,x-y});
    ExactBoxType domain{{-1.0_x,1.0_x},{-2.0_x,2.0_x}};
    ValidatedConstrainedImageSet set(domain,mapping);
    set.new_parameter_constraint(x>=0);
    set.new_space_constraint(x<=ValidatedNumber(2));
    ARIADNE_TEST_EQUALS(set.constraints().size(),2u);
    ARIADNE_TEST_EQUALS(set.constraint(0u).function().argument_size(),2u);
    ARIADNE_TEST_EQUALS(set.constraint_function().result_size(),2u);
    ARIADNE_TEST_EQUALS(set.constraint_bounds().dimension(),2u);
    ARIADNE_TEST_EQUALS(set.number_of_parameters(),2u);
    ARIADNE_TEST_EQUALS(set.dimension(),2u);
    std::unique_ptr<ValidatedConstrainedImageSet> copy(set.clone());
    ARIADNE_TEST_EQUALS(copy->number_of_constraints(),2u);
    auto pieces=set.split();
    ARIADNE_TEST_EQUALS(pieces.first.number_of_constraints(),2u);
    ARIADNE_TEST_EQUALS(pieces.second.number_of_constraints(),2u);
    ARIADNE_TEST_ASSERT(same(pieces.first.domain(),ExactBoxType{{-1.0_x,1.0_x},{-2.0_x,0.0_x}}));
    auto explicit_pieces=set.split(0u);
    ARIADNE_TEST_ASSERT(same(explicit_pieces.second.domain(),ExactBoxType{{0.0_x,1.0_x},{-2.0_x,2.0_x}}));
    auto restricted=set.restriction(ExactBoxType{{0.0_x,0.5_x},{-1.0_x,1.0_x}});
    ARIADNE_TEST_EQUALS(restricted.number_of_constraints(),2u);
    auto affine=set.affine_approximation();
    auto enclosing_affine=set.affine_over_approximation();
    ARIADNE_TEST_EQUALS(affine.dimension(),2u);
    ARIADNE_TEST_EQUALS(enclosing_affine.dimension(),2u);
    set.reduce();
    ARIADNE_TEST_ASSERT(set.reduced_domain()[0].lower_bound()>=FloatDP(0,dp));
    ARIADNE_TEST_ASSERT(!definitely(set.is_empty()));
    ValidatedConstrainedImageSet impossible(domain,mapping,List<ValidatedConstraint>{x>=3});
    ARIADNE_TEST_ASSERT(definitely(impossible.is_empty()));
    ValidatedConstrainedImageSet box_set(domain);
    ARIADNE_TEST_ASSERT(definitely(box_set.inside(ExactBoxType{{-3.0_x,3.0_x},{-3.0_x,3.0_x}})));
    ARIADNE_TEST_ASSERT(definitely(box_set.separated(ExactBoxType{{3.0_x,4.0_x},{3.0_x,4.0_x}})));
    ARIADNE_TEST_ASSERT(definitely(box_set.overlaps(ExactBoxType{{-0.5_x,0.5_x},{-0.5_x,0.5_x}})));
    auto transformed=image(box_set,mapping);
    ARIADNE_TEST_EQUALS(transformed.dimension(),2u);
    ARIADNE_TEST_EQUALS(transformed.bounding_box()[0].lower_bound().raw(),FloatDP(-3,dp));
    ARIADNE_TEST_EQUALS(transformed.bounding_box()[0].upper_bound().raw(),FloatDP(3,dp));
    transformed.apply(ValidatedVectorMultivariateFunction::identity(2u));
    ARIADNE_TEST_EQUALS(transformed.dimension(),2u);
    auto wrong=ValidatedScalarMultivariateFunction::coordinate(1u,0u);
    ARIADNE_TEST_FAIL(ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x}},mapping));
    ARIADNE_TEST_FAIL(ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x}},mapping,List<ValidatedConstraint>{}));
    ARIADNE_TEST_FAIL(set.new_parameter_constraint(wrong<=ValidatedNumber(1)));
    ARIADNE_TEST_FAIL(set.new_space_constraint(wrong<=ValidatedNumber(1)));
    auto taylor=ValidatedVectorMultivariateTaylorFunctionModelDP(domain,mapping,ThresholdSweeper<FloatDP>(dp,1e-12));
    ValidatedConstrainedImageSet model_set(domain,ValidatedVectorMultivariateFunctionModelDP(taylor));
    ARIADNE_TEST_EQUALS(model_set.dimension(),2u);
    ARIADNE_TEST_FAIL(ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x}},ValidatedVectorMultivariateFunctionModelDP(taylor)));
    Grid grid(2u);
    auto outer=box_set.outer_approximation(grid,0u);
    ARIADNE_TEST_ASSERT(!outer.is_empty());
    GridTreePaving another_outer(grid);
    box_set.adjoin_outer_approximation_to(another_outer,0u);
    ARIADNE_TEST_ASSERT(!another_outer.is_empty());
    NullCanvas canvas;
    set.draw(canvas,Projection2d(2u,0u,1u));
    std::ostringstream repr;
    repr<<set;
    ARIADNE_TEST_ASSERT(repr.str().find("ValidatedConstrainedImageSet")!=String::npos);
}

Void test_one_sided_satisfaction() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    EffectiveVectorMultivariateFunction identity({x});
    ConstrainedImageSet effective(RealBox{{-2,2}},identity,List<EffectiveConstraint>{-0.5_x<=x<=0.5_x});
    ValidatedConstrainedImageSet validated(ExactBoxType{{-2.0_x,2.0_x}},identity,List<ValidatedConstraint>{-0.5_x<=x<=0.5_x});
    ARIADNE_TEST_ASSERT(definitely(effective.satisfies(x>=-1,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(validated.satisfies(x>=-1)));
    ARIADNE_TEST_ASSERT(definitely(effective.satisfies(x<=1,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(validated.satisfies(x<=1)));
    ARIADNE_TEST_ASSERT(definitely(effective.satisfies(-1<=x<=1,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(validated.satisfies(-1<=x<=1)));
    ARIADNE_TEST_ASSERT(definitely(!effective.satisfies(x>=1,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(!validated.satisfies(x>=1)));
    ARIADNE_TEST_ASSERT(definitely(!effective.satisfies(x<=-1,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(!validated.satisfies(x<=-1)));
    ARIADNE_TEST_ASSERT(definitely(effective.satisfies(x<=3,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(validated.satisfies(x<=3)));
    ConstrainedImageSet box(RealBox{{-2,2}});
    auto two_variables=EffectiveVectorMultivariateFunction::identity(2u);
    ARIADNE_TEST_FAIL(ConstrainedImageSet(RealBox{{0,1}},two_variables));
    ARIADNE_TEST_FAIL(ConstrainedImageSet(RealBox{{0,1}},two_variables,List<EffectiveConstraint>{}));
    ARIADNE_TEST_FAIL(box.new_parameter_constraint(EffectiveScalarMultivariateFunction::coordinate(2u,0u)<=1));
    ARIADNE_TEST_FAIL(box.new_space_constraint(EffectiveScalarMultivariateFunction::coordinate(2u,0u)<=1));
    std::unique_ptr<ConstrainedImageSet> clone(box.clone());
    ARIADNE_TEST_EQUALS(clone->dimension(),1u);
    auto split_box=box.split();
    ARIADNE_TEST_EQUALS(split_box.first.dimension(),1u);
    ARIADNE_TEST_ASSERT(definitely(box.inside(ExactBoxType{{-3.0_x,3.0_x}},Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(box.inside(ExactBoxType{{-3.0_x,3.0_x}}),Effort(1)));
    ARIADNE_TEST_EQUALS(image(box,identity).dimension(),1u);
}
Int main() {
    ARIADNE_TEST_CALL(test_default_image_sets());
    ARIADNE_TEST_CALL(test_validated_image_operations());
    ARIADNE_TEST_CALL(test_one_sided_satisfaction());
    return ARIADNE_TEST_FAILURES;
}
