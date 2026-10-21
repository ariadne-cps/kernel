/***************************************************************************
 *            test_constraint_set.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "geometry/function_set.hpp"
#include "geometry/affine_set.hpp"
#include "drawing/null_graphics.hpp"
#include "macro/testing.hpp"
#include <memory>
#include <sstream>
using namespace Ariadne;

Void check_predicate_clone(const LowerKleenean& predicate) {
    ARIADNE_TEST_ASSERT(definitely(predicate,Effort(1)));
    std::unique_ptr<Detail::LogicalInterface> copied(predicate.repr().reference()._copy());
    ARIADNE_TEST_ASSERT(Detail::definitely(copied->_check(Effort(1))));
    std::ostringstream repr;
    predicate.repr().reference()._write(repr);
    copied->_write(repr);
    ARIADNE_TEST_ASSERT(!repr.str().empty());
}

Void test_constraint_sets() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    EffectiveVectorMultivariateFunction identity({x,y});
    ConstraintSet constraint(identity,RealBox{{-1,1},{-1,1}});
    ConstraintSet list_constraint(List<EffectiveConstraint>{-1<=x<=1});
    ARIADNE_TEST_FAIL(ConstraintSet(identity,RealBox{{-1,1}}));
    ARIADNE_TEST_FAIL(ConstraintSet(List<EffectiveConstraint>{x<=1,EffectiveScalarMultivariateFunction::coordinate(1u,0u)<=1}));
    ARIADNE_TEST_FAIL(BoundedConstraintSet(RealBox{{0,1}},identity,RealBox{{-1,1},{-1,1}}));
    ConstraintSet empty_constraint(List<EffectiveConstraint>{});
    ARIADNE_TEST_EQUALS(empty_constraint.dimension(),0u);
    ARIADNE_TEST_EQUALS(constraint.dimension(),2u);
    ARIADNE_TEST_EQUALS(constraint.number_of_constraints(),2u);
    ARIADNE_TEST_EQUALS(constraint.constraints().size(),2u);
    ARIADNE_TEST_EQUALS(constraint.constraint(0u).function().argument_size(),2u);
    ARIADNE_TEST_EQUALS(constraint.codomain().dimension(),2u);
    ARIADNE_TEST_EQUALS(constraint.constraint_function().result_size(),2u);
    std::unique_ptr<ConstraintSet> clone(constraint.clone());
    ARIADNE_TEST_EQUALS(clone->dimension(),2u);
    ExactBoxType centre{{-0.25_x,0.25_x},{-0.25_x,0.25_x}};
    ExactBoxType far{{3.0_x,4.0_x},{3.0_x,4.0_x}};
    check_predicate_clone(constraint.covers(centre));
    check_predicate_clone(constraint.overlaps(centre));
    check_predicate_clone(constraint.separated(far));
    ARIADNE_TEST_ASSERT(definitely(constraint.covers(centre,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(constraint.overlaps(centre,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(constraint.separated(far,Effort(1))));
    ARIADNE_TEST_EQUALS(intersection(constraint,list_constraint).number_of_constraints(),3u);
    RealBox domain{{-2,2},{-3,3}};
    BoundedConstraintSet box_only(domain);
    BoundedConstraintSet bounded(domain,identity,RealBox{{-1,1},{-1,1}});
    BoundedConstraintSet bounded_list(domain,List<EffectiveConstraint>{-1<=x<=1});
    ARIADNE_TEST_EQUALS(bounded.number_of_constraints(),2u);
    ARIADNE_TEST_EQUALS(bounded.constraints().size(),2u);
    ARIADNE_TEST_EQUALS(bounded.constraint(1u).function().argument_size(),2u);
    ARIADNE_TEST_EQUALS(bounded.codomain().dimension(),2u);
    ARIADNE_TEST_EQUALS(bounded.constraint_function().result_size(),2u);
    std::unique_ptr<BoundedConstraintSet> bounded_clone(bounded.clone());
    ARIADNE_TEST_EQUALS(bounded_clone->dimension(),2u);
    check_predicate_clone(bounded.covers(centre));
    check_predicate_clone(bounded.overlaps(centre));
    check_predicate_clone(bounded.separated(far));
    ExactBoxType large{{-4.0_x,4.0_x},{-4.0_x,4.0_x}};
    check_predicate_clone(bounded.inside(large));
    ARIADNE_TEST_ASSERT(!definitely(bounded.covers(far,Effort(1))));
    ARIADNE_TEST_ASSERT(!definitely(bounded.overlaps(far,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(bounded.separated(ExactBoxType{{1.5_x,1.75_x},{1.5_x,1.75_x}},Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(box_only.overlaps(centre,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(box_only.inside(large,Effort(1))));
    ARIADNE_TEST_ASSERT(definitely(subset(ExactBoxType{{-2.0_x,2.0_x},{-3.0_x,3.0_x}},cast_exact_box(bounded.bounding_box()))));
    RealBox clipped_domain{{-0.5_x,0.5_x},{-2,2}};
    auto clipped=intersection(bounded_list,clipped_domain);
    auto reverse_clipped=intersection(clipped_domain,bounded_list);
    ARIADNE_TEST_EQUALS(clipped.dimension(),2u);
    ARIADNE_TEST_ASSERT(decide(clipped.domain()==clipped_domain));
    ARIADNE_TEST_ASSERT(decide(reverse_clipped.domain()==clipped_domain));
    auto combined=intersection(bounded_list,BoundedConstraintSet(clipped_domain));
    ARIADNE_TEST_ASSERT(decide(combined.domain()==clipped_domain));
    ARIADNE_TEST_EQUALS(combined.number_of_constraints(),1u);
    auto combined_constraints=intersection(bounded_list,constraint);
    auto reverse_combined_constraints=intersection(constraint,bounded_list);
    ARIADNE_TEST_ASSERT(decide(combined_constraints.domain()==domain));
    ARIADNE_TEST_ASSERT(decide(reverse_combined_constraints.domain()==domain));
    ARIADNE_TEST_EQUALS(combined_constraints.number_of_constraints(),3u);
    ARIADNE_TEST_EQUALS(intersection(constraint,domain).dimension(),2u);
    ARIADNE_TEST_EQUALS(intersection(domain,constraint).dimension(),2u);
    auto mapped=image(bounded_list,identity);
    ARIADNE_TEST_EQUALS(mapped.number_of_constraints(),1u);
    ARIADNE_TEST_EQUALS(mapped.dimension(),2u);
    ConstrainedImageSet converted(bounded_list);
    ARIADNE_TEST_EQUALS(converted.number_of_constraints(),1u);
    ARIADNE_TEST_EQUALS(class_name<ConstraintSet>(),String("EffectiveConstraintSet"));
    ARIADNE_TEST_EQUALS(class_name<BoundedConstraintSet>(),String("EffectiveBoundedConstraintSet"));
    ARIADNE_TEST_EQUALS(class_name<ConstrainedImageSet>(),String("EffectiveConstrainedImageSet"));
    ARIADNE_TEST_EQUALS(class_name<ValidatedConstrainedImageSet>(),String("ValidatedConstrainedImageSet"));
    ARIADNE_TEST_EQUALS(class_name<ValidatedAffineConstrainedImageSet>(),String("ValidatedAffineConstrainedImageSet"));
    NullCanvas canvas;
    bounded.draw(canvas,Projection2d(2u,0u,1u));
    std::ostringstream repr;
    repr<<constraint<<bounded<<bounded_list;
    ARIADNE_TEST_ASSERT(repr.str().find("BoundedConstraintSet")!=String::npos);
}

Int main() {
    ARIADNE_TEST_CALL(test_constraint_sets());
    return ARIADNE_TEST_FAILURES;
}
