/***************************************************************************
 *            test_set_adapters.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "geometry/list_set.hpp"
#include "geometry/set_wrapper.hpp"
#include "geometry/union_of_intervals.hpp"
#include "function/taylor_multifunction.hpp"
#include "drawing/null_graphics.hpp"
#include "macro/testing.hpp"
#include <memory>
#include <sstream>
using namespace Ariadne;

Void test_list_set() {
    ExactBoxType a{{0.0_x,1.0_x},{0.0_x,1.0_x}}, b{{2.0_x,3.0_x},{2.0_x,3.0_x}};
    ListSet<ExactBoxType> sets(2u);
    ARIADNE_TEST_ASSERT(sets.empty());
    ARIADNE_TEST_EQUALS(sets.dimension(),0u);
    ARIADNE_TEST_EQUALS(sets.bounding_box().dimension(),0u);
    ARIADNE_TEST_ASSERT(sets.bounding_boxes().empty());
    std::ostringstream empty_text; empty_text<<sets;
    ARIADNE_TEST_EQUALS(empty_text.str(),String("ListSet( )"));
    sets.adjoin(a); sets.adjoin(b);
    ARIADNE_TEST_EQUALS(sets.size(),2u);
    ARIADNE_TEST_ASSERT(sets.capacity()>=sets.size());
    ARIADNE_TEST_EQUALS(sets.dimension(),2u);
    ARIADNE_TEST_EQUALS(sets[1],b);
    ARIADNE_TEST_ASSERT(same(cast_exact_box(sets.bounding_box()),ExactBoxType{{0.0_x,3.0_x},{0.0_x,3.0_x}}));
    ARIADNE_TEST_EQUALS(sets.bounding_boxes().size(),2u);
    std::unique_ptr<ListSet<ExactBoxType>> copy(sets.clone());
    ARIADNE_TEST_EQUALS(*copy,sets);
    auto summary=sets.summary();
    ARIADNE_TEST_EQUALS(summary.size,2u);
    ARIADNE_TEST_EQUALS(summary.dimension,2u);
    std::ostringstream text; sets._write(text); text<<summary;
    ARIADNE_TEST_ASSERT(text.str().find("s=2, d=2")!=String::npos);
    const auto& const_sets=sets;
    ListSet<ExactBoxType> from_iterators(const_sets.begin(),const_sets.end());
    ARIADNE_TEST_EQUALS(from_iterators,sets);
    ListSet<ExactBoxType> from_vector(std::vector<ExactBoxType>{a,b});
    ARIADNE_TEST_EQUALS(from_vector,sets);
    ListSet<UpperBoxType> converted(sets);
    ARIADNE_TEST_EQUALS(converted.size(),2u);
    copy->erase(copy->begin());
    ARIADNE_TEST_EQUALS(copy->pop(),b);
    ARIADNE_TEST_ASSERT(copy->empty());
    copy->push_back(a); copy->clear();
    copy->adjoin(sets);
    copy->adjoin(ListSet<ExactBoxType>());
    copy->adjoin(ListSet<ExactBoxType>(a));
    ARIADNE_TEST_EQUALS(copy->size(),3u);
    ARIADNE_TEST_FAIL(sets.adjoin(ExactBoxType{{0.0_x,1.0_x}}));
    ARIADNE_TEST_FAIL(sets.adjoin(ListSet<ExactBoxType>(ExactBoxType{{0.0_x,1.0_x}})));
    NullCanvas canvas;
    sets.draw(canvas,Projection2d(2u,0u,1u));
}

template<class C> Void check_set_adapter() {
    ExactBoxType shape{{0.0_x,1.0_x}}, inner{{0.25_x,0.75_x}}, outside{{2.0_x,3.0_x}}, enclosing{{-1.0_x,2.0_x}};
    using W=WrapperType<ExactBoxType,C>;
    W wrapper(shape);
    C handle(shape);
    C cloned=wrap<C>(wrapper);
    ARIADNE_TEST_EQUALS(handle.dimension(),1u);
    ARIADNE_TEST_EQUALS(cloned.dimension(),1u);
    if constexpr(requires { handle.overlaps(inner); }) {
        ARIADNE_TEST_ASSERT(definitely(handle.overlaps(inner)));
        ARIADNE_TEST_ASSERT(!definitely(handle.overlaps(outside)));
    }
    if constexpr(requires { handle.covers(inner); }) {
        ARIADNE_TEST_ASSERT(definitely(handle.covers(inner)));
        ARIADNE_TEST_ASSERT(!definitely(handle.covers(outside)));
    }
    if constexpr(requires { handle.separated(outside); }) {
        ARIADNE_TEST_ASSERT(definitely(handle.separated(outside)));
        ARIADNE_TEST_ASSERT(!definitely(handle.separated(inner)));
    }
    if constexpr(requires { handle.inside(enclosing); }) {
        ARIADNE_TEST_ASSERT(definitely(handle.inside(enclosing)));
        ARIADNE_TEST_ASSERT(!definitely(handle.inside(outside)));
        ARIADNE_TEST_ASSERT(same(cast_exact_box(handle.bounding_box()),shape));
    }
    std::ostringstream expected,actual,direct;
    expected<<shape; actual<<handle; direct<<wrapper;
    ARIADNE_TEST_EQUALS(actual.str(),expected.str());
    ARIADNE_TEST_EQUALS(direct.str(),expected.str());
}

Void test_set_adapters() {
    ARIADNE_TEST_CALL((check_set_adapter<OpenSet<ValidatedTag,RealVector>>()));
    ARIADNE_TEST_CALL((check_set_adapter<ClosedSet<ValidatedTag,RealVector>>()));
    ARIADNE_TEST_CALL((check_set_adapter<OvertSet<ValidatedTag,RealVector>>()));
    ARIADNE_TEST_CALL((check_set_adapter<CompactSet<ValidatedTag,RealVector>>()));
    ARIADNE_TEST_CALL((check_set_adapter<RegularSet<ValidatedTag,RealVector>>()));
    ARIADNE_TEST_CALL((check_set_adapter<LocatedSet<ValidatedTag,RealVector>>()));
    ARIADNE_TEST_CALL((check_set_adapter<RegularLocatedSet<ValidatedTag,RealVector>>()));
    EuclideanSpace space(3u);
    ARIADNE_TEST_EQUALS(space.dimension(),3u);
    UnionOfIntervals<FloatDP> intervals({{0,1}},dp);
    RegularSet<ValidatedTag,Real> regular(intervals);
    ARIADNE_TEST_ASSERT(definitely(regular.covers(ExactIntervalType(0.25_x,0.75_x))));
    ARIADNE_TEST_ASSERT(definitely(regular.separated(ExactIntervalType(2,3))));
    auto open=wrap_open(intervals);
    ARIADNE_TEST_EQUALS(open.dimension(),1u);
    ARIADNE_TEST_ASSERT(definitely(open.overlaps(ExactIntervalType(0.25_x,0.75_x))));
    ARIADNE_TEST_ASSERT(definitely(open.covers(ExactIntervalType(0.25_x,0.75_x))));
    std::unique_ptr<SetInterfaceBase<Real>> cloned(open.reference().clone());
    ARIADNE_TEST_EQUALS(cloned->dimension(),1u);
    IntervalSet<FloatDPUpperBound> interval(UpperIntervalType(0,1));
    CompactSet<ValidatedTag,Real> compact(interval);
    ARIADNE_TEST_EQUALS(compact.dimension(),1u);
    ARIADNE_TEST_ASSERT(definitely(compact.inside(ExactIntervalType(-1,2))));
    ARIADNE_TEST_ASSERT(definitely(compact.separated(ExactIntervalType(2,3))));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(compact.bounding_box()),ExactIntervalType(0,1)));
    std::unique_ptr<SetInterfaceBase<Real>> compact_clone(compact.reference().clone());
    ARIADNE_TEST_EQUALS(compact_clone->dimension(),1u);
    std::ostringstream text; text<<compact;
    ARIADNE_TEST_ASSERT(!text.str().empty());
}
Int main() {
    ARIADNE_TEST_CALL(test_list_set());
    ARIADNE_TEST_CALL(test_set_adapters());
    return ARIADNE_TEST_FAILURES;
}
