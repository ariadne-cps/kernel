/***************************************************************************
 *            test_geometry.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "geometry/geometry.hpp"
#include "macro/testing.hpp"
using namespace Ariadne;
using Located=LocatedSetInterface<EffectiveTag,RealVector>;
using Regular=RegularSetInterface<EffectiveTag,RealVector>;

Void test_geometric_relations() {
    const FloatDP tolerance(0.015625_x,dp);
    ExactBoxSetType source(ExactBoxType{{0.0_x,1.0_x}});
    ExactBoxSetType enclosing(ExactBoxType{{-1.0_x,2.0_x}});
    ExactBoxSetType disjoint_set(ExactBoxType{{2.0_x,3.0_x}});
    ExactBoxSetType touching(ExactBoxType{{1.0_x,2.0_x}});
    const Located& located=source;
    const Regular& enclosing_regular=enclosing;
    const Regular& disjoint_regular=disjoint_set;
    const Regular& touching_regular=touching;
    ARIADNE_TEST_ASSERT(definitely(overlap(located,enclosing_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(inside(located,enclosing_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(!separated(located,enclosing_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(!overlap(located,disjoint_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(!inside(located,disjoint_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(separated(located,disjoint_regular,tolerance)));
    ARIADNE_TEST_ASSERT(is_indeterminate(overlap(located,touching_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(!inside(located,touching_regular,tolerance)));
    ARIADNE_TEST_ASSERT(is_indeterminate(separated(located,touching_regular,tolerance)));
    ExactBoxSetType empty(ExactBoxType(1u));
    const Located& empty_located=empty;
    ARIADNE_TEST_ASSERT(definitely(!overlap(empty_located,enclosing_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(inside(empty_located,enclosing_regular,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(separated(empty_located,enclosing_regular,tolerance)));
    ExactBoxType away{{4.0_x,5.0_x}};
    ARIADNE_TEST_ASSERT(definitely(!overlap(located,enclosing_regular,away,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(inside(located,enclosing_regular,away,tolerance)));
    ARIADNE_TEST_ASSERT(definitely(separated(located,enclosing_regular,away,tolerance)));

    ExactBoxType query{{0.0_x,1.0_x}};
    const ExactBoxType supports[]={ExactBoxType{{0.125_x,0.25_x}},ExactBoxType{{0.75_x,0.875_x}},ExactBoxType{{0.25_x,0.75_x}}};
    const ExactBoxType targets[]={ExactBoxType{{0.0625_x,0.3125_x}},ExactBoxType{{0.6875_x,0.9375_x}},ExactBoxType{{0.125_x,0.875_x}}};
    for(SizeType i=0u;i!=3u;++i) {
        ARIADNE_TEST_PRINT(supports[i]);
        ExactBoxSetType small_source(supports[i]);
        ExactBoxSetType small_target(targets[i]);
        const Located& small_located=small_source;
        const Regular& small_regular=small_target;
        ARIADNE_TEST_ASSERT(definitely(overlap(small_located,small_regular,query,tolerance)));
        ARIADNE_TEST_ASSERT(definitely(inside(small_located,small_regular,query,tolerance)));
        ARIADNE_TEST_ASSERT(definitely(!separated(small_located,small_regular,query,tolerance)));
    }
    // The query touches the source only on its boundary, so it supplies
    // no witness in the query's interior, even though the target covers it.
    ExactBoxSetType boundary_source(ExactBoxType{{0.0_x,0.25_x}});
    ExactBoxSetType boundary_target(ExactBoxType{{0.125_x,0.5_x}});
    const Located& boundary_located=boundary_source;
    const Regular& boundary_regular=boundary_target;
    ExactBoxType boundary_query{{0.25_x,0.375_x}};
    ARIADNE_TEST_ASSERT(is_indeterminate(overlap(boundary_located,boundary_regular,boundary_query,tolerance)));
    ARIADNE_TEST_ASSERT(is_indeterminate(separated(boundary_located,boundary_regular,boundary_query,tolerance)));
    ARIADNE_TEST_ASSERT(is_indeterminate(inside(boundary_located,disjoint_regular,boundary_query,tolerance)));
}
Int main() {
    ARIADNE_TEST_CALL(test_geometric_relations());
    return ARIADNE_TEST_FAILURES;
}
