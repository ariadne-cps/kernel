/***************************************************************************
 *            test_union_of_intervals.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "geometry/union_of_intervals.hpp"
#include "geometry/set.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

Void test_interval_union() {
    using U=UnionOfIntervals<Dyadic>;
    U empty({});
    ARIADNE_TEST_EQUALS(empty.size(),0u);
    ARIADNE_TEST_EQUALS(empty.measure(),0);
    ARIADNE_TEST_FAIL(empty.bounding_box());
    U left({{0,1},{3,4}});
    U right({{2,3.5_dy}});
    ARIADNE_TEST_EQUALS(intersection(left,right).measure(),0.5_dy);
    ARIADNE_TEST_EQUALS(intersection(right,left).measure(),0.5_dy);
    ARIADNE_TEST_EQUALS(intersection(empty,left).size(),0u);
    ARIADNE_TEST_EQUALS(intersection(left,empty).size(),0u);
    ARIADNE_TEST_ASSERT(!intersect(empty,left));
    ARIADNE_TEST_ASSERT(!intersect(left,empty));
    ARIADNE_TEST_ASSERT(intersect(left,right));
    ARIADNE_TEST_ASSERT(intersect(right,left));
    ARIADNE_TEST_ASSERT(!intersect(U({{0,1}}),U({{2,3}})));
    ARIADNE_TEST_ASSERT(!intersect(U({{2,3}}),U({{0,1}})));
    ARIADNE_TEST_ASSERT(subset(empty,left));
    ARIADNE_TEST_ASSERT(!subset(left,empty));
    ARIADNE_TEST_ASSERT(subset(left,left));
    ARIADNE_TEST_ASSERT(subset(U({{3,3.5_dy}}),left));
    ARIADNE_TEST_ASSERT(!subset(U({{-1,0.5_dy}}),left));
    ARIADNE_TEST_ASSERT(!subset(U({{0.5_dy,1.5_dy}}),left));
    ARIADNE_TEST_EQUALS(left.bounding_box(),DyadicInterval(0,4));
    ARIADNE_TEST_EQUALS(measure(left),2);
    ARIADNE_TEST_EQUALS(left.dimension(),SizeOne());
    ARIADNE_TEST_ASSERT(left.separated(IntervalDomainType(1.5_dy,2.5_dy)));
    ARIADNE_TEST_ASSERT(left.overlaps(IntervalDomainType(0.5_dy,2.5_dy)));
    ARIADNE_TEST_ASSERT(left.covers(IntervalDomainType(0.25_dy,0.75_dy)));
    ARIADNE_TEST_ASSERT(left.inside(IntervalDomainType(-1,5)));
    ARIADNE_TEST_EQUALS(intersection(left,IntervalDomainType(0.5_dy,3.5_dy)).measure(),1);
    U sorted({{3,4},{0,2},{0,1},{0,3},{5,6}});
    ARIADNE_TEST_EQUALS(sorted.size(),2u);
    ARIADNE_TEST_EQUALS(sorted.measure(),5);
    ARIADNE_TEST_EQUALS(sorted.begin()->lower_bound(),0);
    const U& const_sorted=sorted;
    ARIADNE_TEST_EQUALS(const_sorted.begin()->upper_bound(),4);
    ARIADNE_TEST_EQUALS(interior(left,Dyadic(0.25_dy)).measure(),1);
    ARIADNE_TEST_EQUALS(interior(left,Dyadic(1)).size(),0u);
    LowerUpperBoundLess less;
    ARIADNE_TEST_ASSERT(less(DyadicInterval(0,1),DyadicInterval(1,2)));
    ARIADNE_TEST_ASSERT(less(DyadicInterval(0,1),DyadicInterval(0,2)));
    ARIADNE_TEST_ASSERT(!less(DyadicInterval(0,2),DyadicInterval(0,1)));
    ARIADNE_TEST_ASSERT(!less(DyadicInterval(1,2),DyadicInterval(0,2)));
    std::ostringstream repr;
    repr<<left;
    ARIADNE_TEST_ASSERT(!repr.str().empty());
}

Void test_floating_interval_union() {
    UnionOfIntervals<FloatDP> floating({{3,4},{0,2},{1,3}},dp);
    ARIADNE_TEST_EQUALS(floating.size(),1u);
    ARIADNE_TEST_EQUALS(floating.measure().upper().raw(),FloatDP(4,dp));
    ARIADNE_TEST_ASSERT(floating.overlaps(IntervalDomainType(1,2)));
    ARIADNE_TEST_EQUALS(floating.create(IntervalDomainType(1,2)).measure().upper().raw(),FloatDP(1,dp));
    UnionOfIntervals<FloatDP> empty({},dp);
    ARIADNE_TEST_EQUALS(empty.measure().upper().raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(empty.create(IntervalDomainType(1,2)).size(),1u);
    UnionOfIntervals<FloatMP> multiprecision({{0,2}},precision(128));
    ARIADNE_TEST_EQUALS(multiprecision.create(IntervalDomainType(1,2)).size(),1u);
    UnionOfIntervals<FloatMP> empty_mp({},precision(128));
    ARIADNE_TEST_EQUALS(empty_mp.create(IntervalDomainType(1,2)).size(),1u);
    ARIADNE_TEST_FAIL(empty.bounding_box());
    LowerUpperBoundLess less;
    ARIADNE_TEST_ASSERT(less(ExactIntervalType(0,1),ExactIntervalType(1,2)));
    ARIADNE_TEST_ASSERT(less(ExactIntervalType(0,1),ExactIntervalType(0,2)));
    ARIADNE_TEST_ASSERT(!less(ExactIntervalType(0,2),ExactIntervalType(0,1)));
    ARIADNE_TEST_ASSERT(!less(ExactIntervalType(1,2),ExactIntervalType(0,2)));
    // Raw floating comparisons use short-circuit Boolean conditions;
    // exercise every outcome in this instantiation as well.
    auto make=[](InitializerList<DyadicInterval> intervals) { return UnionOfIntervals<FloatDP>(intervals,dp); };
    auto first=make({{0,1},{3,4}});
    auto second=make({{2,3.5_dy}});
    auto no_intervals=make({});
    ARIADNE_TEST_EQUALS(first.dimension(),SizeOne());
    ARIADNE_TEST_EQUALS(first.bounding_box(),ExactIntervalType(0,4));
    ARIADNE_TEST_ASSERT(first.inside(IntervalDomainType(-1,5)));
    ARIADNE_TEST_ASSERT(first.covers(IntervalDomainType(0.25_dy,0.75_dy)));
    ARIADNE_TEST_ASSERT(first.separated(IntervalDomainType(1.5_dy,2.5_dy)));
    ARIADNE_TEST_ASSERT(intersect(first,second));
    ARIADNE_TEST_ASSERT(intersect(second,first));
    ARIADNE_TEST_ASSERT(!intersect(first,no_intervals));
    ARIADNE_TEST_ASSERT(!intersect(no_intervals,first));
    ARIADNE_TEST_ASSERT(!intersect(make({{0,1}}),make({{2,3}})));
    ARIADNE_TEST_ASSERT(!intersect(make({{2,3}}),make({{0,1}})));
    ARIADNE_TEST_ASSERT(subset(no_intervals,first));
    ARIADNE_TEST_ASSERT(!subset(first,no_intervals));
    ARIADNE_TEST_ASSERT(subset(first,first));
    ARIADNE_TEST_ASSERT(subset(make({{3,3.5_dy}}),first));
    ARIADNE_TEST_ASSERT(!subset(make({{-1,0.5_dy}}),first));
    ARIADNE_TEST_ASSERT(!subset(make({{0.5_dy,1.5_dy}}),first));
    auto survivor=make({{0,1},{0.5_dy,2},{3,4}});
    ARIADNE_TEST_EQUALS(survivor.size(),2u);
    ARIADNE_TEST_EQUALS(survivor.measure().upper().raw(),FloatDP(3,dp));
    std::ostringstream floating_repr;
    floating_repr<<first;
    ARIADNE_TEST_ASSERT(!floating_repr.str().empty());
    UnionOfIntervals<Dyadic> exact({{0,1}});
    UnionOfIntervals<Real> converted(exact);
    ARIADNE_TEST_EQUALS(converted.size(),1u);
}
Int main() {
    ARIADNE_TEST_CALL(test_interval_union());
    ARIADNE_TEST_CALL(test_floating_interval_union());
    return ARIADNE_TEST_FAILURES;
}
