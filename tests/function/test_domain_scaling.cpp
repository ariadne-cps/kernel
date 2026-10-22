/***************************************************************************
 *            test_domain_scaling.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "numeric/numeric.hpp"
#include "function/domain.hpp"
#include "function/scaling.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

Void test_domains() {
    RealDomain real,explicit_real{SizeOne()};
    ARIADNE_TEST_EQUALS(real.dimension(),1u);
    ARIADNE_TEST_ASSERT(real==explicit_real);
    ARIADNE_TEST_ASSERT(intersection(real,explicit_real)==real);
    IntervalDomainType entire=real;
    ARIADNE_TEST_EQUALS(entire.lower_bound(),FloatDP(-infty,dp));
    ARIADNE_TEST_EQUALS(entire.upper_bound(),FloatDP(+infty,dp));
    ARIADNE_TEST_ASSERT(real==entire);
    ARIADNE_TEST_ASSERT(entire==real);
    ARIADNE_TEST_ASSERT(!(real==IntervalDomainType(-1,1)));
    ARIADNE_TEST_ASSERT(!(IntervalDomainType(-1,1)==real));
    EuclideanDomain plane(2u),another_plane(2u,real),line(1u);
    ARIADNE_TEST_EQUALS(plane.dimension(),2u);
    ARIADNE_TEST_ASSERT(plane[0u]==real);
    ARIADNE_TEST_ASSERT(plane==another_plane);
    ARIADNE_TEST_ASSERT(!(plane==line));
    ARIADNE_TEST_ASSERT(intersection(plane,another_plane)==plane);
    ARIADNE_TEST_FAIL(intersection(plane,line));
    ARIADNE_TEST_EQUALS(product(plane,line).dimension(),3u);
    ARIADNE_TEST_EQUALS(product(plane,real).dimension(),3u);
    BoxDomainType entire_plane=plane;
    ARIADNE_TEST_EQUALS(entire_plane.dimension(),2u);
    ARIADNE_TEST_ASSERT(plane==entire_plane);
    ARIADNE_TEST_ASSERT(entire_plane==plane);
    ARIADNE_TEST_ASSERT(!(plane==BoxDomainType::unit_box(2u)));
    ARIADNE_TEST_ASSERT(!(BoxDomainType::unit_box(2u)==plane));
    UnitInterval unit;
    ARIADNE_TEST_EQUALS(unit.dimension(),1u);
    ARIADNE_TEST_ASSERT(unit==UnitInterval());
    ARIADNE_TEST_ASSERT(intersection(unit,unit)==unit);
    ARIADNE_TEST_EQUALS(IntervalDomainType(unit),IntervalDomainType(-1,1));
    // Both conversion APIs must describe the same signed unit interval.
    ARIADNE_TEST_EQUALS(unit.operator IntervalDomainType(),IntervalDomainType(unit));
    UnitBox square(2u),another_square(2u,unit),segment(1u);
    ARIADNE_TEST_EQUALS(square.dimension(),2u);
    ARIADNE_TEST_ASSERT(square[0u]==unit);
    ARIADNE_TEST_ASSERT(square==another_square);
    ARIADNE_TEST_ASSERT(!(square==segment));
    ARIADNE_TEST_ASSERT(intersection(square,another_square)==square);
    ARIADNE_TEST_FAIL(intersection(square,segment));
    ARIADNE_TEST_EQUALS(product(square,segment).dimension(),3u);
    ARIADNE_TEST_EQUALS(product(square,unit).dimension(),3u);
    ARIADNE_TEST_EQUALS(BoxDomainType(square),BoxDomainType::unit_box(2u));
    ARIADNE_TEST_EQUALS(BoxDomainType(UnitBox(0u)).dimension(),0u);
    ARIADNE_TEST_EQUALS(BoxDomainType(EuclideanDomain(0u)).dimension(),0u);
    std::ostringstream text;text<<real<<" "<<plane<<" "<<unit<<" "<<square;
    ARIADNE_TEST_EQUALS(text.str(),String("R R2 [-1:+1] [-1:+1]^2"));
}

Void test_scaling() {
    IntervalDomainType interval(2,6);
    ARIADNE_TEST_EQUALS(med(interval),Dyadic(4));
    ARIADNE_TEST_EQUALS(rad(interval),Dyadic(2));
    ARIADNE_TEST_EQUALS(med_apprx(interval),FloatDPApproximation(4,dp));
    ARIADNE_TEST_EQUALS(rad_apprx(interval),FloatDPApproximation(2,dp));
    ARIADNE_TEST_ASSERT(definitely(med_val(interval)==FloatDPBounds(4,dp)));
    ARIADNE_TEST_ASSERT(definitely(rad_val(interval)==FloatDPBounds(2,dp)));
    ARIADNE_TEST_EQUALS(convert_interval(interval,dp),interval);
    MultiplePrecision precision(128u);
    auto mp_interval=convert_interval(interval,precision);
    ARIADNE_TEST_EQUALS(mp_interval.lower_bound(),FloatMP(2,precision));
    ARIADNE_TEST_EQUALS(mp_interval.upper_bound(),FloatMP(6,precision));
    Scaling scaling(interval);Unscaling unscaling(interval);
    ARIADNE_TEST_ASSERT(scaling.domain()==UnitInterval());
    ARIADNE_TEST_EQUALS(scaling.codomain(),interval);
    ARIADNE_TEST_EQUALS(unscaling.domain(),interval);
    ARIADNE_TEST_ASSERT(unscaling.codomain()==UnitInterval());
    for(Int k=-1;k<=1;++k) {
        FloatDPApproximation x(k,dp);
        ARIADNE_TEST_EQUALS(scaling(x),FloatDPApproximation(4+2*k,dp));
        ARIADNE_TEST_EQUALS(unscaling(scaling(x)),x);
    }
    ARIADNE_TEST_EQUALS(unscale(FloatDPApproximation(3,dp),IntervalDomainType(3,3)),FloatDPApproximation(0,dp));
    BoxDomainType box{{2,6},{-3,-1}};
    VectorScaling vector_scaling(box);VectorUnscaling vector_unscaling(box);
    ARIADNE_TEST_EQUALS(vector_scaling.size(),2u);
    ARIADNE_TEST_EQUALS(vector_unscaling.size(),2u);
    ARIADNE_TEST_EQUALS(vector_scaling.codomain(),box);
    ARIADNE_TEST_EQUALS(vector_unscaling.domain(),box);
    ARIADNE_TEST_EQUALS(vector_scaling[0u].codomain(),box[0u]);
    ARIADNE_TEST_EQUALS(vector_unscaling[1u].domain(),box[1u]);
    Vector<FloatDPApproximation> point({-1,1},dp),expected({2,-1},dp);
    ARIADNE_TEST_EQUALS(vector_scaling(point),expected);
    ARIADNE_TEST_EQUALS(vector_unscaling(expected),point);
    Vector<FloatDPApproximation> wrong_dimension({0},dp);
    ARIADNE_TEST_FAIL(vector_scaling(wrong_dimension));
    ARIADNE_TEST_FAIL(vector_unscaling(wrong_dimension));
    Vector<FloatDPApproximation> empty_point(0u,dp);
    ARIADNE_TEST_EQUALS(VectorScaling(BoxDomainType())(empty_point).size(),0u);
    ARIADNE_TEST_EQUALS(VectorUnscaling(BoxDomainType())(empty_point).size(),0u);
    UnscalingException error("invalid domain ",interval);
    ARIADNE_TEST_EQUALS(error.domain,interval);
    ARIADNE_TEST_ASSERT(String(error.what()).find("invalid domain")!=String::npos);
}
Int main() {
    ARIADNE_TEST_CALL(test_domains());
    ARIADNE_TEST_CALL(test_scaling());
    return ARIADNE_TEST_FAILURES;
}
