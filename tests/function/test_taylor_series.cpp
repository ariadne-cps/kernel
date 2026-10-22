/***************************************************************************
 *            test_taylor_series.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "numeric/numeric.hpp"
#include "algebra/series.hpp"
#include "function/taylor_series.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

Void test_polynomial_series() {
    IntervalDomainType domain(-1,1);
    TaylorSeries<FloatDPBounds> empty_series(domain,3u);
    ARIADNE_TEST_EQUALS(empty_series.degree(),3u);
    ARIADNE_TEST_EQUALS(empty_series.expansion().size(),4u);
    ARIADNE_TEST_EQUALS(empty_series.error().raw(),FloatDP(0,dp));
    TaylorSeries<FloatDPBounds> square(Sqr(),domain,FloatDP(0,dp),2u);
    ARIADNE_TEST_EQUALS(square.degree(),2u);
    ARIADNE_TEST_EQUALS(square[0u],FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(square[1u],FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(square[2u],FloatDP(1,dp));
    ARIADNE_TEST_EQUALS(square.error().raw(),FloatDP(0,dp));
    square.sweep(FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(square[2u],FloatDP(1,dp));
    ARIADNE_TEST_EQUALS(square.error().raw(),FloatDP(0,dp));
    square.sweep(FloatDP(1,dp));
    ARIADNE_TEST_EQUALS(square[2u],FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(square.error().raw(),FloatDP(1,dp));
    TaylorSeries<FloatDPBounds> constant(Sqr(),IntervalDomainType(2,2),FloatDP(2,dp),0u);
    ARIADNE_TEST_EQUALS(constant[0u],FloatDP(4,dp));
    ARIADNE_TEST_EQUALS(constant.error().raw(),FloatDP(0,dp));
    std::ostringstream text;text<<square;
    ARIADNE_TEST_ASSERT(text.str().find("TS(")==0u);
}
Int main() {
    ARIADNE_TEST_CALL(test_polynomial_series());
    return ARIADNE_TEST_FAILURES;
}
