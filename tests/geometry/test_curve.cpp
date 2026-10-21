/***************************************************************************
 *            test_curve.cpp
 *
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 ****************************************************************************/

#include <cmath>
#include <memory>
#include <sstream>

#include "geometry/curve.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

Void test_function_curve() {
    auto x=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction f({x,x*x});
    Curve curve(f);
    ARIADNE_TEST_EQUALS(curve.dimension(),2u);
    ARIADNE_TEST_ASSERT(curve.smoothness()==SMOOTH);

    Curve::ParameterType s(2.0_x,dp);
    auto pt=curve.value(s);
    ARIADNE_TEST_EQUALS(pt.size(),2u);
    ARIADNE_TEST_ASSERT(std::abs(pt[0].raw().get_d()-2.0)<1e-8);
    ARIADNE_TEST_ASSERT(std::abs(pt[1].raw().get_d()-4.0)<1e-8);

    auto tangent=curve.tangent(s);
    ARIADNE_TEST_EQUALS(tangent.size(),2u);
    ARIADNE_TEST_ASSERT(std::abs(tangent[0].raw().get_d()-1.0)<1e-8);
    ARIADNE_TEST_ASSERT(std::abs(tangent[1].raw().get_d()-4.0)<1e-8);

    Curve copy(curve);
    ARIADNE_TEST_EQUALS(copy.dimension(),2u);
    ARIADNE_TEST_EQUALS(copy.value(s).size(),2u);
    std::unique_ptr<Curve> cloned(curve.clone());
    ARIADNE_TEST_EQUALS(cloned->dimension(),2u);
    ARIADNE_TEST_ASSERT(std::abs(cloned->value(s)[1].raw().get_d()-4.0)<1e-8);

    std::ostringstream repr;
    repr<<curve;
    ARIADNE_TEST_ASSERT(repr.str().find("Curve( function=")!=std::string::npos);
}

Void test_interpolated_curve() {
    using P=InterpolatedCurve::PointType;
    P a({FloatDPApproximation(1,dp),FloatDPApproximation(2,dp)});
    P b({FloatDPApproximation(3,dp),FloatDPApproximation(4,dp)});
    P c({FloatDPApproximation(-1,dp),FloatDPApproximation(5,dp)});

    InterpolatedCurve curve(a);
    ARIADNE_TEST_EQUALS(curve.size(),1u);
    ARIADNE_TEST_EQUALS(curve.dimension(),2u);
    curve.insert(InterpolatedCurve::ParameterType(1.0_x,dp),b);
    curve.insert(Dyadic(2),c);
    ARIADNE_TEST_EQUALS(curve.size(),3u);

    auto box=curve.bounding_box();
    ARIADNE_TEST_EQUALS(box.dimension(),2u);
    ARIADNE_TEST_PRINT(box);
    ARIADNE_TEST_ASSERT(definitely(subset(box,ExactBoxType({{-2.0_x,4.0_x},{1.0_x,6.0_x}}))));

    std::unique_ptr<InterpolatedCurve> cloned(curve.clone());
    ARIADNE_TEST_EQUALS(cloned->size(),3u);
    ARIADNE_TEST_EQUALS(cloned->dimension(),2u);

    std::ostringstream repr;
    repr<<curve;
    ARIADNE_TEST_ASSERT(repr.str().find("InterpolatedCurve( size=3")!=std::string::npos);
}

Int main() {
    ARIADNE_TEST_CALL(test_function_curve());
    ARIADNE_TEST_CALL(test_interpolated_curve());
    return ARIADNE_TEST_FAILURES;
}
