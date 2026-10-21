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
#include "drawing/null_graphics.hpp"

namespace Ariadne { template<> double numeric_cast<double>(const CheckedAssignable<FloatMPApproximation>&); }
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
    // All three points must contribute to the range, including the final one.
    ARIADNE_TEST_ASSERT(std::abs(box[0].lower_bound().raw().get_d()+1.0)<1e-8);
    ARIADNE_TEST_ASSERT(std::abs(box[0].upper_bound().raw().get_d()-3.0)<1e-8);
    ARIADNE_TEST_ASSERT(std::abs(box[1].lower_bound().raw().get_d()-2.0)<1e-8);
    ARIADNE_TEST_ASSERT(std::abs(box[1].upper_bound().raw().get_d()-5.0)<1e-8);

    std::unique_ptr<InterpolatedCurve> cloned(curve.clone());
    ARIADNE_TEST_EQUALS(cloned->size(),3u);
    ARIADNE_TEST_EQUALS(cloned->dimension(),2u);

    std::ostringstream repr;
    repr<<curve;
    ARIADNE_TEST_ASSERT(repr.str().find("InterpolatedCurve( size=3")!=std::string::npos);
}

class CurveCanvas : public NullCanvas {
  public:
    Nat moves=0, lines=0, fills=0, strokes=0;
    double last_x=0, last_y=0;
    Void move_to(double x, double y) override { ++moves; last_x=x; last_y=y; }
    Void line_to(double x, double y) override { ++lines; last_x=x; last_y=y; }
    Void fill() override { ++fills; }
    Void stroke() override { ++strokes; }
};

Void test_curve_overloads_and_drawing() {
    using P=InterpolatedCurve::PointType;
    P a({FloatDPApproximation(1,dp),FloatDPApproximation(2,dp)});
    P b({FloatDPApproximation(3,dp),FloatDPApproximation(4,dp)});
    InterpolatedCurve empty;
    ARIADNE_TEST_EQUALS(empty.size(),0u);
    empty.insert(Dyadic(0),a);
    ARIADNE_TEST_EQUALS(empty.size(),1u);
    InterpolatedCurve parameter(FloatDP(2,dp),a);
    ARIADNE_TEST_EQUALS(parameter.begin()->first,FloatDP(2,dp));
    InterpolatedCurve dyadic(Dyadic(3),a);
    ARIADNE_TEST_EQUALS(dyadic.begin()->first,FloatDP(3,dp));
    Vector<FloatDP> raw({1.0_x,2.0_x},dp);
    InterpolatedCurve from_raw(FloatDP(0,dp),raw);
    from_raw.insert(FloatDP(1,dp),raw);
    ARIADNE_TEST_EQUALS(from_raw.size(),2u);
    P wrong_dimension({FloatDPApproximation(1,dp)});
    ARIADNE_TEST_FAIL(from_raw.insert(FloatDP(2,dp),wrong_dimension));
    ARIADNE_TEST_FAIL(from_raw.insert(FloatDP(2,dp),Vector<FloatDP>({1.0_x},dp)));
    ARIADNE_TEST_ASSERT(decide(from_raw.begin()->second==a));

    CurveCanvas canvas;
    InterpolatedCurve open(a,b);
    open.draw(canvas,Projection2d(2,0,1));
    ARIADNE_TEST_EQUALS(canvas.lines,1u);
    ARIADNE_TEST_EQUALS(canvas.strokes,1u);
    ARIADNE_TEST_EQUALS(canvas.last_x,3.0);
    ARIADNE_TEST_EQUALS(canvas.last_y,4.0);
    open.insert(Dyadic(2),a);
    open.draw(canvas,Projection2d(2,0,1));
    ARIADNE_TEST_EQUALS(canvas.fills,1u);
    ARIADNE_TEST_EQUALS(canvas.lines,3u);
    from_raw.draw(canvas,Projection2d(2,0,1));
    ARIADNE_TEST_EQUALS(canvas.fills,2u);
    std::ostringstream repr;
    open._write(repr);
    ARIADNE_TEST_ASSERT(repr.str().find("InterpolatedCurve")!=String::npos);
    auto x=EffectiveScalarUnivariateFunction::coordinate();
    Curve function_curve(EffectiveVectorUnivariateFunction({x,x*x}));
    const CurveInterface& interface=function_curve;
    repr<<interface;
    ARIADNE_TEST_ASSERT(repr.str().find("Curve( function=")!=String::npos);
    Vector<FloatMPApproximation> mp({1.5_x},precision(128));
    ARIADNE_TEST_EQUALS(numeric_cast<double>(mp[0]),1.5);
}

Int main() {
    ARIADNE_TEST_CALL(test_curve_overloads_and_drawing());
    ARIADNE_TEST_CALL(test_function_curve());
    ARIADNE_TEST_CALL(test_interpolated_curve());
    return ARIADNE_TEST_FAILURES;
}
