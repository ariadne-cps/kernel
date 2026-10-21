/***************************************************************************
 *            test_point.cpp
 *
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 ****************************************************************************/

#include "geometry/point.hpp"
#include "geometry/point.tpl.hpp"
#include "drawing/null_graphics.hpp"
#include <cmath>
#include <memory>
#include <sstream>
#include "utility/string.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

// These explicit specializations provide the stable human-readable names
// used by diagnostics and other type-aware clients.
Void test_point_class_names() {
    ARIADNE_TEST_EQUALS(class_name<Point<Real>>(), String("Point<Real>"));
    ARIADNE_TEST_EQUALS(class_name<Point<FloatDP>>(), String("Point<FloatDP>"));
    ARIADNE_TEST_EQUALS(class_name<Point<FloatDPBounds>>(), String("Point<FloatDPBounds>"));
    ARIADNE_TEST_EQUALS(class_name<Point<FloatDPApproximation>>(), String("Point<FloatDPApproximation>"));
}

class PointCanvas : public NullCanvas {
  public:
    double x=0,y=0;
    Nat dots=0,strokes=0;
    Void dot(double px,double py) override { x=px; y=py; ++dots; }
    Void stroke() override { ++strokes; }
};

Void test_point_geometry() {
    Point<FloatDPApproximation> pt({1.0_x,2.0_x},dp);
    Point<FloatDPApproximation> zero;
    ARIADNE_TEST_EQUALS(zero.dimension(),0u);
    Point<FloatDP> generated(2u,[](SizeType i){return FloatDP(static_cast<Int>(i),dp);});
    Point<Real> generated_real(2u,[](SizeType i){return Real(static_cast<Int>(i));});
    ARIADNE_TEST_EQUALS(generated[1],FloatDP(1,dp));
    ARIADNE_TEST_ASSERT(decide(generated_real[1]==1));
    PointCanvas canvas;
    pt.draw(canvas,Projection2d(2u,1u,0u));
    ARIADNE_TEST_EQUALS(canvas.dots,1u);
    ARIADNE_TEST_EQUALS(canvas.strokes,1u);
    ARIADNE_TEST_EQUALS(canvas.x,2.0);
    ARIADNE_TEST_EQUALS(canvas.y,1.0);
    auto bounds=pt.bounding_box();
    ARIADNE_TEST_EQUALS(bounds.dimension(),2u);
    ARIADNE_TEST_PRINT(bounds);
    ARIADNE_TEST_ASSERT(definitely(subset(ExactBoxType{{1.0_x,1.0_x},{2.0_x,2.0_x}},cast_exact_box(bounds))));
    ARIADNE_TEST_EQUALS(zero.bounding_box().dimension(),0u);
    std::unique_ptr<Point<FloatDPApproximation>> clone(pt.clone());
    ARIADNE_TEST_ASSERT(decide(*clone==pt));
    ARIADNE_TEST_ASSERT(decide(clone->centre()==pt.vector()));
    std::ostringstream repr;
    pt._write(repr);
    ARIADNE_TEST_ASSERT(!repr.str().empty());
    auto joined=product(pt,pt);
    ARIADNE_TEST_EQUALS(joined.dimension(),4u);
    ARIADNE_TEST_ASSERT(decide(project(joined,Projection2d(4u,2u,3u))==pt));
    auto origin=Point<FloatDPApproximation>::origin(2u,dp);
    ARIADNE_TEST_ASSERT(std::abs(distance(pt,origin).raw().get_d()-std::sqrt(5.0))<1e-12);
    Point<FloatDPApproximation> line({3.0_x},dp);
    Point<FloatDPApproximation> line_origin({0.0_x},dp);
    ARIADNE_TEST_EQUALS(distance(line,line_origin).raw().get_d(),3.0);
    ARIADNE_TEST_FAIL(distance(pt,line));
    ARIADNE_TEST_FAIL(distance(zero,zero));
    ARIADNE_TEST_FAIL(project(pt,Projection2d(3u,0u,1u)));
}

Int main() {
    ARIADNE_TEST_CALL(test_point_geometry());
    ARIADNE_TEST_CALL(test_point_class_names());
    return ARIADNE_TEST_FAILURES;
}
