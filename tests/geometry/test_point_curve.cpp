/***************************************************************************
 *            test_point_curve.cpp
 *
 *  Copyright  2005-21  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/module.hpp"
#include "geometry/point.hpp"
#include "labelling/point.hpp"
#include "labelling/curve.hpp"
#include "geometry/curve.hpp"
#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "drawing/figure.hpp"
#include "drawing/null_graphics.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;

class TestCurveCanvas : public NullCanvas {
  public:
    Void move_to(double, double) override { moved=true; }
    Void line_to(double, double) override { lined=true; }

    Bool moved=false;
    Bool lined=false;
};

class TestPointCanvas : public NullCanvas {
  public:
    Void dot(double, double) override { dotted=true; }
    Void stroke() override { stroked=true; }

    Bool dotted=false;
    Bool stroked=false;
};

template<class X> class TestPoint {
    typedef LabelledPoint<X> PointType;
  public:
    void test() {
        ARIADNE_TEST_CALL(test_creation());
        ARIADNE_TEST_CALL(test_projection());
    }

    void test_creation() {
        RealVariable x("x"), y("y");
        LabelledPoint<Real> rpt({x=3,y=2.1_dec});
        PointType pt(rpt,dp);

        ARIADNE_TEST_EQUALS(class_name<LabelledPoint<Real>>(),String("LabelledPoint<Real>"));
        ARIADNE_TEST_EQUALS(class_name<LabelledPoint<FloatDPApproximation>>(),String("LabelledPoint<FloatDPApproximation>"));

        ARIADNE_TEST_EQUALS(pt.dimension(),2);
        ARIADNE_TEST_EQUAL(pt.space(),RealSpace({x,y}));
        ARIADNE_TEST_EQUAL(pt.variables(),Set<RealVariable>(RealSpace({x,y}).variables()));
        ARIADNE_TEST_PRINT(pt.bounding_box());

        std::unique_ptr<PointType> clone(pt.labelled_clone());
        ARIADNE_TEST_EQUALS(clone->dimension(),pt.dimension());
        ARIADNE_TEST_EQUAL(clone->space(),pt.space());
        ARIADNE_TEST_EQUAL(clone->variables(),pt.variables());

        TestPointCanvas canvas;
        pt.draw(canvas,Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.dotted);
        ARIADNE_TEST_ASSERT(canvas.stroked);
    }

    void test_projection() {
        RealVariable x("x"), y("y"), z("z");
        LabelledPoint<Real> rpt({x=3,y=2.1_dec,z=1/5_q});
        PointType pt(rpt,dp);
        ARIADNE_TEST_PRINT(pt);
        ARIADNE_TEST_EQUALS(pt.dimension(),3);

        auto pt_pr = project(pt,Variables2d(z,x));
        ARIADNE_TEST_PRINT(pt_pr);
        ARIADNE_TEST_EQUALS(pt_pr.dimension(),2);
        ARIADNE_TEST_EQUAL(pt.at(0),pt_pr.at(1));
        ARIADNE_TEST_EQUAL(pt.at(2),pt_pr.at(0));
    }
};

class TestCurve {
    typedef LabelledInterpolatedCurve CurveType;
    typedef CurveType::PointType PointType;
  public:
    void test() {
        ARIADNE_TEST_CALL(test_creation());
    }

    void test_creation() {
        RealVariable x("x"), y("y");
        LabelledPoint<Real> rpt({x=3,y=2.1_dec});
        PointType pt(rpt,dp);
        CurveType curve(pt);
        ARIADNE_TEST_EQUALS(curve.dimension(),2);
        ARIADNE_TEST_EQUAL(curve.space(),pt.space());
        ARIADNE_TEST_EQUAL(curve.variables(),pt.variables());
        ARIADNE_TEST_EQUALS(curve.size(),1);
        curve.insert(FloatDP(1.0_x,dp),pt);
        ARIADNE_TEST_EQUALS(curve.size(),2);

        std::unique_ptr<CurveType> clone(curve.labelled_clone());
        ARIADNE_TEST_EQUALS(clone->dimension(),curve.dimension());
        ARIADNE_TEST_EQUAL(clone->space(),curve.space());
        ARIADNE_TEST_EQUALS(clone->size(),curve.size());

        TestCurveCanvas canvas;
        curve.draw(canvas,Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.moved);
        ARIADNE_TEST_ASSERT(canvas.lined);
    }
};


int main() {
    TestPoint<FloatDPApproximation>().test();
    TestCurve().test();

    return ARIADNE_TEST_FAILURES;
}
