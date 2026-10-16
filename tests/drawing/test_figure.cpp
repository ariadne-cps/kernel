/***************************************************************************
 *            test_graphics.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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
#include <sstream>

#include "macro/testing.hpp"

#include "function/function.hpp"
#include "drawing/figure.hpp"
#include "drawing/embedded.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include "drawing/cairo.hpp"
#endif
#include "drawing/null_graphics.hpp"
#include "drawing/drawer.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "geometry/curve.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_paving.hpp"
#include "logging/logging.hpp"

using namespace Ariadne;

class TestDrawable final : public Drawable2dInterface {
  public:
    explicit TestDrawable(DimensionType dimension) : _dimension(dimension) { }
    TestDrawable* clone() const override { return new TestDrawable(*this); }
    Void draw(CanvasInterface&, const Projection2d&) const override { }
    DimensionType dimension() const override { return _dimension; }
  private:
    DimensionType _dimension;
};

class RecordingCanvas final : public NullCanvas {
  public:
    Void initialise(StringType x, StringType y, double, double, double, double) override {
        x_label=x;
        y_label=y;
    }

    StringType x_label;
    StringType y_label;
};

class TestGraphics {
  public:
    void test() {
        ARIADNE_TEST_CALL(test_point2d())
        ARIADNE_TEST_CALL(test_geometry2d())
        ARIADNE_TEST_CALL(test_colour())
        ARIADNE_TEST_CALL(test_graphics_properties())
        ARIADNE_TEST_CALL(test_figure_header_helpers())
        ARIADNE_TEST_CALL(test_graphics_interface())
        ARIADNE_TEST_CALL(test_graphics_base())
        ARIADNE_TEST_CALL(test_embedded())
#ifdef ARIADNE_HAS_CAIRO
        ARIADNE_TEST_CALL(test_cairo())
#endif
        ARIADNE_TEST_CALL(test_construct_figure())
        ARIADNE_TEST_CALL(test_figure_projection())
        ARIADNE_TEST_CALL(test_figure_public_api())
        ARIADNE_TEST_CALL(test_figure_edge_paths())
        ARIADNE_TEST_CALL(test_default_drawer())
    }

    void test_point2d() {
        ARIADNE_TEST_PRINT(Point2d(1,1))
    }

    void test_geometry2d() {
        Vector2d v1(1.0,2.0);
        Vector2d v2(3.0,4.0);

        Vector2d neg = -v1;
        ARIADNE_TEST_EQUALS(neg.x,-1.0)
        ARIADNE_TEST_EQUALS(neg.y,-2.0)

        Vector2d sum = v1+v2;
        ARIADNE_TEST_EQUALS(sum.x,4.0)
        ARIADNE_TEST_EQUALS(sum.y,6.0)

        Vector2d difference = v2-v1;
        ARIADNE_TEST_EQUALS(difference.x,2.0)
        ARIADNE_TEST_EQUALS(difference.y,2.0)

        Vector2d scaled = 2.0*v1;
        ARIADNE_TEST_EQUALS(scaled.x,2.0)
        ARIADNE_TEST_EQUALS(scaled.y,4.0)

        std::ostringstream vector_stream;
        vector_stream << v1;
        ARIADNE_TEST_EQUALS(vector_stream.str(),"[1,2]")

        Point2d default_point;
        default_point.x=0.0;
        default_point.y=0.0;
        ARIADNE_TEST_EQUALS(default_point.x,0.0)
        ARIADNE_TEST_EQUALS(default_point.y,0.0)

        Point2d point(1.0,2.0);
        point[0]=3.0;
        Point2d const& const_point=point;
        ARIADNE_TEST_EQUALS(const_point[0],3.0)

        Point2d equal_point(3.0,2.0);
        Point2d different_y(3.0,4.0);
        Point2d different_x(4.0,2.0);
        ARIADNE_TEST_ASSERT(point==equal_point)
        ARIADNE_TEST_ASSERT(!(point==different_y))
        ARIADNE_TEST_ASSERT(!(point==different_x))

        point-=Vector2d(1.0,1.0);
        ARIADNE_TEST_EQUALS(point.x,2.0)
        ARIADNE_TEST_EQUALS(point.y,1.0)

        Box2d box;
        box.xl=-1.0; box.xu=1.0; box.yl=-2.0; box.yu=2.0;
        std::ostringstream box_stream;
        box_stream << box;
        ARIADNE_TEST_EQUALS(box_stream.str(),"[-1,1]x[-2,2]")

        Polytope2d polytope(List<Point2d>{Point2d(0.0,0.0),Point2d(1.0,1.0)});
        ARIADNE_TEST_EQUALS(polytope.number_of_vertices(),2u)
        polytope.vertices()[0]=Point2d(2.0,3.0);
        Polytope2d const& const_polytope=polytope;
        ARIADNE_TEST_EQUALS(const_polytope.vertices().size(),2u)
        ARIADNE_TEST_EQUALS(const_polytope.vertex(0).x,2.0)
        ARIADNE_TEST_EQUALS(const_polytope.vertex(0).y,3.0)
    }

    void test_colour() {
        ARIADNE_TEST_PRINT(Colour(1,1,1,0))
        ARIADNE_TEST_PRINT(Colour("white",1,1,1))
    }

    void test_graphics_properties() {
        GraphicsProperties p;
        p.set_dot_radius(2.0)
         .set_line_colour(0,0,1)
         .set_fill_opacity(0.5)
         .set_fill_colour(0,1,0);
        ARIADNE_TEST_EQUALS(p.dot_radius,2.0)
        ARIADNE_TEST_PRINT(p)
    }

    void test_figure_header_helpers() {
        LineStyle line_style_value = line_style(false);
        ARIADNE_TEST_EQUALS(static_cast<Bool>(line_style_value),false)

        LineWidth line_width_value = line_width(2.0);
        ARIADNE_TEST_EQUALS(static_cast<Dbl>(line_width_value),2.0)

        LineColour line_colour_value = line_colour(black);
        LineColour direct_line_colour(0.1,0.2,0.3);

        FillStyle fill_style_value = fill_style(false);
        ARIADNE_TEST_EQUALS(static_cast<Bool>(fill_style_value),false)

        FillOpacity fill_opacity_value = fill_opacity(0.25);
        ARIADNE_TEST_EQUALS(static_cast<Dbl>(fill_opacity_value),0.25)

        FillColour direct_fill_colour(0.4,0.5,0.6);

        Animated animated_value = set_animated(true);
        ARIADNE_TEST_EQUALS(static_cast<Bool>(animated_value),true)

        GraphicsProperties properties(false,2.0,3.0,black,false,white);
        ARIADNE_TEST_EQUALS(properties.line_style,false)
        ARIADNE_TEST_EQUALS(properties.line_width,2.0)
        ARIADNE_TEST_EQUALS(properties.dot_radius,3.0)
        ARIADNE_TEST_EQUALS(properties.fill_style,false)

        Interval<Double> double_interval(1.0,2.0);
        ARIADNE_TEST_EQUALS(double_interval.lower_bound(),1.0)
        ARIADNE_TEST_EQUALS(double_interval.upper_bound(),2.0)

        ApproximateDoubleInterval approximate_interval(1.0,2.0);
        ARIADNE_TEST_PRINT(approximate_interval.lower_bound())
        ARIADNE_TEST_PRINT(approximate_interval.upper_bound())

        GraphicsBoundingBoxType bounding_box({
            {FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)},
            {FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)}
        });
        Figure figure(bounding_box,0,1);
        figure << line_style_value
               << line_width_value
               << direct_line_colour
               << fill_style_value
               << fill_opacity_value
               << direct_fill_colour
               << animated_value;

    }

    void test_graphics_interface() {
        GraphicsBoundingBoxType bx({{FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)},
                                    {FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)}});
        Figure figure(bx,0,1);
        Polytope2d shape(List<Point2d>{Point2d(0.0,0.0)});
        FigureInterface& figure_interface = figure;
        Drawable2dInterface const& drawable = shape;
        draw(figure_interface,drawable);
    }

    void test_graphics_base() {
        EmbeddedCanvas canvas(ImageSize2d(32u,32u));
        List<Point2d> boundary{Point2d(0.0,0.0)};
        canvas.fill_boundary(boundary);
    }

    void test_embedded() {
        ARIADNE_TEST_FAIL(ImageSize2d(0,1))
        ARIADNE_TEST_FAIL(ImageSize2d(1,0))
        EmbeddedCanvas canvas(ImageSize2d(32u,32u));
        canvas.initialise("x","y",-1.0,1.0,-1.0,1.0);
        canvas.set_dot_radius(4.0);
        canvas.circle(0.0,0.0,0.25);
        canvas.dot(0.5,0.5);
        canvas.stroke();
        canvas.finalise();
        canvas.write("test_embedded-explicit.png");
    }

#ifdef ARIADNE_HAS_CAIRO
    void test_cairo() {
        ARIADNE_TEST_FAIL(ImageSize2d(0,1))
        ARIADNE_TEST_FAIL(ImageSize2d(1,0))

        CairoCanvas canvas(ImageSize2d(32u,32u));
        canvas.initialise("x","y",-1.0,1.0,-1.0,1.0);
        canvas.set_dot_radius(4.0);
        canvas.circle(0.0,0.0,0.25);
        canvas.dot(0.5,0.5);
        canvas.finalise();
        canvas.write("test_cairo-explicit.png");
    }

#endif

    void test_construct_figure() {
        GraphicsBoundingBoxType bx3({{FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)},
                                    {FloatDPApproximation(-1.0,dp),FloatDPApproximation(0.0,dp)},
                                    {FloatDPApproximation(-2.0,dp),FloatDPApproximation(2.0,dp)}});
        ARIADNE_TEST_EXECUTE(Figure(bx3,0,1))
        ARIADNE_TEST_EXECUTE(Figure(bx3,{0,1}))
        ARIADNE_TEST_EXECUTE(Figure(bx3,Projection2d(3,0,1)))
        GraphicsBoundingBoxType bx1({{FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)}});
        ARIADNE_TEST_FAIL(Figure(bx1,0,1))
        ARIADNE_TEST_FAIL(Figure(bx1,{0,1}))
        ARIADNE_TEST_FAIL(Figure(bx1,Projection2d(2,0,1)))
        Figure fig(bx3,0,1);
        ARIADNE_TEST_PRINT(fig.get_bounding_box())
        auto const& p = fig.properties();
        ARIADNE_TEST_PRINT(p)
    }

    void test_figure_projection() {
        GraphicsBoundingBoxType bx({{FloatDPApproximation(0.0,dp),FloatDPApproximation(1.0,dp)},
                                    {FloatDPApproximation(-1.0,dp),FloatDPApproximation(0.0,dp)},
                                    {FloatDPApproximation(-2.0,dp),FloatDPApproximation(2.0,dp)}});
        Figure f(bx,0,1);
        ARIADNE_TEST_FAIL(f.set_projection(4,0,1))
        ARIADNE_TEST_EXECUTE(f.set_projection(3,0,1))
        ARIADNE_TEST_FAIL(f.set_projection_map(Projection2d(4,0,1)))
        ARIADNE_TEST_EXECUTE(f.set_projection_map(Projection2d(3,0,1)))
        ARIADNE_TEST_PRINT(f.get_2d_projection_map())
    }

    void test_figure_public_api() {
        GraphicsBoundingBoxType bx({
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)},
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)}
        });
        Figure figure(bx,0,1);

        figure.set_dot_radius(3.0)
              .set_line_style(false)
              .set_fill_style(false)
              .set_animated(true);

        ARIADNE_TEST_EQUALS(figure.get_dot_radius(),3.0)
        ARIADNE_TEST_EQUALS(figure.get_line_style(),false)
        ARIADNE_TEST_EQUALS(figure.get_line_width(),1.0)
        ARIADNE_TEST_EQUALS(figure.get_fill_style(),false)
        ARIADNE_TEST_EQUALS(figure.get_fill_opacity(),1.0)
        ARIADNE_TEST_EQUALS(figure.get_line_colour().red,0.0)
        ARIADNE_TEST_EQUALS(figure.get_fill_colour().red,1.0)

        Figure const& const_figure=figure;
        ARIADNE_TEST_EQUALS(const_figure.properties().is_animated,true)

        RealBox real_box({{-0.5_dec,0.5_dec},{-0.25_dec,0.25_dec}});
        draw(figure,real_box);
        figure.draw(real_box);

        figure.set_axis_labels("horizontal","vertical");
        RecordingCanvas canvas;
        figure._paint_all(canvas);
        ARIADNE_TEST_EQUALS(canvas.x_label,StringType("horizontal"))
        ARIADNE_TEST_EQUALS(canvas.y_label,StringType("vertical"))
    }

void test_figure_edge_paths() {
        GraphicsBoundingBoxType bx2({
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)},
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)}
        });
        GraphicsBoundingBoxType bx3({
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)},
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)},
            {FloatDPApproximation(-1.0,dp),FloatDPApproximation(1.0,dp)}
        });
        NullCanvas canvas;

        Figure fallback(bx2,0,1);
        fallback.set_bounding_box(GraphicsBoundingBoxType(0));
        fallback._paint_all(canvas);

        Figure mismatched_box(bx2,0,1);
        mismatched_box.set_bounding_box(bx3);
        ARIADNE_TEST_FAIL(mismatched_box._paint_all(canvas))

        Figure invalid_x(bx2,Projection2d(2,2,1));
        ARIADNE_TEST_FAIL(invalid_x._paint_all(canvas))

        Figure invalid_y(bx2,Projection2d(2,0,2));
        ARIADNE_TEST_FAIL(invalid_y._paint_all(canvas))

        Figure empty_shape_figure(bx2,0,1);
        TestDrawable empty_shape(0u);
        empty_shape_figure.draw(empty_shape);
        empty_shape_figure._paint_all(canvas);

        Figure mismatched_shape_figure(bx2,0,1);
        TestDrawable mismatched_shape(3u);
        mismatched_shape_figure.draw(mismatched_shape);
        ARIADNE_TEST_FAIL(mismatched_shape_figure._paint_all(canvas))
    }

void test_default_drawer()
    {
        ExactBoxType bx1(2); bx1[0]=ExactIntervalType(-0.2_pr,0.2_pr); bx1[1]=ExactIntervalType(-0.1_pr,0.10_pr);
        ExactBoxType bx2(2); bx2[0]=ExactIntervalType(0.1_pr,0.3_pr); bx2[1]=ExactIntervalType(0.05_pr,0.15_pr);
        ExactBoxType bx3(2); bx3[0]=ExactIntervalType(0.2_pr,0.4_pr); bx3[1]=ExactIntervalType(0.10_pr,0.25_pr);
        ExactBoxType bx4(2); bx4[0]=ExactIntervalType(0.25_pr,0.5_pr); bx4[1]=ExactIntervalType(0.20_pr,0.50_pr);
        ExactBoxType bx5(2); bx5[0]=ExactIntervalType(0.4_pr,0.8_pr); bx5[1]=ExactIntervalType(0.40_pr,1.1_pr);

        Real p(0.5_x);
        EffectiveVectorMultivariateFunction x=EffectiveVectorMultivariateFunction::identity(3);
        EffectiveVectorMultivariateFunction afn1={0.05_dec*x[0]+0.05_dec*x[2]+0.15_dec,0.05_dec*x[1]+0.05_dec*x[2]+0.6_dec};
        ValidatedConstrainedImageSet s1(ExactBoxType::unit_box(3),afn1);
        ApproximateBoxType bbx1=widen(s1.bounding_box(),0.25_x);

        EffectiveVectorMultivariateFunction rf(1u, sqr(x[0])+sqr(x[1])-sqr(p));
        ConstraintSet cs1(rf,RealBox(1u,RealInterval(-1,0)));

        {
            ExactDouble h=10000;
            Figure g(ExactBoxType{{-1.0_x,+1.0_x},{-h,+h}},Projection2d(2,0,1));
            g.set_fill_colour(0.5,1.0,1.0);
            g.set_line_width(10);
            g << ExactBoxType({{-0.5_x,+0.0_x},{-h/two, +h/two}});
            g.set_line_width(1);
            g << ExactBoxType({{0.25_x,+0.75_x},{-h/two, +h/two}});
            g.write("test_graphics-canvas");
        }

        Figure g(ApproximateBoxType({{-1,+1},{-1,+1}}),Projection2d(2,0,1));
        g << fill_colour(0.5,1.0,1.0)
          << line_colour(0.0,0.0,0.0)
          << bx1
          << bx2
          << bx3
          << bx4
          << bx5;
        g.write("test_graphics-bx1");
        g.clear();

        g.set_fill_colour(1.0,0.5,1.0);
        g.draw(bx1);
        g.draw(bx2);
        g.set_fill_colour(magenta);
        g.draw(bx5);
        g.write("test_graphics-bx2");
        g.clear();


        ExactBoxType bx2d(2); bx2d[0]=ExactIntervalType(0.2_pr,0.4_pr); bx2d[1]=ExactIntervalType(0.2_pr,0.5_pr);
        ExactBoxType bx3d(3); bx3d[0]=ExactIntervalType(0.2_pr,0.4_pr); bx3d[1]=ExactIntervalType(0.2_pr,0.5_pr); bx3d[2]=ExactIntervalType(0.2_pr,0.7_pr);
        g.set_bounding_box(bx3d.bounding_box());
        g.set_projection(3,0,1);
        g.draw(bx3d);
        g.write("test_graphics-bx3");
        g.clear();

        plot("test_graphics_set",Projection2d(2,0,1),bbx1,{{orange,s1}});

        InterpolatedCurve cv(0,Point<FloatDP>(2,FloatDP(0.0_x,dp)));
        for(Int i=1; i<=10; ++i) {
            Point<FloatDP> pt(2,dp); pt[0]=FloatDP(cast_exact(i/10.),dp); pt[1]=FloatDP(cast_exact(i*i/100.),dp);
            cv.insert(i,pt);
        }
        g.set_bounding_box(cv.bounding_box());
        g.set_projection(2,0,1);
        g.set_line_colour(1,0,0);
        g.draw(cv);
        g.set_line_colour(0,0,0);
        g.write("test_graphics-curve");
        g.clear();

        GridTreePaving gts(2);
        gts.adjoin_outer_approximation(bx1, 4);
        gts.adjoin_outer_approximation(bx2, 5);
        gts.adjoin_outer_approximation(bx3, 6);
        gts.adjoin_outer_approximation(bx4, 7);
        gts.adjoin_outer_approximation(bx5, 8);
        gts.recombine();

        ExactBoxType bbox(2); bbox[0]=ExactIntervalType(-2,2); bbox[1]=ExactIntervalType(-2,2);
        g.set_bounding_box(bbox);
        g.set_projection(2,0,1);
        g << gts;
        g.write("test_graphics-paving");
        g.clear();
    }
};

Int main()
{
    TestGraphics().test();

    return ARIADNE_TEST_FAILURES;
}
