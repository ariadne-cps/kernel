/***************************************************************************
 *            tests/drawing/test_default_graphics_robustness.cpp
 *
 *  Copyright  2026  Luca Geretti
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

/*! \file tests/drawing/test_default_graphics_robustness.cpp
 *  \brief Regression checks for invalid input, opacity and plot clipping.
 */

#include "drawing/default_graphics.hpp"
#include "default_graphics_png.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace Ariadne;

namespace {

//! \brief Require the supplied operation to reject an invalid argument.
template<class Function>
void must_reject(Function function,const char* scenario) {
    try {
        function();
    } catch(const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(scenario);
}

//! \brief Signed image dimensions must fail the ImageSize2d precondition.
template<class Function>
void must_fail_precondition(Function function,const char* scenario) {
    try {
        function();
    } catch(const std::runtime_error& error) {
        // Other runtime failures are not evidence that the size guard ran.
        const std::string message=error.what();
        if(message.find("Precondition")!=std::string::npos &&
           message.find("x > 0 && y > 0")!=std::string::npos)
            return;
        throw;
    }
    throw std::runtime_error(scenario);
}

//! \brief Append and fill a rectangular contour in plot coordinates.
void rectangle(DefaultGraphicsCanvas& canvas,double xl,double yl,
               double xu,double yu) {
    canvas.move_to(xl,yl);
    canvas.line_to(xu,yl);
    canvas.line_to(xu,yu);
    canvas.line_to(xl,yu);
    canvas.fill();
}

//! \brief Verify that each embedded font glyph is visibly rendered.
void test_embedded_axis_glyphs() {
    constexpr unsigned plot_width=520u;
    constexpr unsigned plot_height=100u;
    constexpr int advance=18;
    const auto render_label=[=](std::string const& label,const char* filename) {
        DefaultGraphicsCanvas canvas(ImageSize2d(plot_width,plot_height));
        canvas.initialise(label,"",-1.0,1.0,-1.0,1.0);
        canvas.write(filename);
        return DrawingTest::read_default_graphics_png(filename);
    };
    const auto verify_cells=[=](DrawingTest::Image const& text,
                                DrawingTest::Image const& blank,
                                std::string const& label) {
        const int x0=160+static_cast<int>(plot_width)/2-
                     static_cast<int>(label.size())*advance/2-3;
        const int y0=10+static_cast<int>(plot_height)+4;
        for(std::size_t i=0;i<label.size();++i) {
            bool painted=false;
            for(int y=y0;y<y0+21;++y)
                for(int x=x0+static_cast<int>(i)*advance;
                    x<x0+static_cast<int>(i)*advance+15;++x)
                    if(text.pixel(static_cast<unsigned>(x),static_cast<unsigned>(y)) !=
                       blank.pixel(static_cast<unsigned>(x),static_cast<unsigned>(y)))
                        painted=true;
            DrawingTest::require(painted==(label[i]!=' '),
                                 "Embedded font glyph missing at index "+std::to_string(i));
        }
    };

    const std::string lower="abcdefghijklmnopqrstuvwxyz";
    const std::string upper="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const auto blank_letters=render_label(std::string(lower.size(),' '),
                                          "default_graphics-blank-letters.png");
    const auto lower_image=render_label(lower,"default_graphics-lowercase.png");
    const auto upper_image=render_label(upper,"default_graphics-uppercase.png");
    verify_cells(lower_image,blank_letters,lower);
    DrawingTest::require(lower_image.pixels==upper_image.pixels,
                         "Uppercase and lowercase glyphs render differently");

    // Digits, punctuation, whitespace and unsupported '?' fallback.
    const std::string symbols="0123456789+-.,:_/()? ";
    const auto blank_symbols=render_label(std::string(symbols.size(),' '),
                                          "default_graphics-blank-symbols.png");
    const auto symbol_image=render_label(symbols,"default_graphics-symbols.png");
    verify_cells(symbol_image,blank_symbols,symbols);
}

//! \brief Reinitialising must erase pixels and pending drawing paths.
void test_canvas_reset_and_offscreen_geometry() {
    DefaultGraphicsCanvas reused(ImageSize2d(64u,48u));
    reused.initialise("","",-1.0,1.0,-1.0,1.0);
    reused.set_fill_colour(1.0,0.0,0.0);
    reused.move_to(-0.8,-0.8);
    reused.line_to(0.8,-0.8);
    reused.line_to(0.8,0.8);
    reused.line_to(-0.8,0.8);
    reused.fill();
    // An unconsumed path must also be discarded on reinitialisation.
    reused.move_to(-0.9,-0.9);
    reused.line_to(0.9,0.9);
    reused.initialise("","",-1.0,1.0,-1.0,1.0);

    DefaultGraphicsCanvas clean(ImageSize2d(64u,48u));
    clean.initialise("","",-1.0,1.0,-1.0,1.0);
    clean.write("default_graphics-clean-canvas.png");
    const auto expected=DrawingTest::read_default_graphics_png(
        "default_graphics-clean-canvas.png");

    // Rectangle fast path and generic polygon path must both clip away
    // completely off-screen geometry, without modifying the plot or margins.
    reused.set_line_width(0.0);
    reused.move_to(2.0,2.0);
    reused.line_to(3.0,2.0);
    reused.line_to(3.0,3.0);
    reused.line_to(2.0,3.0);
    reused.fill();
    reused.move_to(2.0,2.0);
    reused.line_to(3.0,2.5);
    reused.line_to(2.5,3.0);
    reused.fill();
    reused.write("default_graphics-reset-and-clipped.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-reset-and-clipped.png");
    DrawingTest::require(actual.pixels==expected.pixels,
                         "Canvas reinitialisation or off-screen fill changed pixels");
}

//! \brief A five-vertex contour is not a rectangle with implicit closure.
void test_nonrectangular_five_vertex_fill() {
    DefaultGraphicsCanvas canvas(ImageSize2d(120u,100u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.set_fill_opacity(1.0);
    // The fifth vertex cuts away a wedge from the left of the rectangle.
    // Treating the first four vertices as a rectangle would paint that wedge.
    canvas.move_to(-0.8,-0.8);
    canvas.line_to(0.8,-0.8);
    canvas.line_to(0.8,0.8);
    canvas.line_to(-0.8,0.8);
    canvas.line_to(0.0,0.0);
    canvas.fill();
    canvas.write("default_graphics-five-vertex-fill.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-five-vertex-fill.png");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.4,-0.5)==
            std::array<unsigned char,3>{255,0,0},
        "Five-vertex contour failed to fill the interior");
    DrawingTest::require(
        DrawingTest::world_pixel(image,-0.7,0.0)==
            std::array<unsigned char,3>{247,247,247},
        "Five-vertex contour incorrectly used the rectangle fast path");
}

//! \brief Rectangles entirely beyond each plot edge must not paint pixels.
void test_clipped_rectangles_on_all_sides() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-clipped-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-clipped-blank.png");

    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    const auto rectangle=[&](double xl,double yl,double xu,double yu) {
        canvas.move_to(xl,yl);
        canvas.line_to(xu,yl);
        canvas.line_to(xu,yu);
        canvas.line_to(xl,yu);
        canvas.fill();
    };
    rectangle(-2.0,-0.4,-1.2,0.4); // Left.
    rectangle(1.2,-0.4,2.0,0.4);   // Right.
    rectangle(-0.4,-2.0,0.4,-1.2); // Bottom.
    rectangle(-0.4,1.2,0.4,2.0);   // Top.
    canvas.write("default_graphics-clipped-all-sides.png");
    const auto clipped=DrawingTest::read_default_graphics_png(
        "default_graphics-clipped-all-sides.png");
    DrawingTest::require(clipped.pixels==blank.pixels,
                         "Off-screen rectangles unexpectedly painted the canvas");
}

//! \brief Two overlapping rectangles must retain the nonzero winding fill.
void test_overlapping_rectangular_contours() {
    struct Overlap {
        double xl,yl,xu,yu;
        double probe_x,probe_y;
    };
    const std::array<Overlap,4> cases{{
        {-0.8,-0.25,-0.2,0.25,-0.7,0.0}, // Inner extends left.
        {0.2,-0.25,0.8,0.25,0.7,0.0},    // Inner extends right.
        {-0.25,-0.8,0.25,-0.2,0.0,-0.7}, // Inner extends below.
        {-0.25,0.2,0.25,0.8,0.0,0.7}     // Inner extends above.
    }};
    for(std::size_t i=0;i<cases.size();++i) {
        DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
        canvas.initialise("","",-1.0,1.0,-1.0,1.0);
        canvas.set_line_width(0.0);
        canvas.set_fill_colour(1.0,0.0,0.0);
        canvas.set_fill_opacity(1.0);
        const auto contour=[&](double xl,double yl,double xu,double yu) {
            canvas.move_to(xl,yl);
            canvas.line_to(xu,yl);
            canvas.line_to(xu,yu);
            canvas.line_to(xl,yu);
        };
        contour(-0.6,-0.6,0.6,0.6);
        const auto& r=cases[i];
        contour(r.xl,r.yl,r.xu,r.yu);
        canvas.fill();
        const auto filename="default_graphics-overlap-"+std::to_string(i)+".png";
        canvas.write(filename.c_str());
        const auto image=DrawingTest::read_default_graphics_png(filename);
        DrawingTest::require(
            DrawingTest::world_pixel(image,r.probe_x,r.probe_y)==
                std::array<unsigned char,3>{255,0,0},
            "Overlapping nonzero-winding contours lost the outside portion");
        DrawingTest::require(
            DrawingTest::world_pixel(image,0.0,0.0)==
                std::array<unsigned char,3>{255,0,0},
            "Overlapping nonzero-winding contours lost the outer interior");
        DrawingTest::require(
            DrawingTest::world_pixel(image,0.95,0.95)==
                std::array<unsigned char,3>{247,247,247},
            "Overlapping contours painted beyond their geometry");
    }
}

//! \brief A redundant collinear vertex must not change a nested-hole fill.
void test_rectangle_fast_path_matches_generic_fill() {
    const auto render=[](bool redundant_vertex,const char* filename) {
        DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
        canvas.initialise("","",-1.0,1.0,-1.0,1.0);
        canvas.set_line_width(0.0);
        canvas.set_fill_colour(1.0,0.0,0.0);
        canvas.set_fill_opacity(1.0);
        canvas.move_to(-0.75,-0.75);
        if(redundant_vertex) canvas.line_to(0.0,-0.75);
        canvas.line_to(0.75,-0.75);
        canvas.line_to(0.75,0.75);
        canvas.line_to(-0.75,0.75);
        canvas.move_to(-0.25,-0.25);
        canvas.line_to(-0.25,0.25);
        canvas.line_to(0.25,0.25);
        canvas.line_to(0.25,-0.25);
        canvas.fill();
        canvas.write(filename);
        return DrawingTest::read_default_graphics_png(filename);
    };
    const auto fast=render(false,"default_graphics-hole-fast.png");
    const auto generic=render(true,"default_graphics-hole-general.png");
    DrawingTest::require(fast.pixels==generic.pixels,
                         "Nested-hole fill changed when using the generic rasterizer");
}

//! \brief Empty, two-point and collinear contours never cover an area.
void test_degenerate_fills_leave_pixels_unchanged() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-degenerate-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-degenerate-blank.png");
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.fill();
    canvas.move_to(-0.5,-0.5);
    canvas.line_to(0.5,0.5);
    canvas.fill();
    canvas.move_to(-0.6,0.0);
    canvas.line_to(0.0,0.0);
    canvas.line_to(0.6,0.0);
    canvas.fill();
    canvas.write("default_graphics-degenerate-fills.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-degenerate-fills.png");
    DrawingTest::require(image.pixels==blank.pixels,
                         "Degenerate path unexpectedly painted the canvas");
}

//! \brief Distinct paint colours retain their values after cache eviction.
void test_paint_cache_survives_eviction() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(0.0);
    canvas.set_fill_opacity(1.0);
    const std::array<std::array<double,3>,8> colours{{
        {{1.0,0.0,0.0}}, {{0.0,1.0,0.0}}, {{0.0,0.0,1.0}},
        {{1.0,1.0,0.0}}, {{1.0,0.0,1.0}}, {{0.0,1.0,1.0}},
        {{1.0,0.0,0.0}}, {{0.0,0.0,1.0}}
    }};
    for(std::size_t i=0;i<colours.size();++i) {
        const double cx=-0.84+0.24*static_cast<double>(i);
        canvas.set_fill_colour(colours[i][0],colours[i][1],colours[i][2]);
        rectangle(canvas,cx-0.08,-0.5,cx+0.08,0.5);
    }
    canvas.write("default_graphics-paint-cache.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-paint-cache.png");
    for(std::size_t i=0;i<colours.size();++i) {
        const double cx=-0.84+0.24*static_cast<double>(i);
        const std::array<unsigned char,3> expected{{
            static_cast<unsigned char>(255.0*colours[i][0]),
            static_cast<unsigned char>(255.0*colours[i][1]),
            static_cast<unsigned char>(255.0*colours[i][2])
        }};
        DrawingTest::require(DrawingTest::world_pixel(image,cx,0.0)==expected,
                             "Cached paint colour changed at index "+std::to_string(i));
    }
}

//! \brief A rectangle narrower than the subpixel sampling grid paints nothing.
void test_subpixel_empty_rectangle() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-subpixel-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-subpixel-baseline.png");
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    rectangle(canvas,-0.00002,-0.5,0.00002,0.5);
    canvas.write("default_graphics-subpixel-empty.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-subpixel-empty.png");
    DrawingTest::require(image.pixels==blank.pixels,
                         "Rectangle with no subpixel samples unexpectedly painted");
}

//! \brief Overflowed device-space circle geometry is ignored.
void test_overflowed_circle_geometry() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-overflow-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-overflow-baseline.png");
    canvas.circle(0.0,0.0,std::numeric_limits<double>::max());
    canvas.circle(std::numeric_limits<double>::max(),0.0,0.25);
    canvas.set_dot_radius(-1.0);
    canvas.dot(0.0,0.0);
    canvas.stroke();
    canvas.write("default_graphics-overflow-circle.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-overflow-circle.png");
    DrawingTest::require(image.pixels==blank.pixels,
                         "Overflowed circle geometry unexpectedly painted");
}

//! \brief Non-finite paint channels and opacity have defined, safe semantics.
void test_nonfinite_paint_values() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(0.0);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double positive_infinity=std::numeric_limits<double>::infinity();

    // Undefined or infinite colour channels normalize to zero.
    canvas.set_fill_colour(nan,positive_infinity,1.0);
    canvas.set_fill_opacity(1.0);
    rectangle(canvas,-0.8,-0.5,-0.2,0.5);

    // An undefined alpha contributes zero coverage.
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.set_fill_opacity(nan);
    rectangle(canvas,0.2,-0.5,0.8,0.5);
    canvas.write("default_graphics-nonfinite-paint.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-nonfinite-paint.png");
    DrawingTest::require(
        DrawingTest::world_pixel(image,-0.5,0.0)==
            std::array<unsigned char,3>{0,0,255},
        "Non-finite RGB channels were not normalized to zero");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.5,0.0)==
            std::array<unsigned char,3>{247,247,247},
        "Non-finite fill opacity painted the background");
}

//! \brief The first line_to starts a path when move_to has not been called.
void test_line_to_without_move_to() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(3.0);
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.line_to(-0.75,0.0);
    canvas.line_to(0.75,0.0);
    canvas.stroke();
    canvas.write("default_graphics-line-without-move.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-line-without-move.png");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.0,0.0)==
            std::array<unsigned char,3>{255,0,0},
        "line_to failed to start and extend an empty path");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.0,0.75)==
            std::array<unsigned char,3>{247,247,247},
        "line_to painted pixels outside its stroke");
}

//! \brief Strokes entirely outside any of the four plot edges paint nothing.
void test_clipped_strokes_outside_plot() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-stroke-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-stroke-blank.png");
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.set_line_width(4.0);
    const auto segment=[&](double x0,double y0,double x1,double y1) {
        canvas.move_to(x0,y0);
        canvas.line_to(x1,y1);
        canvas.stroke();
    };
    segment(-4.0,-0.8,-3.0,0.8);
    segment(3.0,-0.8,4.0,0.8);
    segment(-0.8,3.0,0.8,4.0);
    segment(-0.8,-3.0,0.8,-4.0);
    segment(-3.0,-0.5,-3.0,0.5);
    segment(-0.5,3.0,0.5,3.0);
    canvas.write("default_graphics-stroke-offscreen.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-stroke-offscreen.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "Completely off-screen strokes changed the raster");
}

//! \brief Filling across every boundary is equivalent to filling its clip.
void test_clipped_rectangle_matches_plot_bounds() {
    const auto render=[](double low,double high,const char* file) {
        DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
        canvas.initialise("","",-1.0,1.0,-1.0,1.0);
        canvas.set_line_width(0.0);
        canvas.set_fill_colour(0.2,0.6,0.9);
        canvas.set_fill_opacity(0.65);
        rectangle(canvas,low,low,high,high);
        canvas.write(file);
        return DrawingTest::read_default_graphics_png(file);
    };
    const auto clipped=render(-3.0,3.0,"default_graphics-rectangle-overflow.png");
    const auto bounded=render(-1.0,1.0,"default_graphics-rectangle-bounded.png");
    DrawingTest::require(clipped.pixels==bounded.pixels,
                         "Rectangular fill clipping changes interior or margins");
}

//! \brief A near-reversing stroke must obey the bounded miter limit.
void test_acute_stroke_miter_limit() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(8.0);
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.move_to(-0.7,0.0);
    canvas.line_to(0.0,0.0);
    canvas.line_to(-0.7,0.02);
    canvas.stroke();
    canvas.write("default_graphics-bounded-miter.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-bounded-miter.png");
    // Butt caps end at the geometric vertex (device x=224), while the
    // first subpixel sample in that column lies at x=224.125. Check
    // fully covered stroke interior instead of requiring an opaque endpoint.
    DrawingTest::require(
        DrawingTest::world_pixel(image,-0.1,0.0)==
            std::array<unsigned char,3>{255,0,0},
        "Near-reversing join lost the stroke approaching its vertex");

    // For this nearly reversed path the unlimited miter would project
    // horizontally to the right, not vertically above the vertex.
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.8,0.0)==
            std::array<unsigned char,3>{247,247,247},
        "Near-reversing join exceeded the miter length limit");
}

//! \brief Degenerate rectangle-like paths must not gain artificial area.
void test_rectangle_recognition_rejects_degenerate_contours() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-rectangle-degenerate-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-rectangle-degenerate-baseline.png");
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    rectangle(canvas,-0.5,0.0,0.5,0.0);
    rectangle(canvas,0.0,-0.5,0.0,0.5);

    // The second and fourth vertices coincide.  A four-corner fast path
    // must not treat this retraced contour as an ordinary rectangle.
    canvas.move_to(-0.7,-0.7);
    canvas.line_to(0.7,-0.7);
    canvas.line_to(0.7,0.7);
    canvas.line_to(0.7,-0.7);
    canvas.fill();
    canvas.write("default_graphics-rectangle-degenerate.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-rectangle-degenerate.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "Degenerate or retraced contour unexpectedly filled");
}

//! \brief Non-axis quadrilaterals and mixed subpaths use the general filler.
void test_mixed_nonrectangular_contours() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.set_line_width(0.0);
    canvas.move_to(-0.8,-0.4);
    canvas.line_to(-0.2,-0.4);
    canvas.line_to(-0.2,0.4);
    canvas.line_to(-0.8,0.4);
    canvas.move_to(0.15,-0.5);
    canvas.line_to(0.85,-0.5);
    canvas.line_to(0.5,0.5);
    canvas.fill();
    canvas.write("default_graphics-mixed-subpaths.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-mixed-subpaths.png");
    DrawingTest::require(
        DrawingTest::world_pixel(image,-0.5,0.0)==
            std::array<unsigned char,3>{255,0,0},
        "Rectangle in mixed path was not filled");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.5,-0.2)==
            std::array<unsigned char,3>{255,0,0},
        "Triangular subpath was not filled");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.0,0.7)==
            std::array<unsigned char,3>{247,247,247},
        "Mixed subpaths filled outside their contours");
}

//! \brief The paint-cache key includes opacity, even for identical RGB.
void test_same_rgb_different_opacity() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.set_fill_opacity(0.25);
    rectangle(canvas,-0.8,-0.5,-0.2,0.5);
    canvas.set_fill_opacity(0.75);
    rectangle(canvas,0.2,-0.5,0.8,0.5);
    canvas.write("default_graphics-opacity-cache.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-opacity-cache.png");
    DrawingTest::require(
        DrawingTest::world_pixel(image,-0.5,0.0)==
            std::array<unsigned char,3>{249,185,185},
        "Quarter-opacity red used an incorrect cached paint");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.5,0.0)==
            std::array<unsigned char,3>{253,62,62},
        "Three-quarter-opacity red reused the wrong cached paint");
}

//! \brief Empty paths and invalid widths cannot paint a stroke.
void test_empty_and_nonfinite_strokes() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-empty-stroke-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-empty-stroke-baseline.png");
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.set_line_width(3.0);
    canvas.move_to(0.0,0.0);
    canvas.stroke(); // One point is insufficient to form a segment.
    canvas.set_line_width(std::numeric_limits<double>::infinity());
    canvas.move_to(-0.8,0.0);
    canvas.line_to(0.8,0.0);
    canvas.stroke();
    canvas.write("default_graphics-empty-stroke.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-empty-stroke.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "An empty path or infinite width painted a stroke");
}

//! \brief General scanline spans outside either horizontal clip boundary
//! must leave the raster unmodified even when vertically visible.
void test_offscreen_generic_fill_spans() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-offscreen-span-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-offscreen-span-blank.png");
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    const auto triangle=[&](double x0,double x1,double x2,
                            double y0,double y1,double y2) {
        canvas.move_to(x0,y0);
        canvas.line_to(x1,y1);
        canvas.line_to(x2,y2);
        canvas.fill();
    };
    triangle(-3.0,-2.0,-2.5,-0.4,-0.4,0.4);
    triangle(2.0,3.0,2.5,-0.4,-0.4,0.4);
    triangle(-0.4,0.4,0.0,-3.0,-3.0,-2.0);
    canvas.write("default_graphics-offscreen-spans.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-offscreen-spans.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "Offscreen scanline spans painted the plot");
}

//! \brief Finite extreme stroke widths must not overflow the axial fast path
//! into an incorrect image or unsafe memory access.
void test_extreme_axial_stroke_clipping() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-extreme-stroke-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-extreme-stroke-blank.png");

    // The transformed vertical line has device x near 9.6e307. Its
    // 8.9e307 half-width makes xmax overflow to infinity, while the
    // segment length and cross-product bound remain finite. The generic
    // rasterizer must still reject the segment as horizontally offscreen.
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.set_line_width(1.78e308);
    canvas.move_to(1.5e306,-0.015);
    canvas.line_to(1.5e306,0.015);
    canvas.stroke();
    canvas.write("default_graphics-extreme-stroke.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-extreme-stroke.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "Extreme finite axial stroke incorrectly painted the plot");
}

//! \brief Long labels may cross both raster bounds, but never the plot.
void test_long_axis_labels_clip_safely() {
    const auto render=[](StringType const& x_label,const char* filename) {
        DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
        canvas.initialise(x_label,"",-1.0,1.0,-1.0,1.0);
        canvas.write(filename);
        return DrawingTest::read_default_graphics_png(filename);
    };
    const auto blank=render("","default_graphics-short-label.png");
    const auto wide=render(std::string(64,'W'),
                           "default_graphics-overlong-label.png");
    DrawingTest::require(wide.width==blank.width &&
                         wide.height==blank.height,
                         "Long axis labels changed the PNG dimensions");
    DrawingTest::require(wide.pixels!=blank.pixels,
                         "Long axis label was not rendered in the margin");
    for(unsigned y=10;y<106;++y)
        for(unsigned x=160;x<288;++x)
            DrawingTest::require(wide.pixel(x,y)==blank.pixel(x,y),
                                 "Long axis label unexpectedly painted inside the plot");
}

//! \brief Repeated vertices and circles of zero radius have no stroke area.
void test_zero_length_strokes_are_empty() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-zero-stroke-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-zero-stroke-blank.png");
    canvas.set_line_width(4.0);
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.move_to(0.0,0.0);
    canvas.line_to(0.0,0.0);
    canvas.line_to(0.0,0.0);
    canvas.stroke();
    canvas.circle(0.0,0.0,0.0);
    canvas.stroke();
    canvas.set_dot_radius(0.0);
    canvas.dot(0.0,0.0);
    canvas.stroke();
    canvas.write("default_graphics-zero-strokes.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-zero-strokes.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "Zero-length strokes unexpectedly changed the image");
}

//! \brief A slanted quadrilateral must use nonzero-winding rasterisation.
void test_slanted_quadrilateral_fill() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.set_fill_opacity(1.0);
    canvas.move_to(-0.8,0.0);
    canvas.line_to(0.0,0.8);
    canvas.line_to(0.8,0.0);
    canvas.line_to(0.0,-0.8);
    canvas.fill();
    canvas.write("default_graphics-slanted-quad.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-slanted-quad.png");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.0,0.0)==
            std::array<unsigned char,3>{255,0,0},
        "Slanted quadrilateral did not fill its interior");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.9,0.9)==
            std::array<unsigned char,3>{247,247,247},
        "Slanted quadrilateral incorrectly filled its exterior");
}

//! \brief Before initialise(), raster coordinates are not clipped to the
//! plotting area, and a later initialise() must reset those pixels.
void test_uninitialised_canvas_fill_and_reset() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    rectangle(canvas,30.0,20.0,90.0,55.0);
    canvas.write("default_graphics-before-initialise.png");
    const auto preinit=DrawingTest::read_default_graphics_png(
        "default_graphics-before-initialise.png");
    DrawingTest::require(
        preinit.pixel(60u,35u)==std::array<unsigned char,3>{255,0,0},
        "Uninitialised canvas failed to paint in device coordinates");
    DrawingTest::require(
        preinit.pixel(120u,35u)==std::array<unsigned char,3>{255,255,255},
        "Uninitialised canvas incorrectly painted outside its rectangle");
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    DefaultGraphicsCanvas fresh(ImageSize2d(128u,96u));
    fresh.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-preinit-reset.png");
    fresh.write("default_graphics-preinit-fresh.png");
    DrawingTest::require(
        DrawingTest::read_default_graphics_png(
            "default_graphics-preinit-reset.png").pixels==
        DrawingTest::read_default_graphics_png(
            "default_graphics-preinit-fresh.png").pixels,
        "Initialise did not clear pre-initialisation pixels");
}

//! \brief Horizontal-first and vertical-first slanted quadrilaterals
//! must be treated as general nonzero-winding contours.
void test_slanted_quadrilateral_edge_orientations() {
    const auto render=[](bool vertical_first,const char* filename) {
        DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
        canvas.initialise("","",-1.0,1.0,-1.0,1.0);
        canvas.set_line_width(0.0);
        canvas.set_fill_colour(1.0,0.0,0.0);
        canvas.move_to(-0.8,-0.6);
        if(vertical_first) {
            canvas.line_to(-0.8,0.6);
            canvas.line_to(0.5,0.8);
            canvas.line_to(0.5,-0.8);
        } else {
            canvas.line_to(0.5,-0.6);
            canvas.line_to(0.8,0.6);
            canvas.line_to(-0.5,0.6);
        }
        canvas.fill();
        canvas.write(filename);
        return DrawingTest::read_default_graphics_png(filename);
    };
    const auto horizontal=render(false,"default_graphics-slanted-hfirst.png");
    const auto vertical=render(true,"default_graphics-slanted-vfirst.png");
    for(auto const* image:{&horizontal,&vertical}) {
        DrawingTest::require(
            DrawingTest::world_pixel(*image,0.0,0.0)==
                std::array<unsigned char,3>{255,0,0},
            "Slanted quadrilateral failed to paint the interior");
        DrawingTest::require(
            DrawingTest::world_pixel(*image,0.9,0.9)==
                std::array<unsigned char,3>{247,247,247},
            "Slanted quadrilateral painted the exterior");
    }
}

//! \brief A nearly horizontal triangle with enormous finite x coordinates
//! must not generate an infinite-slope edge or paint unrepresented samples.
void test_extreme_thin_triangle_fill() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-thin-triangle-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-thin-triangle-baseline.png");
    canvas.set_line_width(0.0);
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.move_to(-1e306,0.0);
    canvas.line_to(1e306,1e-14);
    canvas.line_to(0.0,2e-14);
    canvas.fill();
    canvas.write("default_graphics-thin-triangle.png");
    const auto actual=DrawingTest::read_default_graphics_png(
        "default_graphics-thin-triangle.png");
    DrawingTest::require(actual.pixels==blank.pixels,
                         "Subpixel-thin triangle unexpectedly painted the raster");
}

//! \brief Repeating a consecutive vertex does not change a triangle fill.
void test_duplicate_consecutive_vertices_preserve_fill() {
    const auto render=[](int duplicate,const char* filename) {
        DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
        canvas.initialise("","",-1.0,1.0,-1.0,1.0);
        canvas.set_line_width(0.0);
        canvas.set_fill_colour(1.0,0.0,0.0);
        canvas.move_to(-0.7,-0.6);
        if(duplicate==1) canvas.line_to(-0.7,-0.6);
        canvas.line_to(0.7,-0.6);
        if(duplicate==2) canvas.line_to(0.7,-0.6);
        canvas.line_to(0.0,0.7);
        canvas.fill();
        canvas.write(filename);
        return DrawingTest::read_default_graphics_png(filename);
    };
    const auto triangle=render(0,"default_graphics-triangle-original.png");
    const auto repeat_first=render(1,"default_graphics-triangle-repeat-first.png");
    const auto repeat_second=render(2,"default_graphics-triangle-repeat-second.png");
    DrawingTest::require(triangle.pixels==repeat_first.pixels &&
                         triangle.pixels==repeat_second.pixels,
                         "A repeated consecutive vertex changed the filled triangle");
    DrawingTest::require(
        DrawingTest::world_pixel(triangle,0.0,0.0)==
            std::array<unsigned char,3>{255,0,0},
        "Reference triangle did not fill its interior");
}

//! \brief Subpixel-small device rectangles remain empty even if their
//! exact positive area underflows in floating-point cross products.
void test_subnormal_rectangle_cross_underflow() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.write("default_graphics-underflow-blank.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-underflow-blank.png");
    canvas.set_fill_colour(1.0,0.0,0.0);
    canvas.set_line_width(0.0);
    rectangle(canvas,0.0,0.0,1e-200,1e-200);
    canvas.write("default_graphics-underflow-rectangle.png");
    const auto result=DrawingTest::read_default_graphics_png(
        "default_graphics-underflow-rectangle.png");
    DrawingTest::require(result.pixels==blank.pixels,
                         "An unresolvable subpixel rectangle changed the image");
}

//! \brief Separate overflow in a segment's squared length from overflow
//! in its cross-product limit; neither invalid intermediate may paint.
void test_segment_arithmetic_overflow() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-segment-overflow-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-segment-overflow-baseline.png");
    canvas.set_line_colour(1.0,0.0,0.0);

    // Endpoints are finite, but squaring the device-space displacement
    // overflows. There is no representable segment length to rasterize.
    canvas.set_line_width(3.0);
    canvas.move_to(-1e154,-1e154);
    canvas.line_to(1e154,1e154);
    canvas.stroke();

    // Here length squared remains finite, but half-width times length
    // overflows; this is a separate early exit in the segment rasterizer.
    canvas.set_line_width(1e240);
    canvas.move_to(1e98,-0.001);
    canvas.line_to(2e98,0.001);
    canvas.stroke();

    canvas.write("default_graphics-segment-overflows.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-segment-overflows.png");
    DrawingTest::require(image.pixels==blank.pixels,
                         "Overflowed segment intermediates changed the image");
}

//! \brief Finite axial segments can overflow each rectangle bound while
//! remaining far outside the plot. They must be safely clipped.
void test_extreme_axial_strokes_on_all_sides() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-axial-extreme-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-axial-extreme-baseline.png");
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.set_line_width(1.78e308);
    const auto segment=[&](double x0,double y0,double x1,double y1) {
        canvas.move_to(x0,y0);
        canvas.line_to(x1,y1);
        canvas.stroke();
    };
    // A short axis-aligned segment keeps its cross-product limit finite,
    // but its width overflows xmin, xmax, ymin or ymax respectively.
    segment(-2.35e306,-0.01,-2.35e306,0.01);
    segment(2.35e306,-0.01,2.35e306,0.01);
    segment(-0.01,3.13e306,0.01,3.13e306);
    segment(-0.01,-3.13e306,0.01,-3.13e306);
    canvas.write("default_graphics-axial-extreme-four-sides.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-axial-extreme-four-sides.png");
    DrawingTest::require(image.pixels==blank.pixels,
                         "Extreme axial strokes painted outside their clip");
}

//! \brief Offscreen joins exercise triangular clipping in four directions.
void test_offscreen_miter_join_triangles() {
    DefaultGraphicsCanvas canvas(ImageSize2d(128u,96u));
    canvas.initialise("","",-1.0,1.0,-1.0,1.0);
    canvas.write("default_graphics-offscreen-joins-baseline.png");
    const auto blank=DrawingTest::read_default_graphics_png(
        "default_graphics-offscreen-joins-baseline.png");
    canvas.set_line_colour(1.0,0.0,0.0);
    canvas.set_line_width(6.0);
    const auto join=[&](double ax,double ay,double bx,double by,
                        double cx,double cy) {
        canvas.move_to(ax,ay);
        canvas.line_to(bx,by);
        canvas.line_to(cx,cy);
        canvas.stroke();
    };
    join(-4.0,0.0,-3.0,0.0,-3.0,0.7);
    join(3.0,0.0,4.0,0.0,4.0,0.7);
    join(0.0,3.0,0.7,3.0,0.7,4.0);
    join(0.0,-3.0,0.7,-3.0,0.7,-4.0);
    canvas.write("default_graphics-offscreen-miter-joins.png");
    const auto image=DrawingTest::read_default_graphics_png(
        "default_graphics-offscreen-miter-joins.png");
    DrawingTest::require(image.pixels==blank.pixels,
                         "Offscreen join triangles unexpectedly painted the canvas");
}

} // namespace

int main() {
    try {
        // Invalid dimensions must be rejected before any margin arithmetic
        // or potentially huge raster-buffer allocation.
        must_reject([] {
            DefaultGraphicsCanvas canvas(ImageSize2d(0u,8u));
        },"zero plot width accepted");
        must_reject([] {
            DefaultGraphicsCanvas canvas(ImageSize2d(8u,0u));
        },"zero plot height accepted");
        must_reject([] {
            DefaultGraphicsCanvas canvas(ImageSize2d(8u,16385u));
        },"excessive plot height accepted");
        must_fail_precondition([] {
            const ImageSize2d invalid(-1,8);
            (void)invalid;
        },"negative signed plot width accepted");
        must_fail_precondition([] {
            const ImageSize2d invalid(8,-1);
            (void)invalid;
        },"negative signed plot height accepted");
        // The signed overload must also accept valid dimensions and preserve
        // both coordinates, instead of merely rejecting invalid inputs.
        const ImageSize2d valid_signed_size(Int(13),Int(29));
        DrawingTest::require(
            valid_signed_size.nx==Nat(13) && valid_signed_size.ny==Nat(29),
            "Signed image dimensions were not preserved");
        must_reject([] {
            DefaultGraphicsCanvas canvas(ImageSize2d(
                std::numeric_limits<Nat>::max(),8u));
        },"huge plot width accepted");
        must_reject([] {
            DefaultGraphicsCanvas canvas(ImageSize2d(12000u,12000u));
        },"excessive raster area accepted");

        DefaultGraphicsCanvas canvas(ImageSize2d(64u,48u));
        const double nan_axes=std::numeric_limits<double>::quiet_NaN();
        const double infinite_axes=std::numeric_limits<double>::infinity();
        must_reject([&] {
            canvas.initialise("x","y",-1,nan_axes,-1,1);
        },"NaN upper x axis endpoint accepted");
        must_reject([&] {
            canvas.initialise("x","y",-1,1,-1,infinite_axes);
        },"Infinite upper y axis endpoint accepted");
        must_reject([&] {
            canvas.initialise("x","y",-1,1,1,1);
        },"Empty y interval accepted");

        const auto drawing_size=canvas.size_in_pixels();
        DrawingTest::require(drawing_size.nx==64u && drawing_size.ny==48u,
                             "Canvas reported incorrect drawable dimensions");
        const double nan=std::numeric_limits<double>::quiet_NaN();
        const double nonfinite_coordinate=std::numeric_limits<double>::infinity();
        const double tiny=std::numeric_limits<double>::denorm_min();
        must_reject([&] {
            canvas.initialise("x","y",0,0,-1,1);
        },"empty x interval accepted");
        must_reject([&] {
            canvas.initialise("x","y",1,-1,-1,1);
        },"reversed x interval accepted");
        must_reject([&] {
            canvas.initialise("x","y",nan,1,-1,1);
        },"NaN endpoint accepted");
        must_reject([&] {
            canvas.initialise("x","y",-1,1,-nonfinite_coordinate,1);
        },"infinite endpoint accepted");
        must_reject([&] {
            canvas.initialise("x","y",-1e308,1e308,-1,1);
        },"overflowing x interval accepted");
        must_reject([&] {
            canvas.initialise("x","y",-1,1,-1e308,1e308);
        },"overflowing y interval accepted");
        must_reject([&] {
            canvas.initialise("x","y",0,tiny,-1,1);
        },"non-finite x scale accepted");
        must_reject([&] {
            canvas.initialise("x","y",-1,1,0,tiny);
        },"non-finite y scale accepted");

        // The failed initializations must not prevent normal drawing.
        canvas.initialise("x","y",-1,1,-1,1);
        canvas.set_line_width(0.0); // Test pure fill and alpha only.
        canvas.set_fill_colour(1.0,0.0,0.0);
        canvas.set_fill_opacity(0.0);
        rectangle(canvas,-0.95,-0.65,-0.65,0.65); // No paint.
        canvas.set_fill_colour(0.0,1.0,0.0);
        canvas.set_fill_opacity(0.5);
        rectangle(canvas,-0.45,-0.65,-0.05,0.65);
        canvas.set_fill_colour(1.0,0.0,0.0);
        canvas.set_fill_opacity(1.0);
        rectangle(canvas,0.12,-0.65,0.55,0.65);
        canvas.set_fill_colour(0.0,0.0,1.0);
        rectangle(canvas,-4.0,-4.0,-0.8,-0.78); // Clip to plot.
        canvas.set_fill_colour(1.0,1.0,0.0);
        canvas.set_fill_opacity(2.0); // Clamp to opaque.
        rectangle(canvas,0.6,-0.65,0.9,0.65);

        // Non-finite path coordinates/radii are ignored. An invalid stroke
        // width cannot accidentally rasterize the remaining path.
        canvas.move_to(nan,0.0);
        canvas.line_to(nonfinite_coordinate,0.0);
        canvas.stroke();
        canvas.circle(0.0,0.0,nan);
        canvas.move_to(0.0,nonfinite_coordinate);
        canvas.line_to(0.0,nonfinite_coordinate);
        canvas.circle(0.0,nonfinite_coordinate,0.25);

        canvas.circle(nonfinite_coordinate,0.0,0.25);
        canvas.circle(0.0,0.0,-0.25);
        canvas.stroke();
        canvas.set_line_width(nan);
        canvas.move_to(-0.8,0.0);
        canvas.line_to(0.8,0.0);
        canvas.stroke();

        // Invalid PNG destinations must report failure.
        must_reject([&] { canvas.write(nullptr); },
                    "Null PNG filename accepted");
        bool rejected_missing_directory=false;
        try {
            canvas.write("default_graphics_missing_parent/invalid.png");
        } catch(const std::runtime_error&) {
            rejected_missing_directory=true;
        }
        DrawingTest::require(rejected_missing_directory,
                             "PNG writing succeeded in a nonexistent directory");

        test_embedded_axis_glyphs();
        test_canvas_reset_and_offscreen_geometry();
        test_nonrectangular_five_vertex_fill();
        test_clipped_rectangles_on_all_sides();
        test_overlapping_rectangular_contours();
        test_rectangle_fast_path_matches_generic_fill();
        test_degenerate_fills_leave_pixels_unchanged();
        test_paint_cache_survives_eviction();
        test_subpixel_empty_rectangle();
        test_overflowed_circle_geometry();
        test_nonfinite_paint_values();
        test_line_to_without_move_to();
        test_clipped_strokes_outside_plot();
        test_clipped_rectangle_matches_plot_bounds();
        test_acute_stroke_miter_limit();
        test_rectangle_recognition_rejects_degenerate_contours();
        test_mixed_nonrectangular_contours();
        test_same_rgb_different_opacity();
        test_empty_and_nonfinite_strokes();
        test_offscreen_generic_fill_spans();
        test_extreme_axial_stroke_clipping();
        test_long_axis_labels_clip_safely();
        test_zero_length_strokes_are_empty();
        test_slanted_quadrilateral_fill();
        test_uninitialised_canvas_fill_and_reset();
        test_slanted_quadrilateral_edge_orientations();
        test_extreme_thin_triangle_fill();
        test_duplicate_consecutive_vertices_preserve_fill();
        test_subnormal_rectangle_cross_underflow();
        test_segment_arithmetic_overflow();
        test_extreme_axial_strokes_on_all_sides();
        test_offscreen_miter_join_triangles();
        canvas.write("default_graphics-robustness.png");
        const auto image=DrawingTest::read_default_graphics_png(
            "default_graphics-robustness.png");
        DrawingTest::require(image.width==234 && image.height==98,
                             "Unexpected graphics image size");
        auto check=[&](double x,double y,std::array<unsigned char,3> expected,
                       const char* message) {
            DrawingTest::require(DrawingTest::world_pixel(image,x,y)==expected,
                                 message);
        };
        check(-0.73,0.0,{247,247,247},"Transparent fill changed background");
        check(-0.25,0.0,{124,251,124},"Half-opacity fill incorrect");
        check(0.30,0.0,{255,0,0},"Opaque fill incorrect");
        check(0.75,0.0,{255,255,0},"Clamped opacity incorrect");
        check(-0.90,-0.90,{0,0,255},"Clipping failed");
        check(0.0,0.0,{247,247,247},"Invalid geometry changed background");
        DrawingTest::require(image.pixel(230,25)==
                             std::array<unsigned char,3>{255,255,255},
                             "Clipping modified the margin");
        return 0;
    } catch(const std::exception& exception) {
        std::cerr<<"test_default_graphics_robustness: "<<exception.what()<<'\n';
        return 1;
    }
}
