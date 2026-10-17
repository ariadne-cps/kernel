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
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.0,0.0)==
            std::array<unsigned char,3>{255,0,0},
        "Near-reversing join failed to paint its stroke vertex");
    DrawingTest::require(
        DrawingTest::world_pixel(image,0.0,0.8)==
            std::array<unsigned char,3>{247,247,247},
        "Excessive miter spike extended beyond the stroke");
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
