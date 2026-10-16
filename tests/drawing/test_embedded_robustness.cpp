/*
 * Ariadne - free software under the GNU General Public License version 3
 * or (at your option) any later version. No warranty.
 */
#include "drawing/embedded.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace Ariadne;

namespace {

template<class Function>
void must_reject(Function function,const char* scenario) {
    try {
        function();
    } catch(const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(scenario);
}

void rectangle(EmbeddedCanvas& canvas,double xl,double yl,
               double xu,double yu) {
    canvas.move_to(xl,yl);
    canvas.line_to(xu,yl);
    canvas.line_to(xu,yu);
    canvas.line_to(xl,yu);
    canvas.fill();
}

} // namespace

int main() {
    try {
        // Invalid dimensions must be rejected before any margin arithmetic
        // or potentially huge raster-buffer allocation.
        must_reject([] {
            EmbeddedCanvas canvas(ImageSize2d(0u,8u));
        },"zero plot width accepted");
        must_reject([] {
            EmbeddedCanvas canvas(ImageSize2d(8u,0u));
        },"zero plot height accepted");
        must_reject([] {
            EmbeddedCanvas canvas(ImageSize2d(
                std::numeric_limits<Nat>::max(),8u));
        },"huge plot width accepted");
        must_reject([] {
            EmbeddedCanvas canvas(ImageSize2d(12000u,12000u));
        },"excessive raster area accepted");

        EmbeddedCanvas canvas(ImageSize2d(64u,48u));
        const double nan=std::numeric_limits<double>::quiet_NaN();
        const double inf=std::numeric_limits<double>::infinity();
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
            canvas.initialise("x","y",-1,1,-inf,1);
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
        canvas.line_to(inf,0.0);
        canvas.stroke();
        canvas.circle(0.0,0.0,nan);
        canvas.circle(inf,0.0,0.25);
        canvas.circle(0.0,0.0,-0.25);
        canvas.stroke();
        canvas.set_line_width(nan);
        canvas.move_to(-0.8,0.0);
        canvas.line_to(0.8,0.0);
        canvas.stroke();

        canvas.write("embedded-robustness.png");
        return 0;
    } catch(const std::exception& exception) {
        std::cerr<<"test_embedded_robustness: "<<exception.what()<<'\n';
        return 1;
    }
}
