#include "drawing/default_backend.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include "drawing/cairo.hpp"
#endif
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <stdexcept>
#include "default_backend_png.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include <cairo.h>
#endif

using namespace Ariadne;

static void render(CanvasInterface& canvas, const char* filename) {
    canvas.initialise("x", "y", -1.0, 1.0, -1.0, 1.0);
    canvas.set_line_width(2.0);
    canvas.set_line_colour(0.1, 0.2, 0.7);
    canvas.set_fill_colour(0.8, 0.3, 0.1);
    canvas.set_fill_opacity(0.6);
    canvas.move_to(-0.6, -0.5);
    canvas.line_to(0.7, -0.5);
    canvas.line_to(0.0, 0.8);
    canvas.line_to(-0.6, -0.5);
    canvas.fill();
    canvas.set_line_colour(0.1, 0.6, 0.2);
    canvas.circle(0.0, 0.0, 0.35);
    canvas.stroke();
    canvas.write(filename);
}

// Geometric regression cases: shared path winding, joins, clipping, arcs.
static void geometry(CanvasInterface& c, const char* filename, int scene) {
    c.initialise("x", "y", -1.0, 1.0, -1.0, 1.0);
    c.set_line_width(3.0);
    c.set_line_colour(0.05, 0.1, 0.4);
    c.set_fill_colour(0.75, 0.25, 0.1);
    c.set_fill_opacity(0.6);
    if(scene == 0) {
        // Outer CCW and inner CW: nonzero winding leaves a hole.
        c.move_to(-0.8,-0.8); c.line_to(0.8,-0.8);
        c.line_to(0.8,0.8); c.line_to(-0.8,0.8);
        c.line_to(-0.8,-0.8);
        c.move_to(-0.4,-0.4); c.line_to(-0.4,0.4);
        c.line_to(0.4,0.4); c.line_to(0.4,-0.4);
        c.line_to(-0.4,-0.4);
        c.fill();
    } else if(scene == 1) {
        c.set_line_width(18.0);
        c.move_to(-0.85,-0.65);
        c.line_to(0.0,0.65);
        c.line_to(0.85,-0.65);
        c.stroke();
        c.set_fill_opacity(0.3);
        c.set_line_width(3.0);
        c.move_to(-0.3,0.0); c.line_to(0.3,0.0);
        c.line_to(0.3,-0.5); c.line_to(-0.3,-0.5);
        c.fill();
    } else if(scene == 2) {
        c.set_line_width(5.0);
        c.move_to(-1.6,-1.6);
        c.line_to(1.6,1.6);
        c.stroke();
        c.circle(0,0,0.6);
        c.fill();
    } else if(scene == 3) {
        c.set_line_width(4.0);
        c.move_to(-0.8,0.0);
        c.circle(0.0,0.0,0.35);
        c.line_to(0.8,0.0);
        c.stroke();
    } else if(scene == 4) {
        // Both explicit and implicit closing edges, at fractional device
        // coordinates. Also check winding direction and clipping at the edge.
        c.set_line_width(0.8);
        c.set_fill_opacity(0.4);
        for(int row=0;row<7;++row) for(int col=0;col<10;++col) {
            double x=-1.12+0.238*col;
            double y=-1.08+0.336*row;
            double w=0.192,h=0.271;
            if((row+col)%2==0) {
                c.move_to(x,y);
                c.line_to(x+w,y);
                c.line_to(x+w,y+h);
                c.line_to(x,y+h);
                c.line_to(x,y);  // Explicit closure.
            } else {
                c.move_to(x,y);
                c.line_to(x,y+h);
                c.line_to(x+w,y+h);
                c.line_to(x+w,y); // Implicit closure on fill.
            }
            c.fill();
        }
    } else if(scene == 5) {
        // Butt caps, subpixel-width strokes and clipped segments at
        // multiple orientations, including exact horizontal/vertical lines.
        for(int k=0;k<12;++k) {
            c.set_line_width(0.75+0.28*(k%5));
            double y=-0.90+0.16*k;
            c.move_to(-1.2,y);
            c.line_to(1.2,y+0.12*std::sin(k));
            c.stroke();
        }
        c.set_line_width(3.25);
        c.move_to(0.0,-1.2); c.line_to(0.0,1.2); c.stroke();
        c.move_to(-1.2,0.0); c.line_to(1.2,0.0); c.stroke();
    } else if(scene == 6) {
        // Concave multi-edge fill must NOT take the rectangle fast path.
        c.move_to(-0.8,-0.8); c.line_to(0.8,-0.8);
        c.line_to(0.2,-0.2); c.line_to(0.8,0.8);
        c.line_to(-0.8,0.8); c.line_to(-0.2,0.1);
        c.fill();
    } else if(scene == 7 || scene == 8) {
        // Nested rectangles: reversed subpath order, same/opposite winding,
        // fractional coordinates and clipping. Also exercise the generic
        // path when nested contours touch or otherwise fail the fast path.
        c.set_line_width(1.25);
        c.set_fill_opacity(0.55);
        auto rectangle=[&c](double xl,double yl,double xu,double yu,
                             bool reverse) {
            if(reverse) {
                c.move_to(xl,yl); c.line_to(xl,yu);
                c.line_to(xu,yu); c.line_to(xu,yl);
                c.line_to(xl,yl);
            } else {
                c.move_to(xl,yl); c.line_to(xu,yl);
                c.line_to(xu,yu); c.line_to(xl,yu);
                c.line_to(xl,yl);
            }
        };
        if(scene == 7) {
            // Opposite winding; both possible contour submission orders.
            rectangle(-0.95,-0.69,-0.38,0.68,false);
            rectangle(-0.81,-0.42,-0.52,0.39,true);
            c.fill();
            rectangle(-0.19,-0.41,0.13,0.42,true);
            rectangle(-0.33,-0.69,0.27,0.68,false);
            c.fill();
            rectangle(0.43,-0.38,0.81,0.41,false);
            rectangle(0.36,-0.69,0.96,0.68,true);
            c.fill();
        } else {
            // Same winding fills the inner region too. Nested contours
            // touching the outer boundary must instead use general filling.
            rectangle(-1.08,-0.68,-0.17,0.65,false);
            rectangle(-0.83,-0.41,-0.37,0.43,false);
            c.fill();
            rectangle(0.17,-0.68,1.08,0.65,false);
            rectangle(0.44,-0.41,0.83,0.43,true);
            c.fill();
            rectangle(-0.52,-0.90,0.52,-0.75,false);
            rectangle(-0.30,-0.90,0.30,-0.78,true);
            c.fill();
        }
    } else if(scene == 9) {
        // Fractional horizontal and vertical strokes: short segments, caps,
        // thick lines and lines clipped by each edge of the plotting area.
        for(int index=0;index<10;++index) {
            const double offset=-0.94+index*0.21;
            c.set_line_width(0.25+index*0.43);
            c.move_to(-1.18,offset);
            c.line_to(1.18,offset);
            c.stroke();
            c.move_to(offset,-1.20);
            c.line_to(offset,1.20);
            c.stroke();
        }
        c.set_line_width(4.0);
        c.move_to(-0.95,-0.96);
        c.line_to(-0.37,-0.96);
        c.stroke();
        c.move_to(0.77,0.94);
        c.line_to(0.77,0.78);
        c.stroke();
    }
    c.write(filename);
}


static void verify_labels(DrawingTest::Image const& image) {
    using DrawingTest::require;
    auto dark=[&](int x0,int y0,int x1,int y1) {
        unsigned count=0;
        for(int y=y0;y<y1;++y) for(int x=x0;x<x1;++x) {
            const auto colour=image.pixel(static_cast<unsigned>(x),
                                           static_cast<unsigned>(y));
            if(colour[0]<90 && colour[1]<90 && colour[2]<90) ++count;
        }
        return count;
    };
    const int width=static_cast<int>(image.width);
    const int height=static_cast<int>(image.height);
    const int plot_height=height-50;
    const int x_center=160+(width-170)/2;
    require(dark(x_center-25,height-37,x_center+15,height-2)>=45,
            "Horizontal axis label is not legible");
    require(dark(125,10+plot_height/2-4,159,10+plot_height/2+36)>=65,
            "Vertical axis label is not legible");
    require(dark(127,8,158,45)>=45,"Upper numeric axis label is not legible");
}

static void verify_rectangles(DrawingTest::Image const& holes,
                               DrawingTest::Image const& mixed) {
    using DrawingTest::require;
    auto expect=[](DrawingTest::Image const& image,double x,double y,bool filled) {
        const auto colour=DrawingTest::world_pixel(image,x,y);
        const bool actual=colour[1]<210 && colour[2]<210;
        DrawingTest::require(actual==filled,"Rectangular winding mismatch");
    };
    require(holes.width==490 && holes.height==290,"Scene dimensions changed");
    require(mixed.width==490 && mixed.height==290,"Scene dimensions changed");
    for(const double x:{-0.66,-0.03,0.62}) expect(holes,x,0,false);
    expect(holes,-0.91,0,true);
    expect(holes,-0.28,0,true);
    expect(mixed,-0.60,0,true);
    expect(mixed,0.64,0,false);
    expect(mixed,0.94,0,true);
    expect(mixed,0,-0.82,false);
    expect(mixed,-0.45,-0.82,true);
}

#ifdef ARIADNE_HAS_CAIRO
static void compare_optional_reference(std::string const& name,
                                       DrawingTest::Image const& embedded) {
    cairo_surface_t* surface=cairo_image_surface_create_from_png(name.c_str());
    if(cairo_surface_status(surface)!=CAIRO_STATUS_SUCCESS) {
        cairo_surface_destroy(surface);
        throw std::runtime_error("Cannot read reference PNG: "+name);
    }
    cairo_surface_flush(surface);
    DrawingTest::require(cairo_image_surface_get_width(surface)==
                         static_cast<int>(embedded.width) &&
                         cairo_image_surface_get_height(surface)==
                         static_cast<int>(embedded.height),
                         "Reference image dimensions differ");
    const int stride=cairo_image_surface_get_stride(surface);
    const auto* data=cairo_image_surface_get_data(surface);
    double error=0;
    std::size_t outliers=0;
    std::size_t pixels=0;
    for(unsigned y=10;y<embedded.height-40;++y)
        for(unsigned x=160;x<embedded.width-10;++x) {
            const auto* ref=data+static_cast<std::size_t>(y)*stride+
                            static_cast<std::size_t>(x)*4;
            const auto colour=embedded.pixel(x,y);
            const std::array<unsigned char,3> reference{ref[2],ref[1],ref[0]};
            int max_delta=0;
            for(int k=0;k<3;++k) {
                const int delta=std::abs(int(colour[static_cast<std::size_t>(k)])-
                                         int(reference[static_cast<std::size_t>(k)]));
                error+=delta;
                max_delta=std::max(max_delta,delta);
            }
            if(max_delta>32) ++outliers;
            ++pixels;
        }
    cairo_surface_destroy(surface);
    const double mean=error/(3.0*255.0*static_cast<double>(pixels));
    const double fraction=static_cast<double>(outliers)/static_cast<double>(pixels);
    DrawingTest::require(mean<=0.015 && fraction<=0.035,
                         "Reference backend pixel comparison failed");
}
#endif

int main() {
    try {
    DefaultGraphicsCanvas default_backend(ImageSize2d(320u, 240u));
    render(default_backend, "comparison-default_backend.png");
    verify_labels(DrawingTest::read_default_backend_png("comparison-default_backend.png"));
#ifdef ARIADNE_HAS_CAIRO
    CairoCanvas cairo(ImageSize2d(320u, 240u));
    render(cairo, "comparison-cairo.png");
#endif
    for(int scene=0;scene<10;++scene) {
        DefaultGraphicsCanvas e(ImageSize2d(320u,240u));
        std::string ename="comparison-scene"+std::to_string(scene)+"-default_backend.png";
        geometry(e,ename.c_str(),scene);
        const auto generated=DrawingTest::read_default_backend_png(ename);

#ifdef ARIADNE_HAS_CAIRO
        CairoCanvas c(ImageSize2d(320u,240u));
        std::string cname="comparison-scene"+std::to_string(scene)+"-cairo.png";
        geometry(c,cname.c_str(),scene);
        compare_optional_reference(cname,generated);
#endif
    }
    verify_rectangles(
        DrawingTest::read_default_backend_png("comparison-scene7-default_backend.png"),
        DrawingTest::read_default_backend_png("comparison-scene8-default_backend.png"));
    return 0;
    } catch(std::exception const& error) {
        std::cerr<<"default graphics backend comparison: "<<error.what()<<'\\n';
        return 1;
    }
}
