#include "drawing/embedded.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include "drawing/cairo.hpp"
#endif
#include <cmath>
#include <memory>
#include <string>

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
    }
    c.write(filename);
}

int main() {
    EmbeddedCanvas embedded(ImageSize2d(320u, 240u));
    render(embedded, "comparison-embedded.png");
#ifdef ARIADNE_HAS_CAIRO
    CairoCanvas cairo(ImageSize2d(320u, 240u));
    render(cairo, "comparison-cairo.png");
#endif
    for(int scene=0;scene<7;++scene) {
        EmbeddedCanvas e(ImageSize2d(320u,240u));
        std::string ename="comparison-scene"+std::to_string(scene)+"-embedded.png";
        geometry(e,ename.c_str(),scene);
#ifdef ARIADNE_HAS_CAIRO
        CairoCanvas c(ImageSize2d(320u,240u));
        std::string cname="comparison-scene"+std::to_string(scene)+"-cairo.png";
        geometry(c,cname.c_str(),scene);
#endif
    }
    return 0;
}
