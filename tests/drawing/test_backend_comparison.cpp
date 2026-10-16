#include "drawing/embedded.hpp"
#ifdef ARIADNE_HAS_CAIRO
#include "drawing/cairo.hpp"
#endif
#include <memory>

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

int main() {
    EmbeddedCanvas embedded(ImageSize2d(320u, 240u));
    render(embedded, "comparison-embedded.png");
#ifdef ARIADNE_HAS_CAIRO
    CairoCanvas cairo(ImageSize2d(320u, 240u));
    render(cairo, "comparison-cairo.png");
#endif
    return 0;
}
