#include <iostream>
#include <string>

#include "ariadne-kernel.hpp"

using namespace Ariadne;
using std::cout;

namespace {

void section(const std::string& title) {
    cout << "\n== " << title << " ==\n";
}

template<class T>
void show(const std::string& label, const T& value) {
    cout << label << ": " << value << "\n";
}

ValidatedConstrainedImageSet make_demo_set() {
    auto u = EffectiveScalarMultivariateFunction::coordinates(3);
    EffectiveVectorMultivariateFunction map({
        1.25_dec * u[0] + 0.35_dec * u[2],
        u[1] + 0.20_dec * u[0] * u[0],
        u[2] - 0.30_dec * u[0]
    });
    return ValidatedConstrainedImageSet(ExactBoxType::unit_box(3), map);
}

} // namespace

int main() {
    section("Graphics backend and drawer policy");

    // Figure is intentionally separated from the rendering backend.  The
    // GraphicsManager owns the process-wide backend and the strategy used to
    // draw nonlinear validated sets.
    GraphicsManager& manager = GraphicsManager::instance();
    manager.set_backend(DefaultGraphicsBackend());
    manager.set_drawer(AffineDrawer(1));

    show("current drawer", manager.drawer());

    section("Colours and reusable graphics properties");

    Colour custom_blue(0.15, 0.35, 0.80, 0.60);
    GraphicsProperties style;
    style.set_line_width(2.0)
         .set_line_colour(black)
         .set_fill_colour(custom_blue)
         .set_fill_opacity(0.60);

    show("custom colour", custom_blue);
    show("graphics properties", style);

    section("Layering a two-dimensional figure");

    // A Figure records both the shape and the current drawing properties.
    // Changing a property affects only objects drawn afterwards, which makes
    // it natural to build a scene layer by layer.
    ApproximateBoxType viewport({
        {-2.5, 2.5},
        {-2.0, 2.0}
    });
    Figure figure(viewport, Projection2d(2, 0, 1));

    RealBox region_a({
        {-1.8_dec, 0.2_dec},
        {-1.0_dec, 0.7_dec}
    });
    RealBox region_b({
        {-0.4_dec, 1.7_dec},
        {-0.2_dec, 1.4_dec}
    });

    figure.set_fill_colour(0.20, 0.60, 0.90)
          .set_fill_opacity(0.45)
          .set_line_colour(0.05, 0.20, 0.35)
          .set_line_width(2.0)
          .draw(region_a);

    figure.set_fill_colour(0.95, 0.55, 0.15)
          .set_fill_opacity(0.35)
          .set_line_colour(0.45, 0.20, 0.05)
          .draw(region_b);

    show("figure bounding box", figure.get_bounding_box());
    show("line width", figure.get_line_width());
    show("fill opacity", figure.get_fill_opacity());

    figure.write("tutorial_drawing_layers");
    figure.write("tutorial_drawing_layers_large", 1000, 700);

    section("Reusing a figure with clear()");

    // clear() removes the queued graphical objects but keeps the viewport,
    // projection and current properties.  It is useful when producing a
    // family of related plots without reconstructing the Figure.
    figure.clear();
    figure.set_fill_colour(green)
          .set_fill_opacity(0.30)
          .set_line_colour(darkgrey)
          .set_line_width(1.0)
          .draw(intersection(region_a, region_b));
    figure.write("tutorial_drawing_intersection");

    section("Projecting higher-dimensional sets");

    ValidatedConstrainedImageSet set3d = make_demo_set();
    ApproximateBoxType viewport3d({
        {-2.0, 2.0},
        {-2.0, 2.0},
        {-2.0, 2.0}
    });

    // The same three-dimensional object can be rendered through different
    // Projection2d maps.  The set itself is not modified.
    Figure xy(viewport3d, Projection2d(3, 0, 1));
    xy.set_fill_colour(cyan).set_fill_opacity(0.40).draw(set3d);
    xy.write("tutorial_drawing_projection_xy");

    Figure xz(viewport3d, Projection2d(3, 0, 2));
    xz.set_fill_colour(magenta).set_fill_opacity(0.40).draw(set3d);
    xz.write("tutorial_drawing_projection_xz");

    show("xy projection", xy.get_2d_projection_map());
    show("xz projection", xz.get_2d_projection_map());

    section("Choosing how validated nonlinear sets are drawn");

    // The drawer is independent from Figure.  BoxDrawer is cheapest and
    // coarsest, AffineDrawer preserves affine correlation and can recursively
    // split the parameter domain, while GridDrawer renders a cell paving.
    ApproximateBoxType set_view({
        {-2.0, 2.0},
        {-2.0, 2.0},
        {-2.0, 2.0}
    });

    manager.set_drawer(BoxDrawer());
    Figure box_view(set_view, Projection2d(3, 0, 1));
    box_view.set_fill_colour(red).set_fill_opacity(0.25).draw(set3d);
    box_view.write("tutorial_drawing_drawer_box");

    manager.set_drawer(AffineDrawer(2));
    Figure affine_view(set_view, Projection2d(3, 0, 1));
    affine_view.set_fill_colour(blue).set_fill_opacity(0.25).draw(set3d);
    affine_view.write("tutorial_drawing_drawer_affine");

    manager.set_drawer(GridDrawer(4));
    Figure grid_view(set_view, Projection2d(3, 0, 1));
    grid_view.set_fill_colour(orange).set_fill_opacity(0.25).draw(set3d);
    grid_view.write("tutorial_drawing_drawer_grid");

    // Restore the library default used by normal applications.
    manager.set_drawer(AffineDrawer(0));

    section("plot(): a one-shot convenience API");

    // plot() is useful when no later mutation of the Figure is needed.
    ApproximateBoxType quick_view({
        {-2.0, 2.0},
        {-2.0, 2.0}
    });
    plot("tutorial_drawing_plot", quick_view,
         blue, region_a,
         orange, region_b);

    section("What the Drawing layer is responsible for");
    cout
        << "Figure: collect drawable objects and their per-layer style.\n"
        << "Projection2d: choose which coordinates of a higher-dimensional object are displayed.\n"
        << "GraphicsProperties / Colour: control line, fill and opacity independently of geometry.\n"
        << "GraphicsManager: choose the rendering backend and nonlinear-set drawer globally.\n"
        << "BoxDrawer / AffineDrawer / GridDrawer: trade rendering cost for geometric detail.\n"
        << "DefaultGraphicsBackend: materialise figures as PNG files.\n";

    return 0;
}
