/*
 * Ariadne - free software under the GNU General Public License version 3
 * or (at your option) any later version. No warranty.
 */
#ifndef ARIADNE_EMBEDDED_HPP
#define ARIADNE_EMBEDDED_HPP

#include <cstdint>
#include <vector>

#include "drawing/figure.hpp"
#include "drawing/image_size.hpp"

namespace Ariadne {

// Software rasterizer: no Cairo, external image or compression library required.
class EmbeddedCanvas : public CanvasBase {
  public:
    explicit EmbeddedCanvas(const ImageSize2d& size);
    ~EmbeddedCanvas() override = default;

    Void initialise(StringType x, StringType y, double xl, double xu,
                    double yl, double yu) override;
    Void write(const char* filename) const override;
    Void finalise() override;

    Void move_to(double x, double y) override;
    Void line_to(double x, double y) override;
    Void circle(double x, double y, double radius) override;
    Void dot(double x, double y) override;
    Void stroke() override;
    Void fill() override;

    Void set_dot_radius(double radius);
    Void set_line_width(double width) override;
    Void set_line_colour(double r, double g, double b) override;
    Void set_fill_opacity(double opacity) override;
    Void set_fill_colour(double r, double g, double b) override;
    ImageSize2d size_in_pixels() const;

  private:
    struct Point { double x, y; };
    struct Path { std::vector<Point> points; bool closed = false; };
    struct Paint { double r, g, b, a; };

    Nat _width, _height;
    std::vector<std::uint8_t> _pixels;
    std::vector<Path> _paths;
    std::vector<std::uint16_t> _coverage;
    std::vector<std::size_t> _dirty_pixels;
    double _sx = 1.0, _sy = 1.0, _tx = 0.0, _ty = 0.0;
    double _line_width = 1.0, _dot_radius = 1.0;
    Paint _line{0.0, 0.0, 0.0, 1.0};
    Paint _fill{1.0, 1.0, 1.0, 1.0};
    bool _clipped = false;

    Point transform(double x, double y) const;
    void blend(int x, int y, Paint paint, double coverage);
    void rectangle(int x0, int y0, int x1, int y1, Paint paint);
    void mark_sample(int x, int y, int sx, int sy);
    void composite_coverage(Paint paint);
    void rasterise_segment(Point a, Point b, double width);
    void rasterise_triangle(Point a, Point b, Point c);
    void rasterise_join(Point a, Point b, Point c, double half_width);
    void rasterise_fill();
    void text(double x, double y, StringType const& label);
};

class EmbeddedGraphicsBackend : public GraphicsBackendInterface {
  public:
    SharedPointer<CanvasInterface> make_canvas(const char* filename, Nat width,
                                               Nat height, Bool animated) const override;
};

} // namespace Ariadne
#endif // ARIADNE_EMBEDDED_HPP
