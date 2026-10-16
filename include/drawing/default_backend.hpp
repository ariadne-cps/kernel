/*
 * Ariadne - free software under the GNU General Public License version 3
 * or (at your option) any later version. No warranty.
 */
#ifndef ARIADNE_DEFAULT_BACKEND_HPP
#define ARIADNE_EMBEDDED_HPP

#include <array>
#include <cstdint>
#include <vector>

#include "drawing/figure.hpp"
#include "drawing/image_size.hpp"

namespace Ariadne {

// Default graphics backend: self-contained software rasterizer.
class DefaultGraphicsCanvas : public CanvasBase {
  public:
    explicit DefaultGraphicsCanvas(const ImageSize2d& size);
    ~DefaultGraphicsCanvas() override = default;

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

    // Four-entry paint cache avoids regenerating 13,056 channel lookup values
    // whenever fill and stroke colours alternate on successive primitives.
    using BlendTable = std::array<std::array<std::array<std::uint8_t,256>,3>,17>;
    struct BlendCacheEntry {
        BlendTable table{};
        Paint paint{0.0,0.0,0.0,0.0};
        std::uint64_t last_used=0;
        bool valid=false;
    };
    std::array<BlendCacheEntry,4> _blend_cache{};
    std::uint64_t _blend_cache_clock=0;

    Point transform(double x, double y) const;
    void blend(int x, int y, Paint paint, double coverage);
    void rectangle(int x0, int y0, int x1, int y1, Paint paint);
    void mark_sample(int x, int y, int sx, int sy);
    BlendTable const& blend_table_for(Paint paint);
    void composite_coverage(Paint paint);
    void rasterise_span(double from, double to, int y, int sy, bool inclusive_end);
    void rasterise_segment(Point a, Point b, double width);
    void rasterise_disc(Point center, double radius);
    void rasterise_triangle(Point a, Point b, Point c);
    void rasterise_join(Point a, Point b, Point c, double half_width);
    void rasterise_fill();
    void text(double x, double y, StringType const& label);
};

class DefaultGraphicsBackend : public GraphicsBackendInterface {
  public:
    SharedPointer<CanvasInterface> make_canvas(const char* filename, Nat width,
                                               Nat height, Bool animated) const override;
};

} // namespace Ariadne
#endif // ARIADNE_EMBEDDED_HPP
