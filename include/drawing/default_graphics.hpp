/***************************************************************************
 *            drawing/default_graphics.hpp
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

/*! \file drawing/default_graphics.hpp
 *  \brief Self-contained raster canvas and graphics backend.
 */

#ifndef ARIADNE_DEFAULT_GRAPHICS_HPP
#define ARIADNE_DEFAULT_GRAPHICS_HPP

#include <array>
#include <cstdint>
#include <vector>

#include "drawing/figure.hpp"
#include "drawing/image_size.hpp"

namespace Ariadne {

//! \ingroup GraphicsModule
//! \brief Software drawing canvas using 4x4 coverage masks and a built-in PNG writer.
//! The canvas owns all raster, coverage, path and paint state. It does not
//! require an external graphics, compression or font library.
class DefaultGraphicsCanvas : public CanvasBase {
  public:
    //! \brief Construct a canvas with the given drawing-area size in pixels.
    explicit DefaultGraphicsCanvas(const ImageSize2d& size);
    ~DefaultGraphicsCanvas() override = default;

    //! \brief Initialise the plot and its affine coordinate mapping.
    Void initialise(StringType x, StringType y, double xl, double xu,
                    double yl, double yu) override;

    //! \brief Write the current opaque raster to an RGBA PNG file.
    Void write(const char* filename) const override;

    //! \brief Draw the outer frame of the plot.
    Void finalise() override;

    //! \brief Begin a subpath at the specified world-space point.
    Void move_to(double x, double y) override;

    //! \brief Append a line segment to the current path.
    Void line_to(double x, double y) override;

    //! \brief Append a circle or ellipse, according to the coordinate transform.
    Void circle(double x, double y, double radius) override;

    //! \brief Append a dot using the configured radius.
    Void dot(double x, double y) override;

    //! \brief Stroke all current paths and consume them.
    Void stroke() override;

    //! \brief Fill the current paths with nonzero winding, then stroke and consume them.
    Void fill() override;

    //! \brief Set the radius used by subsequent dots.
    Void set_dot_radius(double radius);

    //! \brief Set the outline width in device pixels.
    Void set_line_width(double width) override;

    //! \brief Set the outline colour using normalized RGB channels.
    Void set_line_colour(double r, double g, double b) override;

    //! \brief Set the fill opacity, clamped at rasterization time.
    Void set_fill_opacity(double opacity) override;

    //! \brief Set the fill colour using normalized RGB channels.
    Void set_fill_colour(double r, double g, double b) override;

    //! \brief Return the drawing-area dimensions, excluding plot margins.
    ImageSize2d size_in_pixels() const;

  private:
    //! \brief Point in device coordinates.
    struct Point {
        double x, y;
    };

    //! \brief Accumulated path, optionally explicitly closed.
    struct Path {
        std::vector<Point> points;
        bool closed = false;
    };

    //! \brief Normalized RGB colour and independent opacity.
    struct Paint {
        double r, g, b, a;
    };

    Nat _width, _height;
    std::vector<std::uint8_t> _pixels;
    std::vector<Path> _paths;

    // Per-pixel subpixel coverage is cleared after each paint operation.
    std::vector<std::uint16_t> _coverage;
    std::vector<std::size_t> _dirty_pixels;

    // World-to-device affine transform and current drawing properties.
    double _sx = 1.0, _sy = 1.0, _tx = 0.0, _ty = 0.0;
    double _line_width = 1.0, _dot_radius = 1.0;
    Paint _line{0.0, 0.0, 0.0, 1.0};
    Paint _fill{1.0, 1.0, 1.0, 1.0};
    bool _clipped = false;

    // Four-entry paint cache avoids regenerating 13,056 channel lookup values
    // whenever fill and stroke colours alternate on successive primitives.
    using BlendTable = std::array<std::array<std::array<std::uint8_t, 256>, 3>, 17>;

    struct BlendCacheEntry {
        BlendTable table{};
        Paint paint{0.0, 0.0, 0.0, 0.0};
        std::uint64_t last_used = 0;
        bool valid = false;
    };

    std::array<BlendCacheEntry, 4> _blend_cache{};
    std::uint64_t _blend_cache_clock = 0;

    //! \brief Apply the current world-to-device mapping.
    Point transform(double x, double y) const;

    //! \brief Blend a partially covered pixel in the opaque raster.
    void blend(int x, int y, Paint paint, double coverage);

    //! \brief Paint an integer bounding box, used for backgrounds.
    void rectangle(int x0, int y0, int x1, int y1, Paint paint);

    //! \brief Mark one of sixteen samples covered by the current operation.
    void mark_sample(int x, int y, int sx, int sy);

    //! \brief Retrieve or build the lookup table for one paint colour.
    BlendTable const& blend_table_for(Paint paint);

    //! \brief Blend the union coverage once and reset the dirty masks.
    void composite_coverage(Paint paint);

    //! \brief Rasterize a horizontal interval at a subpixel row.
    void rasterise_span(double from, double to, int y, int sy, bool inclusive_end);

    //! \brief Rasterize one butt-capped stroke segment.
    void rasterise_segment(Point a, Point b, double width);

    //! \brief Rasterize a circular glyph endpoint or dot.
    void rasterise_disc(Point center, double radius);

    //! \brief Rasterize a join triangle.
    void rasterise_triangle(Point a, Point b, Point c);

    //! \brief Rasterize a join, including its bounded miter.
    void rasterise_join(Point a, Point b, Point c, double half_width);

    //! \brief Fill the accumulated compound path using nonzero winding.
    void rasterise_fill();

    //! \brief Rasterize small built-in glyphs for axis labels.
    void text(double x, double y, StringType const& label);
};

//! \ingroup GraphicsModule
//! \brief Create self-contained raster canvases through the graphics abstraction.
class DefaultGraphicsBackend : public GraphicsBackendInterface {
  public:
    //! \brief Create a canvas for the requested drawing-area size.
    SharedPointer<CanvasInterface> make_canvas(const char* filename, Nat width,
                                               Nat height, Bool animated) const override;
};

} // namespace Ariadne

#endif // ARIADNE_DEFAULT_GRAPHICS_HPP
