/***************************************************************************
 *            src/drawing/default_graphics.cpp
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

/*! \file src/drawing/default_graphics.cpp
 *  \brief Default graphics rasterization, compositing, labels and PNG output.
 */

#include "drawing/default_graphics.hpp"
#include "utility/stlio.hpp"
#include "utility/standard.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Ariadne {

namespace {
constexpr int LEFT = 160;
constexpr int RIGHT = 10;
constexpr int TOP = 10;
constexpr int BOTTOM = 40;
constexpr double PI = 3.1415926535897932384626433832795;

double unit(double value) {
    return std::isfinite(value) ? std::clamp(value, 0.0, 1.0) : 0.0;
}

std::uint8_t channel(double value) {
    return static_cast<std::uint8_t>(std::lround(255.0 * unit(value)));
}

// Blend in 8-bit channel space with one explicitly rounded operation.
// The normalized expression (old/255)*(1-alpha)+source*alpha can land
// just below an exact half-integer on some targets (e.g. 150.49999999999997
// instead of 150.5), producing different PNG bytes across toolchains.
std::uint8_t composite_channel(std::uint8_t previous, double source_byte, double alpha) {
    const double old_byte=static_cast<double>(previous);
    const double result=std::fma(alpha,source_byte-old_byte,old_byte);
    return static_cast<std::uint8_t>(std::lround(result));
}

// A minimal, built-in 5x7 font. Rows are bit masks, high bit at the left.
std::array<std::uint8_t,7> glyph(char ch) {
    switch (ch) {
        case '0': return {14,17,19,21,25,17,14};
        case '1': return {4,12,4,4,4,4,14};
        case '2': return {14,17,1,2,4,8,31};
        case '3': return {30,1,1,14,1,1,30};
        case '4': return {2,6,10,18,31,2,2};
        case '5': return {31,16,30,1,1,17,14};
        case '6': return {6,8,16,30,17,17,14};
        case '7': return {31,1,2,4,8,8,8};
        case '8': return {14,17,17,14,17,17,14};
        case '9': return {14,17,17,15,1,2,12};
        case '-': return {0,0,0,31,0,0,0};
        case '+': return {0,4,4,31,4,4,0};
        case '.': return {0,0,0,0,0,12,12};
        case ',': return {0,0,0,0,6,6,4};
        case ':': return {0,12,12,0,12,12,0};
        case '_': return {0,0,0,0,0,0,31};
        case '/': return {1,2,2,4,8,8,16};
        case '(': return {2,4,8,8,8,4,2};
        case ')': return {8,4,2,2,2,4,8};
        case ' ': return {};
        default: break;
    }
    // Lower and upper case have the same small bitmap shapes.
    switch (static_cast<char>((ch >= 'A' && ch <= 'Z') ? ch+('a'-'A') : ch)) {
        case 'a': return {0,14,1,15,17,19,13};
        case 'b': return {16,16,30,17,17,17,30};
        case 'c': return {0,14,17,16,16,17,14};
        case 'd': return {1,1,15,17,17,17,15};
        case 'e': return {0,14,17,31,16,17,14};
        case 'f': return {6,9,8,28,8,8,8};
        case 'g': return {0,15,17,17,15,1,14};
        case 'h': return {16,16,30,17,17,17,17};
        case 'i': return {4,0,12,4,4,4,14};
        case 'j': return {2,0,6,2,2,18,12};
        case 'k': return {16,18,20,24,20,18,17};
        case 'l': return {12,4,4,4,4,4,14};
        case 'm': return {0,26,21,21,21,21,21};
        case 'n': return {0,30,17,17,17,17,17};
        case 'o': return {0,14,17,17,17,17,14};
        case 'p': return {0,30,17,17,30,16,16};
        case 'q': return {0,15,17,17,15,1,1};
        case 'r': return {0,22,25,16,16,16,16};
        case 's': return {0,15,16,14,1,1,30};
        case 't': return {8,8,28,8,8,9,6};
        case 'u': return {0,17,17,17,17,19,13};
        case 'v': return {0,17,17,17,17,10,4};
        case 'w': return {0,17,17,21,21,21,10};
        case 'x': return {0,17,10,4,10,17,17};
        case 'y': return {0,17,17,17,15,1,14};
        case 'z': return {0,31,2,4,8,16,31};
        default: return {31,17,5,4,4,0,4};
    }
}

// PNG: CRC-32, Adler-32, and zlib stored (uncompressed) DEFLATE blocks.
void u32(std::vector<std::uint8_t>& out, std::uint32_t n) {
    out.push_back(static_cast<std::uint8_t>(n >> 24));
    out.push_back(static_cast<std::uint8_t>(n >> 16));
    out.push_back(static_cast<std::uint8_t>(n >> 8));
    out.push_back(static_cast<std::uint8_t>(n));
}
std::uint32_t crc32(const std::uint8_t* bytes, std::size_t n) {
    std::uint32_t crc = 0xffffffffu;
    for (std::size_t i = 0; i < n; ++i) {
        crc ^= bytes[i];
        for (int j = 0; j != 8; ++j)
            crc = (crc >> 1) ^ ((crc & 1u) ? 0xedb88320u : 0u);
    }
    return ~crc;
}
void chunk(std::vector<std::uint8_t>& png, const char name[5],
           std::vector<std::uint8_t> const& data) {
    u32(png, static_cast<std::uint32_t>(data.size()));
    std::size_t start = png.size();
    for (int i = 0; i != 4; ++i) png.push_back(static_cast<std::uint8_t>(name[i]));
    png.insert(png.end(), data.begin(), data.end());
    u32(png, crc32(png.data()+start, png.size()-start));
}
std::vector<std::uint8_t> zlib_stored(std::vector<std::uint8_t> const& raw) {
    std::vector<std::uint8_t> z;
    z.reserve(raw.size() + 6 + (raw.size()/65535+1)*5);
    z.push_back(0x78); z.push_back(0x01);
    std::uint32_t a = 1, b = 0;
    for (auto c : raw) {
        a = (a+c) % 65521u;
        b = (b+a) % 65521u;
    }
    for (std::size_t pos = 0; pos < raw.size();) {
        std::size_t len = std::min<std::size_t>(65535, raw.size()-pos);
        std::uint16_t n = static_cast<std::uint16_t>(len);
        z.push_back(pos+len == raw.size() ? 1 : 0);
        z.push_back(static_cast<std::uint8_t>(n));
        z.push_back(static_cast<std::uint8_t>(n>>8));
        z.push_back(static_cast<std::uint8_t>(~n));
        z.push_back(static_cast<std::uint8_t>((~n)>>8));
        z.insert(z.end(), raw.begin()+static_cast<std::ptrdiff_t>(pos),
                 raw.begin()+static_cast<std::ptrdiff_t>(pos+len));
        pos += len;
    }
    u32(z, (b<<16)|a);
    return z;
}
} // namespace

// Construct the default graphics canvas through the backend abstraction.
SharedPointer<CanvasInterface> DefaultGraphicsBackend::make_canvas(
        const char*, Nat width, Nat height, Bool) const {
    return std::make_shared<DefaultGraphicsCanvas>(ImageSize2d(width,height));
}

DefaultGraphicsCanvas::DefaultGraphicsCanvas(ImageSize2d const& size)
    : _width(0), _height(0) {
    // Validate the requested plot dimensions *before* adding margins: Nat
    // arithmetic could otherwise wrap for pathological caller input.
    if (size.nx == 0 || size.ny == 0 ||
        size.nx > 16384 || size.ny > 16384)
        throw std::invalid_argument("Invalid or excessively large default graphics canvas");
    _width=size.nx+static_cast<Nat>(LEFT+RIGHT);
    _height=size.ny+static_cast<Nat>(TOP+BOTTOM);
    if (static_cast<std::uint64_t>(_width)*_height > 100000000ULL)
        throw std::invalid_argument("Invalid or excessively large Embedded canvas");
    _pixels.assign(static_cast<std::size_t>(_width)*_height*4, 255);
    _coverage.assign(static_cast<std::size_t>(_width)*_height, 0);
}

ImageSize2d DefaultGraphicsCanvas::size_in_pixels() const {
    return ImageSize2d(_width-LEFT-RIGHT, _height-TOP-BOTTOM);
}

DefaultGraphicsCanvas::Point DefaultGraphicsCanvas::transform(double x, double y) const {
    return {_tx+_sx*x, _ty+_sy*y};
}

// Only initialise() paints the rectangular plot background.  It is always
// opaque and covers the entire (strictly positive) plotting area.  All other
// rectangles are regular paths handled by rasterise_fill().
void DefaultGraphicsCanvas::paint_plot_background() {
    const auto shade=channel(0.97);
    for(int y=TOP; y<static_cast<int>(_height)-BOTTOM; ++y) {
        auto offset=(static_cast<std::size_t>(y)*_width+LEFT)*4;
        for(int x=LEFT; x<static_cast<int>(_width)-RIGHT; ++x,offset+=4) {
            _pixels[offset]=shade;
            _pixels[offset+1]=shade;
            _pixels[offset+2]=shade;
        }
    }
}

// Render the small embedded font as an antialiased stroke skeleton. Connect
// adjacent lit cells in both axes and along isolated diagonal chains. The
// previous implementation drew only horizontal runs; e.g. 'x' and '1' lost
// their vertical/diagonal strokes and degenerated into isolated dots.
void DefaultGraphicsCanvas::text(double x, double y, StringType const& label) {
    constexpr double cell=3.0;
    constexpr double advance=18.0;
    constexpr double width=2.35;
    constexpr double dot_radius=1.20;
    for (char c:label) {
        auto const rows=glyph(c);
        auto has=[&rows](int row,int col) -> bool {
            return row>=0 && row<7 && col>=0 && col<5
                && (rows[static_cast<std::size_t>(row)] & (1u << (4-col)))!=0;
        };
        auto point=[x,y](int row,int col) -> Point {
            return {x+(static_cast<double>(col)+0.5)*cell,
                    y+(static_cast<double>(row)+0.5)*cell};
        };
        for (int row=0;row<7;++row) {
            for (int col=0;col<5;++col) {
                if (!has(row,col)) continue;
                Point p=point(row,col);
                rasterise_disc(p,dot_radius);
                if (has(row,col+1))
                    rasterise_segment(p,point(row,col+1),width);
                if (has(row+1,col))
                    rasterise_segment(p,point(row+1,col),width);
                // Diagonals belong to the skeleton only when neither
                // orthogonal neighbour already connects these cells.
                if (has(row+1,col+1) && !has(row,col+1) && !has(row+1,col))
                    rasterise_segment(p,point(row+1,col+1),width);
                if (has(row+1,col-1) && !has(row,col-1) && !has(row+1,col))
                    rasterise_segment(p,point(row+1,col-1),width);
            }
        }
        composite_coverage({0,0,0,1});
        x+=advance;
    }
}

Void DefaultGraphicsCanvas::initialise(StringType text_x, StringType text_y,
                                double xl, double xu, double yl, double yu) {
    if (!(std::isfinite(xl) && std::isfinite(xu) && std::isfinite(yl) &&
          std::isfinite(yu) && xu>xl && yu>yl))
        throw std::invalid_argument("Default graphics canvas requires finite nonempty axes");
    const int w=static_cast<int>(_width)-LEFT-RIGHT;
    const int h=static_cast<int>(_height)-TOP-BOTTOM;
    const double sx=static_cast<double>(w)/(xu-xl);
    const double sy=-static_cast<double>(h)/(yu-yl);
    const double tx=LEFT-sx*xl;
    const double ty=TOP-sy*yu;
    if (!std::isfinite(sx) || !std::isfinite(sy) ||
        !std::isfinite(tx) || !std::isfinite(ty) ||
        sx<=0.0 || sy>=0.0)
        throw std::invalid_argument("Default graphics canvas requires finite axis mapping");

    std::fill(_pixels.begin(), _pixels.end(), 255);
    _clipped = false;
    _paths.clear();
    const auto xl_text = to_str(xl), xu_text = to_str(xu);
    const auto yl_text = to_str(yl), yu_text = to_str(yu);
    constexpr int step = 18, glyph_h = 21;
    text(LEFT-2, TOP+h+4, xl_text);
    text(LEFT+w-static_cast<int>(xu_text.size())*step-4, TOP+h+4, xu_text);
    text(LEFT+w/2-static_cast<int>(text_x.size())*step/2-3, TOP+h+4, text_x);
    text(LEFT-static_cast<int>(yl_text.size())*step-6, TOP+h+2-glyph_h, yl_text);
    text(LEFT-static_cast<int>(yu_text.size())*step-6, TOP+2, yu_text);
    text(LEFT-static_cast<int>(text_y.size())*step-6, TOP+h/2+2, text_y);
    paint_plot_background();
    _sx=sx;
    _sy=sy;
    _tx=tx;
    _ty=ty;
    _clipped = true;
}

Void DefaultGraphicsCanvas::move_to(double x, double y) {
    Point p = transform(x,y);
    if (std::isfinite(p.x) && std::isfinite(p.y)) _paths.push_back({{p},false});
}

Void DefaultGraphicsCanvas::line_to(double x, double y) {
    Point p = transform(x,y);
    if (!std::isfinite(p.x) || !std::isfinite(p.y)) return;
    if (_paths.empty()) _paths.push_back({{p},false});
    else _paths.back().points.push_back(p);
}

Void DefaultGraphicsCanvas::circle(double x, double y, double radius) {
    Point center = transform(x,y);
    if (!std::isfinite(center.x) || !std::isfinite(center.y) ||
        !std::isfinite(radius) || radius<0) return;
    double rx = std::abs(radius*_sx), ry = std::abs(radius*_sy);
    if (!std::isfinite(rx) || !std::isfinite(ry)) return;
    constexpr double MAX_CHORD_ERROR = 0.10; // sagitta in device pixels
    constexpr double MAX_SEGMENTS = 16384.0;
    double R = std::max(rx,ry);
    double half_angle = R>MAX_CHORD_ERROR ? std::acos(1.0-MAX_CHORD_ERROR/R) : PI;
    int n = static_cast<int>(std::clamp(std::ceil(PI/half_angle),12.0,MAX_SEGMENTS));
    if (_paths.empty()) _paths.emplace_back();
    Path& path = _paths.back();
    path.points.reserve(path.points.size()+static_cast<std::size_t>(n)+1);
    for (int i=0; i<=n; ++i) {
        double theta = 2*PI*i/n;
        path.points.push_back({center.x+rx*std::cos(theta), center.y+ry*std::sin(theta)});
    }
}

Void DefaultGraphicsCanvas::dot(double x, double y) { circle(x,y,_dot_radius/1000.0); }

Void DefaultGraphicsCanvas::set_dot_radius(double radius) { _dot_radius = radius; }

Void DefaultGraphicsCanvas::set_line_width(double width) { _line_width = width; }

Void DefaultGraphicsCanvas::set_line_colour(double r,double g,double b) { _line = {r,g,b,1}; }

Void DefaultGraphicsCanvas::set_fill_colour(double r,double g,double b) {
    _fill.r=r; _fill.g=g; _fill.b=b;
}

Void DefaultGraphicsCanvas::set_fill_opacity(double opacity) { _fill.a=opacity; }

// Coverage is the union of all path fragments, never a sum of their alpha.
void DefaultGraphicsCanvas::mark_sample(int x, int y, int sx, int sy) {
    // Both callers iterate over pixel ranges clamped to the raster bounds.
    // Plot clipping below is still required after initialise().
    if (_clipped && (x<LEFT || x>=static_cast<int>(_width)-RIGHT ||
                     y<TOP || y>=static_cast<int>(_height)-BOTTOM)) return;
    const auto index=static_cast<std::size_t>(y)*_width+static_cast<std::size_t>(x);
    if (_coverage[index]==0) _dirty_pixels.push_back(index);
    _coverage[index] |= static_cast<std::uint16_t>(1u << (sy*4+sx));
}

DefaultGraphicsCanvas::BlendTable const& DefaultGraphicsCanvas::blend_table_for(Paint paint) {
    const Paint normalized{unit(paint.r),unit(paint.g),
                           unit(paint.b),unit(paint.a)};
    BlendCacheEntry* entry=nullptr;
    for(auto& candidate:_blend_cache) {
        if(candidate.valid &&
           candidate.paint.r==normalized.r && candidate.paint.g==normalized.g &&
           candidate.paint.b==normalized.b && candidate.paint.a==normalized.a) {
            entry=&candidate;
            break;
        }
    }
    if(entry==nullptr) {
        // Prefer an unused slot, then replace the least recently used paint.
        entry=&_blend_cache.front();
        for(auto& candidate:_blend_cache) {
            if(!candidate.valid) {
                entry=&candidate;
                break;
            }
            if(candidate.last_used<entry->last_used) entry=&candidate;
        }
        const double source_byte[3]={
            255.0*normalized.r,255.0*normalized.g,255.0*normalized.b
        };
        for(std::size_t n=0;n<=16;++n) {
            const double alpha=normalized.a*(static_cast<double>(n)/16.0);
            for(std::size_t color=0;color<3;++color)
                for(std::size_t old=0;old<256;++old)
                    entry->table[n][color][old]=composite_channel(
                        static_cast<std::uint8_t>(old),source_byte[color],alpha);
        }
        entry->paint=normalized;
        entry->valid=true;
    }
    entry->last_used=++_blend_cache_clock;
    return entry->table;
}

void DefaultGraphicsCanvas::composite_coverage(Paint paint) {
    if(_dirty_pixels.empty()) return;
    auto const& lookup=blend_table_for(paint);
    for(auto index:_dirty_pixels) {
        // Every dirty pixel has a nonzero accumulated subpixel mask.
        const auto bits=_coverage[index];
        const auto& table=lookup[static_cast<std::size_t>(std::popcount(bits))];
        const std::size_t pos=index*4;
        _pixels[pos]=table[0][_pixels[pos]];
        _pixels[pos+1]=table[1][_pixels[pos+1]];
        _pixels[pos+2]=table[2][_pixels[pos+2]];
        _coverage[index]=0;
    }
    _dirty_pixels.clear();
}
// Mark a half-open range of 4x4 subpixel samples in a scanline. Stroke
// boundaries are closed (inclusive_end); fills use the original half-open
// winding intervals, preserving the previous 4x4 antialiasing convention.
void DefaultGraphicsCanvas::rasterise_span(double from,double to,int y,int sy,
                                    bool inclusive_end) {
    const double left=_clipped ? LEFT : 0;
    const double right=_clipped ? static_cast<double>(_width)-RIGHT : _width;
    if (!std::isfinite(from) || !std::isfinite(to) ||
        to<left || from>=right || y<0 || y>=static_cast<int>(_height))
        return;
    from=std::clamp(from,left,right);
    to=std::clamp(to,left,right);
    const int first=std::max(0,static_cast<int>(std::ceil(from*4.0-0.5)));
    const int last=std::min(static_cast<int>(_width)*4,
        static_cast<int>(inclusive_end ? std::floor(to*4.0-0.5)+1.0 :
                                         std::ceil(to*4.0-0.5)));
    for(int sample=first;sample<last;) {
        const int x=sample/4;
        const int stop=std::min(last,(x+1)*4);
        const unsigned n=static_cast<unsigned>(stop-sample);
        const unsigned shift=static_cast<unsigned>(sy*4+sample%4);
        const auto bits=static_cast<std::uint16_t>(((1u<<n)-1u)<<shift);
        const auto index=static_cast<std::size_t>(y)*_width+
                         static_cast<std::size_t>(x);
        if (_coverage[index]==0) _dirty_pixels.push_back(index);
        _coverage[index]|=bits;
        sample=stop;
    }
}

// A butt-capped stroke is the intersection of nested_pair pairs of half-planes:
//  0 <= dot(pixel-a,b-a) <= length_squared
//  |cross(pixel-a,b-a)| <= (width/2)*length.
// Each condition is linear in x for a fixed subpixel scanline y. Intersect
// the intervals once per scanline instead of evaluating 16 predicates for
// every candidate pixel. The subpixel mask remains bit-exact with the former
// sampling rasterizer (apart from unavoidable floating-point boundary ties).
void DefaultGraphicsCanvas::rasterise_segment(Point a,Point b,double width) {
    if (!(width>0) || !std::isfinite(width)) return;
    const double dx=b.x-a.x,dy=b.y-a.y;
    const double length2=dx*dx+dy*dy;
    if (!std::isfinite(length2) || length2<=0) return;
    const double radius=width/2.0;
    const double cross_limit=radius*std::sqrt(length2);
    if (!std::isfinite(cross_limit)) return;

    if(dx==0.0 || dy==0.0) {
        // A vertical or horizontal butt-capped segment is exactly a
        // rectangle in device space. Build a single 4x4 pixel mask per
        // scanline instead of intersecting half-planes four times per row.
        // Boundaries are inclusive, just like the general segment rasterizer.
        const double xmin=dy==0.0 ? std::min(a.x,b.x) : a.x-radius;
        const double xmax=dy==0.0 ? std::max(a.x,b.x) : a.x+radius;
        const double ymin=dx==0.0 ? std::min(a.y,b.y) : a.y-radius;
        const double ymax=dx==0.0 ? std::max(a.y,b.y) : a.y+radius;
        if(std::isfinite(xmin) && std::isfinite(xmax) &&
           std::isfinite(ymin) && std::isfinite(ymax)) {
            const double left=_clipped ? LEFT : 0;
            const double right=_clipped ? static_cast<double>(_width)-RIGHT
                                        : static_cast<double>(_width);
            const double top=_clipped ? TOP : 0;
            const double bottom=_clipped ? static_cast<double>(_height)-BOTTOM
                                         : static_cast<double>(_height);
            if(xmax<=left || xmin>=right || ymax<=top || ymin>=bottom)
                return;
            const double x_from=std::clamp(xmin,left,right);
            const double x_to=std::clamp(xmax,left,right);
            const int first=static_cast<int>(std::ceil(x_from*4.0-0.5));
            const int last=static_cast<int>(std::floor(x_to*4.0-0.5)+1.0);
            if(first>=last) return;
            const int y0=static_cast<int>(std::floor(
                std::clamp(ymin,top,bottom)));
            const int y1=static_cast<int>(std::ceil(
                std::clamp(ymax,top,bottom)));
            const int x0=first/4;
            const int x1=(last+3)/4;
            for(int y=y0;y<y1;++y) {
                unsigned vertical_mask=0;
                for(int sy=0;sy<4;++sy) {
                    const double sample_y=y+(sy+0.5)/4.0;
                    if(sample_y>=ymin && sample_y<=ymax)
                        vertical_mask|=0xfu<<(sy*4);
                }
                if(!vertical_mask) continue;
                const std::size_t row=static_cast<std::size_t>(y)*_width;
                for(int x=x0;x<x1;++x) {
                    const int begin=std::clamp(first-4*x,0,4);
                    const int end=std::clamp(last-4*x,0,4);
                    const unsigned horizontal_mask=
                        ((1u<<end)-1u)&~((1u<<begin)-1u);
                    // This pixel has a horizontal sample, and the vertical
                    // sample mask is already known to be nonzero.
                    const auto mask=static_cast<std::uint16_t>(
                        (horizontal_mask*0x1111u)&vertical_mask);
                    const auto index=row+static_cast<std::size_t>(x);
                    if(_coverage[index]==0) _dirty_pixels.push_back(index);
                    _coverage[index]|=mask;
                }
            }
            return;
        }
    }

    const double bottom=std::min(a.y,b.y)-radius-1;
    const double top=std::max(a.y,b.y)+radius+1;
    const double clip_bottom=_clipped ? TOP : 0;
    const double clip_top=_clipped ? static_cast<double>(_height)-BOTTOM : _height;
    if (top<=clip_bottom || bottom>=clip_top) return;
    const int y0=static_cast<int>(std::floor(std::clamp(bottom,clip_bottom,clip_top)));
    const int y1=static_cast<int>(std::ceil(std::clamp(top,clip_bottom,clip_top)));
    for (int y=y0;y<y1;++y) for (int sy=0;sy<4;++sy) {
        const double py=y+(sy+0.5)/4.0;
        double from=_clipped ? LEFT : 0;
        double to=_clipped ? static_cast<double>(_width)-RIGHT : _width;
        // Restrict x for coefficient*x+offset in [lower,upper].
        auto restrict_x=[&](double coefficient,double offset,
                            double lower,double upper) {
            if(coefficient==0.0)
                return offset>=lower && offset<=upper;
            double lo=(lower-offset)/coefficient;
            double hi=(upper-offset)/coefficient;
            if(lo>hi) std::swap(lo,hi);
            from=std::max(from,lo);
            to=std::min(to,hi);
            return from<=to;
        };
        if (!restrict_x(dx,dy*(py-a.y)-dx*a.x,0,length2) ||
            !restrict_x(dy,-dy*a.x-dx*(py-a.y),-cross_limit,cross_limit))
            continue;
        rasterise_span(from,to,y,sy,true);
    }
}
// Round endpoints restore readable glyph strokes while preserving the
// butt caps of ordinary drawing paths.
void DefaultGraphicsCanvas::rasterise_disc(Point center,double radius) {
    if (!(radius>0) || !std::isfinite(center.x) || !std::isfinite(center.y)) return;
    const double xmin=center.x-radius, xmax=center.x+radius;
    const double ymin=center.y-radius, ymax=center.y+radius;
    if (xmax<0 || ymax<0 || xmin>=_width || ymin>=_height) return;
    const int x0=static_cast<int>(std::max(0.0,std::floor(xmin)));
    const int x1=static_cast<int>(std::min<double>(_width,std::ceil(xmax)));
    const int y0=static_cast<int>(std::max(0.0,std::floor(ymin)));
    const int y1=static_cast<int>(std::min<double>(_height,std::ceil(ymax)));
    for (int yy=y0;yy<y1;++yy) for (int xx=x0;xx<x1;++xx)
        for (int sy=0;sy<4;++sy) for (int sx=0;sx<4;++sx) {
            const double dx=xx+(sx+0.5)/4.0-center.x;
            const double dy=yy+(sy+0.5)/4.0-center.y;
            if (dx*dx+dy*dy<=radius*radius)
                mark_sample(xx,yy,sx,sy);
        }
}

void DefaultGraphicsCanvas::rasterise_triangle(Point a, Point b, Point c) {
    const double xmin=std::min({a.x,b.x,c.x}),xmax=std::max({a.x,b.x,c.x});
    const double ymin=std::min({a.y,b.y,c.y}),ymax=std::max({a.y,b.y,c.y});
    if(!std::isfinite(xmin)||!std::isfinite(xmax)||!std::isfinite(ymin)||!std::isfinite(ymax) ||
       xmax<0||ymax<0||xmin>=_width||ymin>=_height) return;
    const double area=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    if(std::abs(area)<1e-15) return;
    const int x0=static_cast<int>(std::max(0.0,std::floor(xmin)));
    const int x1=static_cast<int>(std::min<double>(_width,std::ceil(xmax)));
    const int y0=static_cast<int>(std::max(0.0,std::floor(ymin)));
    const int y1=static_cast<int>(std::min<double>(_height,std::ceil(ymax)));
    for(int y=y0;y<y1;++y) for(int x=x0;x<x1;++x)
        for(int sy=0;sy<4;++sy) for(int sx=0;sx<4;++sx) {
            Point p{x+(sx+0.5)/4,y+(sy+0.5)/4};
            double s1=(b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x);
            double s2=(c.x-b.x)*(p.y-b.y)-(c.y-b.y)*(p.x-b.x);
            double s3=(a.x-c.x)*(p.y-c.y)-(a.y-c.y)*(p.x-c.x);
            if((s1>=0&&s2>=0&&s3>=0)||(s1<=0&&s2<=0&&s3<=0))
                mark_sample(x,y,sx,sy);
        }
}

void DefaultGraphicsCanvas::rasterise_join(Point a,Point b,Point c,double half_width) {
    double ux=b.x-a.x,uy=b.y-a.y,vx=c.x-b.x,vy=c.y-b.y;
    double lu=std::hypot(ux,uy),lv=std::hypot(vx,vy);
    if(lu==0||lv==0) return;
    ux/=lu;uy/=lu;vx/=lv;vy/=lv;
    double cross=ux*vy-uy*vx;
    if(std::abs(cross)<1e-12) return;
    double side=(cross>0 ? -1.0:1.0);
    Point p{b.x+side*(-uy)*half_width,b.y+side*ux*half_width};
    Point q{b.x+side*(-vy)*half_width,b.y+side*vx*half_width};
    // The determinant equals -cross: the earlier parallel-segment check
    // already rules out a vanishing denominator. Limit miter joins to 10.
    double determinant=ux*(-vy)-(-vx)*uy;
    double t=((q.x-p.x)*(-vy)-(-vx)*(q.y-p.y))/determinant;
    Point m{p.x+t*ux,p.y+t*uy};
    // Bridge the ends of the adjoining butt-capped segments first.
    // Without this triangle the outer miter tip is visually detached.
    rasterise_triangle(p,b,q);
    if(std::hypot(m.x-b.x,m.y-b.y) <= 10*half_width)
        rasterise_triangle(p,m,q);
}
// Rasterise horizontal intervals selected by a nonzero-winding scanline.
//
// Maintain only edges intersecting the current row. The previous implementation
// tested every edge for every subpixel, which was prohibitively expensive for
// paving plots containing many polygons.
void DefaultGraphicsCanvas::rasterise_fill() {
    struct Rectangle {
        double xmin, ymin, xmax, ymax;
        bool positive_winding;
    };
    // Recognise only complete, strictly axis-aligned rectangles. A path may
    // repeat its starting vertex, but must contain all four distinct corners
    // exactly once. Concave and degenerate contours use the general scanline.
    auto identify_rectangle=[](Path const& path,Rectangle& result) {
        auto const& points=path.points;
        if(points.size()!=4 && points.size()!=5) return false;
        if(points.size()==5 &&
           (points[0].x!=points[4].x || points[0].y!=points[4].y))
            return false;
        double xmin=std::numeric_limits<double>::infinity();
        double ymin=xmin,xmax=-xmin,ymax=-xmin;
        for(std::size_t i=0;i<4;++i) {
            auto const p=points[i];
            if(!std::isfinite(p.x) || !std::isfinite(p.y)) return false;
            xmin=std::min(xmin,p.x);
            xmax=std::max(xmax,p.x);
            ymin=std::min(ymin,p.y);
            ymax=std::max(ymax,p.y);
        }
        if(xmax<=xmin || ymax<=ymin) return false;
        std::array<bool,4> corners{};
        for(std::size_t i=0;i<4;++i) {
            const auto a=points[i],b=points[(i+1)%4];
            if(!((a.x==b.x && a.y!=b.y) ||
                 (a.y==b.y && a.x!=b.x))) return false;
            const auto index=static_cast<std::size_t>(
                (a.x==xmax ? 1 : 0)+(a.y==ymax ? 2 : 0));
            if(corners[index]) return false;
            corners[index]=true;
        }
        const auto a=points[0],b=points[1],c=points[2];
        const double cross=(b.x-a.x)*(c.y-b.y)-
                           (b.y-a.y)*(c.x-b.x);
        if(cross==0 || std::isnan(cross)) return false;
        result={xmin,ymin,xmax,ymax,cross>0};
        return true;
    };

    // A single rectangle and a strictly nested pair (filled frame) need no
    // edge sorting. The second contour either creates a hole if its winding
    // is opposite, or leaves the entire outer rectangle filled.
    Rectangle outer{},inner{};
    const bool one=_paths.size()==1 &&
        identify_rectangle(_paths.front(),outer);
    const bool nested_pair=_paths.size()==2 &&
        identify_rectangle(_paths.front(),outer) &&
        identify_rectangle(_paths.back(),inner);
    if(one || nested_pair) {
        if(nested_pair) {
            // The same fill works regardless of subpath order.
            const double outer_area=(outer.xmax-outer.xmin)*
                                    (outer.ymax-outer.ymin);
            const double inner_area=(inner.xmax-inner.xmin)*
                                    (inner.ymax-inner.ymin);
            if(inner_area>outer_area) std::swap(outer,inner);
        }
        const bool contained=!nested_pair ||
            (inner.xmin>outer.xmin && inner.xmax<outer.xmax &&
             inner.ymin>outer.ymin && inner.ymax<outer.ymax);
        if(contained) {
            const bool hole=nested_pair &&
                outer.positive_winding!=inner.positive_winding;
            const double clip_left=_clipped ? LEFT : 0;
            const double clip_right=_clipped ?
                static_cast<double>(_width)-RIGHT : _width;
            const double clip_top=_clipped ? TOP : 0;
            const double clip_bottom=_clipped ?
                static_cast<double>(_height)-BOTTOM : _height;
            if(outer.xmax<=clip_left || outer.xmin>=clip_right ||
               outer.ymax<=clip_top || outer.ymin>=clip_bottom)
                return;

            // x sample positions are x+(sx+0.5)/4. Match the half-open
            // ranges used by the general winding scanline.
            auto sample_bounds=[clip_left,clip_right](Rectangle const& rect) {
                const double left=std::clamp(rect.xmin,clip_left,clip_right);
                const double right=std::clamp(rect.xmax,clip_left,clip_right);
                return std::pair<int,int>{
                    static_cast<int>(std::ceil(left*4.0-0.5)),
                    static_cast<int>(std::ceil(right*4.0-0.5))
                };
            };
            const auto [first,last]=sample_bounds(outer);
            if(first>=last) return;
            const auto [inner_first,inner_last]=sample_bounds(inner);
            const int y0=static_cast<int>(std::floor(
                std::clamp(outer.ymin,clip_top,clip_bottom)));
            const int y1=static_cast<int>(std::ceil(
                std::clamp(outer.ymax,clip_top,clip_bottom)));
            auto horizontal_bits=[](int first_sample,int last_sample,int x) {
                const int start=std::clamp(first_sample-4*x,0,4);
                const int stop=std::clamp(last_sample-4*x,0,4);
                return ((1u<<stop)-1u) & ~((1u<<start)-1u);
            };

            // Compose directly: a rectangle's final subpixel mask is known
            // in one pass, so it need not be stored in the coverage/dirty
            // buffers for a second traversal.
            auto const& lookup=blend_table_for(_fill);
            for(int y=y0;y<y1;++y) {
                unsigned outer_vertical=0,inner_vertical=0;
                for(int sy=0;sy<4;++sy) {
                    const double sample_y=y+(sy+0.5)/4.0;
                    if(sample_y>=outer.ymin && sample_y<outer.ymax)
                        outer_vertical|=0xfu<<(sy*4);
                    if(hole && sample_y>=inner.ymin && sample_y<inner.ymax)
                        inner_vertical|=0xfu<<(sy*4);
                }
                if(!outer_vertical) continue;
                const int x_begin=first/4;
                const int x_end=(last+3)/4;
                const std::size_t row_offset=
                    static_cast<std::size_t>(y)*_width*4;

                // Every sample of a strictly interior pixel in a full-height
                // row is covered. Only the first and last pixels need a
                // horizontal edge mask; avoid masking and popcount elsewhere.
                if(!hole && outer_vertical==0xffffu && x_end-x_begin>=3) {
                    auto blend_edge=[&](int x) {
                        const unsigned horizontal=
                            horizontal_bits(first,last,x)*0x1111u;
                        if(!horizontal) return;
                        const auto& table=lookup[
                            static_cast<std::size_t>(std::popcount(horizontal))];
                        const std::size_t offset=
                            row_offset+static_cast<std::size_t>(x)*4;
                        _pixels[offset]=table[0][_pixels[offset]];
                        _pixels[offset+1]=table[1][_pixels[offset+1]];
                        _pixels[offset+2]=table[2][_pixels[offset+2]];
                    };
                    blend_edge(x_begin);
                    const auto& solid=lookup[16];
                    std::size_t offset=
                        row_offset+static_cast<std::size_t>(x_begin+1)*4;
                    for(int x=x_begin+1;x<x_end-1;++x,offset+=4) {
                        _pixels[offset]=solid[0][_pixels[offset]];
                        _pixels[offset+1]=solid[1][_pixels[offset+1]];
                        _pixels[offset+2]=solid[2][_pixels[offset+2]];
                    }
                    blend_edge(x_end-1);
                    continue;
                }

                for(int x=x_begin;x<x_end;++x) {
                    const unsigned outer_horizontal=
                        horizontal_bits(first,last,x);
                    const unsigned outside=outer_horizontal*0x1111u;
                    unsigned mask=outside&outer_vertical;
                    if(hole && mask && inner_vertical) {
                        const unsigned inner_horizontal=
                            horizontal_bits(inner_first,inner_last,x);
                        mask&=~((inner_horizontal*0x1111u)&inner_vertical);
                    }
                    if(!mask) continue;
                    const auto& table=lookup[
                        static_cast<std::size_t>(std::popcount(mask))];
                    const std::size_t offset=
                        row_offset+static_cast<std::size_t>(x)*4;
                    _pixels[offset]=table[0][_pixels[offset]];
                    _pixels[offset+1]=table[1][_pixels[offset+1]];
                    _pixels[offset+2]=table[2][_pixels[offset+2]];
                }
            }
            return;
        }
    }
    struct Edge {
        double ymin, ymax, x_at_ymin, dx_dy;
        int delta;
    };
    std::vector<Edge> edges;
    double ymin=std::numeric_limits<double>::infinity();
    double ymax=-std::numeric_limits<double>::infinity();
    for (auto const& path : _paths) {
        if (path.points.size()<3) continue;
        for (std::size_t i=0; i<path.points.size(); ++i) {
            Point p=path.points[i];
            Point q=path.points[(i+1)%path.points.size()];
            if (p.y == q.y || !std::isfinite(p.x) || !std::isfinite(q.x) ||
                !std::isfinite(p.y) || !std::isfinite(q.y)) continue;
            int delta=1;
            if (p.y>q.y) { std::swap(p,q); delta=-1; }
            double slope=(q.x-p.x)/(q.y-p.y);
            if (!std::isfinite(slope)) continue;
            ymin=std::min(ymin,p.y);
            ymax=std::max(ymax,q.y);
            edges.push_back({p.y,q.y,p.x,slope,delta});
        }
    }
    if (edges.empty()) return;
    double clip_top=_clipped ? TOP : 0;
    double clip_bottom=_clipped ? static_cast<double>(_height)-BOTTOM : _height;
    if (ymax<=clip_top || ymin>=clip_bottom) return;
    const int y0=static_cast<int>(std::floor(std::clamp(ymin,clip_top,clip_bottom)));
    const int y1=static_cast<int>(std::ceil(std::clamp(ymax,clip_top,clip_bottom)));
    std::sort(edges.begin(),edges.end(),
              [](Edge const& a,Edge const& b) { return a.ymin<b.ymin; });
    std::vector<Edge const*> active;
    std::vector<std::pair<double,int>> crossings;
    active.reserve(edges.size());
    crossings.reserve(edges.size());
    std::size_t next=0;
    // A single sample row contains 4 horizontally offset samples per pixel.
    // A contiguous interval maps to an OR of bitmasks, not 4 tests per pixel.
    for (int y=y0; y<y1; ++y) for (int sy=0;sy<4;++sy) {
        const double scan_y=y+(sy+0.5)/4.0;
        // Edges are active on the half-open vertical interval [ymin,ymax).
        while(next<edges.size() && edges[next].ymin<=scan_y)
            active.push_back(&edges[next++]);
        active.erase(std::remove_if(active.begin(),active.end(),
                       [scan_y](Edge const* edge){return scan_y>=edge->ymax;}),
                     active.end());
        crossings.clear();
        for (Edge const* edge:active) {
            double x=edge->x_at_ymin+(scan_y-edge->ymin)*edge->dx_dy;
            if (std::isfinite(x)) crossings.emplace_back(x,edge->delta);
        }
        std::sort(crossings.begin(),crossings.end(),
                  [](auto const& a,auto const& b){return a.first<b.first;});
        int winding=0;
        for (std::size_t i=0;i+1<crossings.size();++i) {
            winding+=crossings[i].second;
            if (winding!=0 && crossings[i].first<crossings[i+1].first)
                rasterise_span(crossings[i].first,crossings[i+1].first,y,sy,false);
        }
    }
}

Void DefaultGraphicsCanvas::stroke() {
    if(_line_width>0 && std::isfinite(_line_width)) {
        for(auto const& path:_paths) {
            auto const& pts=path.points;
            if(pts.size()<2) continue;
            for(std::size_t i=1;i<pts.size();++i)
                rasterise_segment(pts[i-1],pts[i],_line_width);
            if(path.closed) rasterise_segment(pts.back(),pts.front(),_line_width);
            for(std::size_t i=1;i+1<pts.size();++i)
                rasterise_join(pts[i-1],pts[i],pts[i+1],_line_width/2);
            if(path.closed) {
                rasterise_join(pts[pts.size()-2],pts.back(),pts.front(),_line_width/2);
                rasterise_join(pts.back(),pts.front(),pts[1],_line_width/2);
            }
        }
        composite_coverage(_line);
    }
    _paths.clear();
}

Void DefaultGraphicsCanvas::fill() {
    if(unit(_fill.a)>0) {
        rasterise_fill();
        composite_coverage(_fill);
    }
    // Preserve the path through filling, then stroke and consume it.
    stroke();
}

Void DefaultGraphicsCanvas::finalise() {
    bool old_clip=_clipped;
    _clipped=false;
    _paths.clear();
    const double x0=LEFT, x1=_width-RIGHT, y0=TOP, y1=_height-BOTTOM;
    _paths.push_back({{{x0,y1},{x1,y1},{x1,y0},{x0,y0}},true});
    auto old_width=_line_width; auto old_line=_line;
    _line_width=2; _line={0,0,0,1};
    stroke();
    _line_width=old_width; _line=old_line; _clipped=old_clip;
}

Void DefaultGraphicsCanvas::write(const char* name) const {
    if (!name) throw std::invalid_argument("Null PNG filename");
    std::string filename(name);
    if (filename.find('.')==std::string::npos) filename += ".png";
    std::vector<std::uint8_t> raw;
    raw.reserve(static_cast<std::size_t>(_height)*(1+static_cast<std::size_t>(_width)*4));
    for (Nat y=0; y<_height; ++y) {
        raw.push_back(0); // PNG filter: None
        auto pos=static_cast<std::size_t>(y)*_width*4;
        raw.insert(raw.end(),_pixels.begin()+static_cast<std::ptrdiff_t>(pos),
                    _pixels.begin()+static_cast<std::ptrdiff_t>(pos+_width*4));
    }
    std::vector<std::uint8_t> png{137,80,78,71,13,10,26,10};
    std::vector<std::uint8_t> ihdr;
    u32(ihdr,_width);u32(ihdr,_height);
    ihdr.insert(ihdr.end(),{8,6,0,0,0}); // RGBA8
    chunk(png,"IHDR",ihdr);
    chunk(png,"IDAT",zlib_stored(raw));
    chunk(png,"IEND",{});
    std::ofstream stream(filename,std::ios::binary);
    if (!stream || !stream.write(reinterpret_cast<const char*>(png.data()),
                                  static_cast<std::streamsize>(png.size())))
        throw std::runtime_error("Unable to write default graphics PNG: "+filename);
}
} // namespace Ariadne
