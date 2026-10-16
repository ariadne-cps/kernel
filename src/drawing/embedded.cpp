/*
 * Ariadne - free software under the GNU General Public License version 3
 * or (at your option) any later version. No warranty.
 */
#include "drawing/embedded.hpp"
#include "utility/stlio.hpp"
#include "utility/standard.hpp"

#include <algorithm>
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

SharedPointer<CanvasInterface> EmbeddedGraphicsBackend::make_canvas(
        const char*, Nat width, Nat height, Bool) const {
    return std::make_shared<EmbeddedCanvas>(ImageSize2d(width,height));
}

EmbeddedCanvas::EmbeddedCanvas(ImageSize2d const& size)
    : _width(size.nx + LEFT + RIGHT), _height(size.ny + TOP + BOTTOM) {
    if (size.nx == 0 || size.ny == 0 ||
        size.nx > 16384 || size.ny > 16384 ||
        static_cast<std::uint64_t>(_width)*_height > 100000000ULL)
        throw std::invalid_argument("Invalid or excessively large Embedded canvas");
    _pixels.assign(static_cast<std::size_t>(_width)*_height*4, 255);
}

ImageSize2d EmbeddedCanvas::size_in_pixels() const {
    return ImageSize2d(_width-LEFT-RIGHT, _height-TOP-BOTTOM);
}

EmbeddedCanvas::Point EmbeddedCanvas::transform(double x, double y) const {
    return {_tx+_sx*x, _ty+_sy*y};
}

void EmbeddedCanvas::blend(int x, int y, Paint paint, double coverage) {
    if (x < 0 || y < 0 || x >= static_cast<int>(_width) || y >= static_cast<int>(_height)) return;
    if (_clipped && (x < LEFT || x >= static_cast<int>(_width)-RIGHT ||
                     y < TOP || y >= static_cast<int>(_height)-BOTTOM)) return;
    double alpha = unit(paint.a) * unit(coverage);
    if (alpha <= 0) return;
    std::size_t index = (static_cast<std::size_t>(y)*_width+x)*4;
    _pixels[index] = channel((1-alpha)*(_pixels[index]/255.0)+alpha*unit(paint.r));
    _pixels[index+1] = channel((1-alpha)*(_pixels[index+1]/255.0)+alpha*unit(paint.g));
    _pixels[index+2] = channel((1-alpha)*(_pixels[index+2]/255.0)+alpha*unit(paint.b));
}

void EmbeddedCanvas::rectangle(int x0, int y0, int x1, int y1, Paint paint) {
    for (int y = std::max(y0,0); y < std::min(y1,static_cast<int>(_height)); ++y)
        for (int x = std::max(x0,0); x < std::min(x1,static_cast<int>(_width)); ++x)
            blend(x,y,paint,1);
}

void EmbeddedCanvas::text(double x, double y, StringType const& label) {
    constexpr int scale = 3;
    double advance = 6*scale;
    for (char c : label) {
        auto rows = glyph(c);
        for (int row = 0; row != 7; ++row)
            for (int col = 0; col != 5; ++col)
                if (rows[row] & (1u << (4-col)))
                    rectangle(static_cast<int>(x+col*scale), static_cast<int>(y+row*scale),
                              static_cast<int>(x+(col+1)*scale), static_cast<int>(y+(row+1)*scale),
                              {0,0,0,1});
        x += advance;
    }
}

Void EmbeddedCanvas::initialise(StringType text_x, StringType text_y,
                                double xl, double xu, double yl, double yu) {
    if (!(std::isfinite(xl) && std::isfinite(xu) && std::isfinite(yl) &&
          std::isfinite(yu) && xu>xl && yu>yl))
        throw std::invalid_argument("Embedded canvas requires finite nonempty axes");
    std::fill(_pixels.begin(), _pixels.end(), 255);
    _clipped = false;
    _paths.clear();
    int w = static_cast<int>(_width)-LEFT-RIGHT;
    int h = static_cast<int>(_height)-TOP-BOTTOM;
    const auto xl_text = to_str(xl), xu_text = to_str(xu);
    const auto yl_text = to_str(yl), yu_text = to_str(yu);
    constexpr int step = 18, glyph_h = 21;
    text(LEFT-2, TOP+h+4, xl_text);
    text(LEFT+w-static_cast<int>(xu_text.size())*step-4, TOP+h+4, xu_text);
    text(LEFT+w/2-static_cast<int>(text_x.size())*step/2-3, TOP+h+4, text_x);
    text(LEFT-static_cast<int>(yl_text.size())*step-6, TOP+h+2-glyph_h, yl_text);
    text(LEFT-static_cast<int>(yu_text.size())*step-6, TOP+2, yu_text);
    text(LEFT-static_cast<int>(text_y.size())*step-6, TOP+h/2+2, text_y);
    rectangle(LEFT, TOP, LEFT+w, TOP+h, {0.97,0.97,0.97,1.0});
    _sx = static_cast<double>(w)/(xu-xl);
    _sy = -static_cast<double>(h)/(yu-yl);
    _tx = LEFT - _sx*xl;
    _ty = TOP - _sy*yu;
    _clipped = true;
}

Void EmbeddedCanvas::move_to(double x, double y) {
    Point p = transform(x,y);
    if (std::isfinite(p.x) && std::isfinite(p.y)) _paths.push_back({{p},false});
}
Void EmbeddedCanvas::line_to(double x, double y) {
    Point p = transform(x,y);
    if (!std::isfinite(p.x) || !std::isfinite(p.y)) return;
    if (_paths.empty()) _paths.push_back({{p},false});
    else _paths.back().points.push_back(p);
}
Void EmbeddedCanvas::circle(double x, double y, double radius) {
    Point center = transform(x,y);
    if (!std::isfinite(center.x) || !std::isfinite(center.y) ||
        !std::isfinite(radius) || radius<0) return;
    double rx = std::abs(radius*_sx), ry = std::abs(radius*_sy);
    if (!std::isfinite(rx) || !std::isfinite(ry)) return;
    int n = std::clamp(static_cast<int>(std::ceil(PI*std::sqrt(std::max(rx,ry)+1.0))),24,1024);
    Path path;
    path.closed = true;
    path.points.reserve(n);
    for (int i=0; i<n; ++i) {
        double theta = 2*PI*i/n;
        path.points.push_back({center.x+rx*std::cos(theta), center.y+ry*std::sin(theta)});
    }
    _paths.push_back(std::move(path));
}
Void EmbeddedCanvas::dot(double x, double y) { circle(x,y,_dot_radius/1000.0); }
Void EmbeddedCanvas::set_dot_radius(double radius) { _dot_radius = radius; }
Void EmbeddedCanvas::set_line_width(double width) { _line_width = width; }
Void EmbeddedCanvas::set_line_colour(double r,double g,double b) { _line = {r,g,b,1}; }
Void EmbeddedCanvas::set_fill_colour(double r,double g,double b) {
    _fill.r=r; _fill.g=g; _fill.b=b;
}
Void EmbeddedCanvas::set_fill_opacity(double opacity) { _fill.a=opacity; }

// 4x4 subpixel coverage for device-space strokes, width measured in pixels.
void EmbeddedCanvas::segment(Point a, Point b, double width, Paint paint) {
    if (!(width>0) || !std::isfinite(width)) return;
    double radius = width/2;
    double dx=b.x-a.x, dy=b.y-a.y, length2=dx*dx+dy*dy;
    if (!std::isfinite(length2)) return;
    double xmin = std::min(a.x,b.x)-radius-1, xmax=std::max(a.x,b.x)+radius+1;
    double ymin = std::min(a.y,b.y)-radius-1, ymax=std::max(a.y,b.y)+radius+1;
    if (xmax<0 || ymax<0 || xmin>=_width || ymin>=_height) return;
    int x0=std::max(0,static_cast<int>(std::floor(std::max(0.0,xmin))));
    int x1=std::min(static_cast<int>(_width),static_cast<int>(std::ceil(std::min<double>(_width,xmax))));
    int y0=std::max(0,static_cast<int>(std::floor(std::max(0.0,ymin))));
    int y1=std::min(static_cast<int>(_height),static_cast<int>(std::ceil(std::min<double>(_height,ymax))));
    for (int y=y0; y<y1; ++y)
        for (int x=x0; x<x1; ++x) {
            int covered=0;
            for (int sy=0;sy<4;++sy) for (int sx=0;sx<4;++sx) {
                double px=x+(sx+0.5)/4, py=y+(sy+0.5)/4;
                double t=length2>0 ? std::clamp(((px-a.x)*dx+(py-a.y)*dy)/length2,0.0,1.0) : 0.0;
                double ex=px-(a.x+t*dx), ey=py-(a.y+t*dy);
                if (ex*ex+ey*ey <= radius*radius) ++covered;
            }
            if (covered) blend(x,y,paint,covered/16.0);
        }
}

// Nonzero winding fill, using 4x4 supersampling and implicit closing edges.
void EmbeddedCanvas::polygon(Path const& path, Paint paint) {
    if (path.points.size()<3 || unit(paint.a)==0) return;
    double minx=std::numeric_limits<double>::infinity(), maxx=-minx;
    double miny=minx, maxy=-minx;
    for (Point p: path.points) {
        minx=std::min(minx,p.x); maxx=std::max(maxx,p.x);
        miny=std::min(miny,p.y); maxy=std::max(maxy,p.y);
    }
    if (!std::isfinite(minx) || !std::isfinite(maxx) || !std::isfinite(miny) || !std::isfinite(maxy) ||
        maxx<0 || maxy<0 || minx>=_width || miny>=_height) return;
    int x0=static_cast<int>(std::max(0.0,std::floor(minx)));
    int x1=static_cast<int>(std::min<double>(_width,std::ceil(maxx)));
    int y0=static_cast<int>(std::max(0.0,std::floor(miny)));
    int y1=static_cast<int>(std::min<double>(_height,std::ceil(maxy)));
    for (int y=y0; y<y1; ++y)
        for (int x=x0; x<x1; ++x) {
            int covered=0;
            for (int sy=0;sy<4;++sy) for (int sx=0;sx<4;++sx) {
                double px=x+(sx+0.5)/4, py=y+(sy+0.5)/4;
                int winding=0;
                for (std::size_t i=0;i<path.points.size();++i) {
                    Point p=path.points[i], q=path.points[(i+1)%path.points.size()];
                    double cross=(q.x-p.x)*(py-p.y)-(q.y-p.y)*(px-p.x);
                    if (p.y<=py && q.y>py && cross>0) ++winding;
                    if (p.y>py && q.y<=py && cross<0) --winding;
                }
                if (winding) ++covered;
            }
            if (covered) blend(x,y,paint,covered/16.0);
        }
}

Void EmbeddedCanvas::stroke() {
    for (auto const& path:_paths) {
        for (std::size_t i=1;i<path.points.size();++i)
            segment(path.points[i-1],path.points[i],_line_width,_line);
        if (path.closed && path.points.size()>1)
            segment(path.points.back(),path.points.front(),_line_width,_line);
    }
    _paths.clear();
}
Void EmbeddedCanvas::fill() {
    for (auto const& path:_paths) polygon(path,_fill);
    // Mirrors Cairo's cairo_fill_preserve() followed by stroke().
    stroke();
}
Void EmbeddedCanvas::finalise() {
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
Void EmbeddedCanvas::write(const char* name) const {
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
        throw std::runtime_error("Unable to write Embedded PNG: "+filename);
}
} // namespace Ariadne
