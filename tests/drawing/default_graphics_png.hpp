/***************************************************************************
 *            tests/drawing/default_graphics_png.hpp
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

/*! \file tests/drawing/default_graphics_png.hpp
 *  \brief Minimal C++ reader for PNGs emitted by the default graphics backend.
 */

#ifndef ARIADNE_TEST_DEFAULT_GRAPHICS_PNG_HPP
#define ARIADNE_TEST_DEFAULT_GRAPHICS_PNG_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

//! \brief Helpers to verify rendered images in standalone C++ tests.
namespace DrawingTest {

//! \brief Decoded opaque RGB image for checking individual pixels.
struct Image {
    unsigned width=0,height=0;
    std::vector<std::array<unsigned char,3>> pixels;
    //! \brief Return one RGB pixel and validate its coordinates.
    std::array<unsigned char,3> pixel(unsigned x,unsigned y) const {
        if(x>=width||y>=height) throw std::runtime_error("PNG sample out of bounds");
        return pixels[static_cast<std::size_t>(y)*width+x];
    }
};

//! \brief Read a big-endian 32-bit integer from PNG bytes.
inline std::uint32_t big32(std::vector<unsigned char> const& b,std::size_t p) {
    if(p+4>b.size()) throw std::runtime_error("Truncated PNG");
    return (std::uint32_t(b[p])<<24)|(std::uint32_t(b[p+1])<<16)|
           (std::uint32_t(b[p+2])<<8)|b[p+3];
}

//! \brief Compute the PNG chunk CRC-32.
inline std::uint32_t crc32(std::vector<unsigned char> const& b,
                           std::size_t p,std::size_t n) {
    std::uint32_t crc=0xffffffffu;
    for(std::size_t i=p;i<p+n;++i) {
        crc^=b[i];
        for(int bit=0;bit<8;++bit)
            crc=(crc>>1)^((crc&1u)?0xedb88320u:0u);
    }
    return ~crc;
}

// Read the self-contained Default graphics backend's stored-block PNG.
// The optional comparison backend has a separate decoder via its own API.
//! \brief Decode the stored-block RGBA8 PNG emitted by Default graphics.
inline Image read_default_graphics_png(std::string const& name) {
    std::ifstream file(name,std::ios::binary);
    if(!file) throw std::runtime_error("Cannot open "+name);
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(file)),
                                      std::istreambuf_iterator<char>());
    const std::array<unsigned char,8> signature{137,80,78,71,13,10,26,10};
    if(bytes.size()<8 || !std::equal(signature.begin(),signature.end(),bytes.begin()))
        throw std::runtime_error("Invalid PNG signature");
    Image image;
    std::vector<unsigned char> compressed;
    for(std::size_t pos=8;pos<bytes.size();) {
        const auto len=big32(bytes,pos);
        if(len>bytes.size()-pos || len+12u>bytes.size()-pos)
            throw std::runtime_error("Invalid PNG chunk size");
        const std::size_t data=pos+8;
        if(big32(bytes,data+len)!=crc32(bytes,pos+4,4+len))
            throw std::runtime_error("PNG CRC mismatch");
        const std::string type(bytes.begin()+static_cast<std::ptrdiff_t>(pos+4),
                               bytes.begin()+static_cast<std::ptrdiff_t>(pos+8));
        if(type=="IHDR") {
            if(len!=13 || bytes[data+8]!=8 || bytes[data+9]!=6 ||
               bytes[data+10]!=0 || bytes[data+11]!=0 || bytes[data+12]!=0)
                throw std::runtime_error("Unsupported PNG format");
            image.width=big32(bytes,data);
            image.height=big32(bytes,data+4);
        } else if(type=="IDAT") {
            compressed.insert(compressed.end(),bytes.begin()+static_cast<std::ptrdiff_t>(data),
                              bytes.begin()+static_cast<std::ptrdiff_t>(data+len));
        }
        pos=data+len+4;
        if(type=="IEND") break;
    }
    if(image.width==0||image.height==0 || compressed.size()<6 ||
       compressed[0]!=0x78)
        throw std::runtime_error("Invalid default backend PNG");
    const std::size_t stride=static_cast<std::size_t>(image.width)*4+1;
    if(static_cast<std::uint64_t>(stride)*image.height>500000000ULL)
        throw std::runtime_error("PNG too large");
    std::vector<unsigned char> raw;
    raw.reserve(stride*image.height);
    std::size_t pos=2;
    bool done=false;
    while(!done) {
        if(pos+5>compressed.size()) throw std::runtime_error("Truncated stored block");
        const unsigned type=compressed[pos++];
        if((type&0xfeu)!=0) throw std::runtime_error("Unsupported DEFLATE block");
        done=(type&1u)!=0;
        const unsigned len=unsigned(compressed[pos])|(unsigned(compressed[pos+1])<<8);
        const unsigned complement=unsigned(compressed[pos+2])|
                                  (unsigned(compressed[pos+3])<<8);
        pos+=4;
        if((len^complement)!=65535u || len>compressed.size()-pos)
            throw std::runtime_error("Invalid DEFLATE length");
        raw.insert(raw.end(),compressed.begin()+static_cast<std::ptrdiff_t>(pos),
                   compressed.begin()+static_cast<std::ptrdiff_t>(pos+len));
        pos+=len;
    }
    if(raw.size()!=stride*image.height) throw std::runtime_error("PNG raster size mismatch");
    image.pixels.reserve(static_cast<std::size_t>(image.width)*image.height);
    for(unsigned y=0;y<image.height;++y) {
        const std::size_t row=static_cast<std::size_t>(y)*stride;
        if(raw[row]!=0) throw std::runtime_error("Unsupported PNG filter");
        for(unsigned x=0;x<image.width;++x) {
            const std::size_t offset=row+1+static_cast<std::size_t>(x)*4;
            if(raw[offset+3]!=255) throw std::runtime_error("Unexpected PNG alpha");
            image.pixels.push_back({raw[offset],raw[offset+1],raw[offset+2]});
        }
    }
    return image;
}

//! \brief Sample the image at a point in normalized plot coordinates.
inline std::array<unsigned char,3> world_pixel(Image const& image,
                                                 double x,double y) {
    const double width=static_cast<double>(image.width)-170;
    const double height=static_cast<double>(image.height)-50;
    const unsigned px=static_cast<unsigned>(std::lround(160+(x+1)*width/2));
    const unsigned py=static_cast<unsigned>(std::lround(10+(1-y)*height/2));
    return image.pixel(px,py);
}

//! \brief Report a failed graphics regression assertion.
inline void require(bool condition,std::string const& message) {
    if(!condition) throw std::runtime_error(message);
}

} // namespace DrawingTest
#endif
