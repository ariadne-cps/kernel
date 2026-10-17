/***************************************************************************
 *            default_graphics_golden.hpp
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

/*! \file default_graphics_golden.hpp
 *  \brief Pixel-exact regression checks against versioned PNG golden images.
 */

#ifndef ARIADNE_TEST_DEFAULT_GRAPHICS_GOLDEN_HPP
#define ARIADNE_TEST_DEFAULT_GRAPHICS_GOLDEN_HPP

#include "default_graphics_png.hpp"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace DrawingTest {

//! \brief Check the entire RGB image against a reviewed PNG in golden/.
//! \details The fixture is always read-only. On failure a red-on-black PPM
//! difference image is written in the test's working directory; each changed
//! pixel is red with intensity proportional to the largest channel error.
inline void require_golden_image(Image const& current, std::string const& filename) {
    const std::string path = std::string(ARIADNE_DRAWING_GOLDEN_DIR) + "/" + filename;
    const Image reference = read_default_graphics_png(path);

    if(current.width != reference.width || current.height != reference.height) {
        std::ostringstream error;
        error << filename << ": image dimensions differ (actual "
              << current.width << 'x' << current.height << ", golden "
              << reference.width << 'x' << reference.height << ')';
        throw std::runtime_error(error.str());
    }

    std::size_t differences = 0;
    std::size_t first_index = 0;
    unsigned largest_error = 0;
    std::vector<unsigned char> diagnostic;
    diagnostic.reserve(current.pixels.size() * 3);

    for(std::size_t i = 0; i < current.pixels.size(); ++i) {
        const auto& actual = current.pixels[i];
        const auto& expected = reference.pixels[i];
        unsigned maximum = 0;
        for(std::size_t channel = 0; channel < 3; ++channel) {
            const int delta = static_cast<int>(actual[channel]) -
                              static_cast<int>(expected[channel]);
            maximum = std::max(maximum, static_cast<unsigned>(std::abs(delta)));
        }
        if(maximum != 0) {
            if(differences == 0) first_index = i;
            ++differences;
            largest_error = std::max(largest_error, maximum);
        }
        diagnostic.push_back(static_cast<unsigned char>(
            maximum == 0 ? 0u : std::max(64u, maximum)));
        diagnostic.push_back(0);
        diagnostic.push_back(0);
    }

    const std::string diff_path = filename + "-golden-diff.ppm";
    if(differences == 0) {
        // Avoid confusing an old diagnostic with a fresh successful run.
        (void)std::remove(diff_path.c_str());
        return;
    }

    std::ofstream out(diff_path, std::ios::binary);
    if(!out) throw std::runtime_error("Cannot write golden diagnostic: " + diff_path);
    out << "P6\n" << current.width << ' ' << current.height << "\n255\n";
    out.write(reinterpret_cast<const char*>(diagnostic.data()),
              static_cast<std::streamsize>(diagnostic.size()));
    if(!out) throw std::runtime_error("Incomplete golden diagnostic: " + diff_path);

    std::ostringstream error;
    error << filename << ": " << differences << " changed RGB pixels; first at ("
          << first_index % current.width << ", " << first_index / current.width
          << "), maximum channel difference " << largest_error
          << "; see " << diff_path;
    throw std::runtime_error(error.str());
}

} // namespace DrawingTest

#endif // ARIADNE_TEST_DEFAULT_GRAPHICS_GOLDEN_HPP
