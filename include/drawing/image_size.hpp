/*
 * Ariadne - free software under the GNU General Public License version 3
 * or (at your option) any later version. No warranty.
 */
#ifndef ARIADNE_IMAGE_SIZE_HPP
#define ARIADNE_IMAGE_SIZE_HPP

#include "drawing/figure.hpp"

namespace Ariadne {

struct ImageSize2d {
    Nat nx, ny;
    ImageSize2d(Nat x, Nat y) : nx(x), ny(y) { }
    ImageSize2d(Int x, Int y) {
        ARIADNE_PRECONDITION(x > 0 && y > 0);
        nx = static_cast<Nat>(x);
        ny = static_cast<Nat>(y);
    }
};

} // namespace Ariadne
#endif // ARIADNE_IMAGE_SIZE_HPP
