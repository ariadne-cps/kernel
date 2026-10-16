/***************************************************************************
 *            drawing/graphics_manager.cpp
 *
 *  Copyright  2008-21  Mirko Albanese, Luca Geretti
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

#include "utility/handle.hpp"
#include "drawing/graphics_backend_interface.hpp"
#include "drawing/graphics_manager.hpp"
#include "drawing/figure.hpp"
#include "drawing/geometry2d.hpp"
#include "drawing/default_graphics.hpp"
#include "drawing/drawer.hpp"
#include "drawing/null_graphics.hpp"

namespace Ariadne {

GraphicsBackend default_graphics() {
    return DefaultGraphicsBackend();
}

GraphicsManager::GraphicsManager() : _backend(default_graphics()), _drawer(AffineDrawer(0)) {
}

GraphicsManager& GraphicsManager::instance() {
    static GraphicsManager instance;
    return instance;
}

GraphicsBackend const& GraphicsManager::backend() const {
    return _backend;
}

Void GraphicsManager::set_backend(GraphicsBackend const& backend) {
    _backend = backend;
}

Drawer const& GraphicsManager::drawer() const {
    return _drawer;
}

Void GraphicsManager::set_drawer(Drawer const& drawer) {
    _drawer = drawer;
}

} // namespace Ariadne


