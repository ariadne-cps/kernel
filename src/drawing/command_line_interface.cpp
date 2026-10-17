/***************************************************************************
 *            drawing/command_line_interface.cpp
 *
 *  Copyright  2007-26  Luca Geretti
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

#include "logging/command_line.hpp"
#include "threading/command_line.hpp"
#include "threading/thread_manager.hpp"

#include "drawer.hpp"
#include "graphics_manager.hpp"
#include "default_graphics.hpp"
#include "null_graphics.hpp"
#include "command_line_interface.hpp"

namespace Ariadne {

namespace {

ArgumentParser make_drawer_argument_parser() {
    return make_valued_argument_parser(
        "d",
        "drawer",
        "Choose the drawer as a <value> in [ box | grid@X | affine@X | none ] with X a non-negative integer, where none disables drawing (default: affine@0)",
        [](String const& value) -> VoidFunction {
            Drawer drawer=NullDrawer();

            if (value != "none") {
                if (value == "box") {
                    drawer=BoxDrawer();
                } else {
                    SizeType const at_pos=value.find('@',0);
                    if (at_pos == String::npos) {
                        throw std::runtime_error("invalid drawer");
                    }

                    Int const accuracy=std::stoi(value.substr(at_pos+1));
                    if (accuracy < 0) {
                        throw std::runtime_error("negative drawer accuracy");
                    }

                    String const drawer_name=value.substr(0,at_pos);
                    if (drawer_name == "grid") {
                        drawer=GridDrawer(static_cast<Nat>(accuracy));
                    } else if (drawer_name == "affine") {
                        drawer=AffineDrawer(static_cast<Nat>(accuracy));
                    } else {
                        throw std::runtime_error("invalid drawer");
                    }
                }
            }

            return [drawer] {
                GraphicsManager::instance().set_drawer(drawer);
            };
        });
}

ArgumentParser make_graphics_argument_parser() {
    return make_valued_argument_parser(
        "g",
        "graphics",
        "Choose the graphics backend as a <value> in [ default | none ] (default: default)",
        [](String const& value) -> VoidFunction {
            if (value == "default") {
                return [] { GraphicsManager::instance().set_backend(DefaultGraphicsBackend()); };
            }
            if (value == "none") {
                return [] {
                    GraphicsManager::instance().set_backend(
                        NullGraphicsBackend());
                };
            }
            throw std::runtime_error("invalid graphics backend");
        });
}

} // namespace

CommandLineContribution drawing_command_line() {
    return CommandLineContribution(
        "kernel-drawing",
        {
            make_drawer_argument_parser(),
            make_graphics_argument_parser()
        });
}

CommandLineInterface make_kernel_command_line_interface() {
    // ThreadManager owns the thread registry used by the logging schedulers.
    ThreadManager::instance();

    CommandLineInterface command_line;
    command_line
        .add(logging_command_line())
        .add(threading_command_line())
        .add(drawing_command_line());
    return command_line;
}

CommandLineInterface& kernel_command_line_interface() {
    static CommandLineInterface command_line=make_kernel_command_line_interface();
    return command_line;
}

} // namespace Ariadne
