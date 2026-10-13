/***************************************************************************
 *            io/command_line_interface.hpp
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

/*! \file io/command_line_interface.hpp
 *  \brief Kernel command-line composition.
 */

#ifndef ARIADNE_COMMAND_LINE_INTERFACE_HPP
#define ARIADNE_COMMAND_LINE_INTERFACE_HPP

#include "utility/command_line_interface.hpp"

namespace Ariadne {

//! \brief Command-line arguments owned by the kernel IO layer.
//! \details Contributes drawer and graphics backend options only.
CommandLineContribution io_command_line();

//! \brief Build the complete kernel command-line profile.
//! \details Composes logging, threading and kernel IO contributions.
CommandLineInterface make_kernel_command_line_interface();

//! \brief Shared kernel command-line profile used by applications and bindings.
CommandLineInterface& kernel_command_line_interface();

} // namespace Ariadne

#endif // ARIADNE_COMMAND_LINE_INTERFACE_HPP
