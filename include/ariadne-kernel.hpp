#ifndef ARIADNE_ARIADNE_KERNEL_HPP
#define ARIADNE_ARIADNE_KERNEL_HPP

#include "ariadne-algebra.hpp"

#include "function/function.hpp"
#include "function/procedure.hpp"
#include "function/formula.hpp"
#include "function/constraint.hpp"
#include "function/function_model.hpp"
#include "function/function_patch.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "function/symbolic_function.hpp"

#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "geometry/set.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid.hpp"
#include "geometry/grid_paving.hpp"
#include "geometry/affine_set.hpp"

#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/assignment.hpp"

#include "solving/solver.hpp"
#include "solving/constraint_solver.hpp"
#include "solving/integrator.hpp"
#include "solving/linear_programming.hpp"
#include "solving/nonlinear_programming.hpp"

#include "io/figure.hpp"
#include "io/graphics_manager.hpp"
#include "io/command_line_interface.hpp"

#endif
