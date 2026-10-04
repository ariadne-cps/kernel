#include "pybind11.hpp"

void foundation_submodule(pybind11::module& module);
void numeric_submodule(pybind11::module& module);
void interval_submodule(pybind11::module& module);
void algebra_submodule(pybind11::module& module);
void linear_algebra_submodule(pybind11::module& module);
void polynomial_submodule(pybind11::module& module);
void differentiation_submodule(pybind11::module& module);

void function_submodule(pybind11::module& module);
void calculus_submodule(pybind11::module& module);
void geometry_submodule(pybind11::module& module);
void io_submodule(pybind11::module& module);
void solver_submodule(pybind11::module& module);
void optimization_submodule(pybind11::module& module);
void storage_submodule(pybind11::module& module);
void symbolic_submodule(pybind11::module& module);
void graphics_submodule(pybind11::module& module);

PYBIND11_MODULE(pyariadne, module) {
    foundation_submodule(module);
    numeric_submodule(module);
    interval_submodule(module);
    algebra_submodule(module);
    linear_algebra_submodule(module);
    polynomial_submodule(module);
    differentiation_submodule(module);

    function_submodule(module);
    calculus_submodule(module);
    geometry_submodule(module);
    io_submodule(module);
    solver_submodule(module);
    optimization_submodule(module);
    storage_submodule(module);
    symbolic_submodule(module);
    graphics_submodule(module);
}
