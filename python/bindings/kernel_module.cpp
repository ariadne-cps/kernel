#include "pybind11.hpp"

void function_submodule(pybind11::module& module);
void geometry_submodule(pybind11::module& module);
void drawing_submodule(pybind11::module& module);
void solving_submodule(pybind11::module& module);

void kernel_module(pybind11::module& module) {
    function_submodule(module);
    geometry_submodule(module);
    drawing_submodule(module);
    solving_submodule(module);
}
