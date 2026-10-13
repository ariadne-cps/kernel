#include "pybind11.hpp"

void algebra_module(pybind11::module& module);
void bdd_module(pybind11::module& module);

void function_submodule(pybind11::module& module);
void geometry_submodule(pybind11::module& module);
void io_submodule(pybind11::module& module);
void solving_submodule(pybind11::module& module);
void symbolic_submodule(pybind11::module& module);

void kernel_module(pybind11::module& module) {
    algebra_module(module);
    bdd_module(module);

    function_submodule(module);
    geometry_submodule(module);
    io_submodule(module);
    solving_submodule(module);
    symbolic_submodule(module);
}
