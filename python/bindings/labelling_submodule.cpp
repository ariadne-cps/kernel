#include "pybind11.hpp"

void export_labelling(pybind11::module& module);

void labelling_submodule(pybind11::module& module) {
    export_labelling(module);
}
