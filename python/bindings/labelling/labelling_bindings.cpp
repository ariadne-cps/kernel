#include "pybind11.hpp"
#include "kernel-utilities.hpp"

#include "geometry/box.hpp"
#include "labelling/figure.hpp"
#include "labelling/set.hpp"

using namespace Ariadne;

void export_labelling(pybind11::module& module) {
    static constexpr auto reference_internal = pybind11::return_value_policy::reference_internal;

    pybind11::class_<Variables2d> variables2d_class(module,"Variables2d");
    variables2d_class.def(pybind11::init<RealVariable,RealVariable>());
    variables2d_class.def("x", &Variables2d::x);
    variables2d_class.def("y", &Variables2d::y);

    module.def("projection",&projection);

    pybind11::class_<Axes2d> axes2d_class(module,"Axes2d");
    axes2d_class.def(pybind11::init<double,RealVariable,double,double,RealVariable,double>());
    axes2d_class.def(pybind11::init<Real,RealVariable,Real,Real,RealVariable,Real>());

    using LabelledRealBox = LabelledBox<RealInterval>;
    pybind11::class_<LabelledRealBox> labelled_real_box_class(module,"LabelledRealBox");
    labelled_real_box_class.def(pybind11::init<RealSpace,RealBox>());
    labelled_real_box_class.def("dimension", &LabelledRealBox::dimension);
    labelled_real_box_class.def("__str__", &__cstr__<LabelledRealBox>);
    labelled_real_box_class.def("__repr__", &__cstr__<LabelledRealBox>);

    using LabelledFloatDPUpperBox = LabelledBox<FloatDPUpperInterval>;
    pybind11::class_<LabelledFloatDPUpperBox> labelled_floatdp_upper_box_class(module,"LabelledFloatDPUpperBox");
    labelled_floatdp_upper_box_class.def("__getitem__", &LabelledFloatDPUpperBox::operator[]);
    labelled_floatdp_upper_box_class.def("__str__", &__cstr__<LabelledFloatDPUpperBox>);

    pybind11::class_<LabelledFigure> labelled_figure_class(module,"LabelledFigure");
    labelled_figure_class.def(pybind11::init<Axes2d>());
    labelled_figure_class.def("set_axes",(Void(LabelledFigure::*)(const Axes2d&)) &LabelledFigure::set_axes, reference_internal);
    labelled_figure_class.def("properties",(GraphicsProperties&(LabelledFigure::*)())&LabelledFigure::properties, reference_internal);
    labelled_figure_class.def("draw",(LabelledFigure&(LabelledFigure::*)(const LabelledDrawable2dInterface&))&LabelledFigure::draw, reference_internal);
    labelled_figure_class.def("draw", [](LabelledFigure& fig, LabelledRealBox const& box) -> LabelledFigure& {
        LabelledDrawableWrapper<RealBox> wrapper(box);
        return fig.draw(wrapper);
    }, reference_internal);
    labelled_figure_class.def("clear",&LabelledFigure::clear, reference_internal);
    labelled_figure_class.def("set_animated", &LabelledFigure::set_animated, reference_internal);
    labelled_figure_class.def("write",(Void(LabelledFigure::*)(const Char*)const)&LabelledFigure::write);
    labelled_figure_class.def("write",(Void(LabelledFigure::*)(const Char*,Nat,Nat)const)&LabelledFigure::write);
}
