#include "labelled/point.tpl.hpp"

namespace Ariadne {

template<> String class_name<LabelledPoint<Real>>() { return "LabelledPoint<Real>"; }
template<> String class_name<LabelledPoint<FloatDPApproximation>>() { return "LabelledPoint<FloatDPApproximation>"; }

template class LabelledPoint<Real>;
template class LabelledPoint<FloatDPApproximation>;

} // namespace Ariadne
