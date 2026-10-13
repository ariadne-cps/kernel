#include "labelling/curve.hpp"
#include "labelling/figure.hpp"

namespace Ariadne {

LabelledInterpolatedCurve::LabelledInterpolatedCurve(const PointType& pt)
    : InterpolatedCurve(CurveInterface::PointType(pt)), _state_variables(pt.state_variables()) { }

RealSpace LabelledInterpolatedCurve::state_space() const {
    return RealSpace(this->_state_variables);
}

LabelledInterpolatedCurve* LabelledInterpolatedCurve::labelling_clone() const {
    return new LabelledInterpolatedCurve(*this);
}

Void LabelledInterpolatedCurve::draw(CanvasInterface& c, const Variables2d& v) const {
    Projection2d proj=projection(this->state_space(),v);
    InterpolatedCurve::draw(c,proj);
}

} // namespace Ariadne
