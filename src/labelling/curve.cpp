#include "labelling/curve.hpp"
#include "labelling/figure.hpp"

namespace Ariadne {

LabelledInterpolatedCurve::LabelledInterpolatedCurve(const PointType& pt)
    : InterpolatedCurve(CurveInterface::PointType(pt)), Labelling(pt.space()) { }

LabelledInterpolatedCurve* LabelledInterpolatedCurve::labelled_clone() const {
    return new LabelledInterpolatedCurve(*this);
}

Void LabelledInterpolatedCurve::draw(CanvasInterface& c, const Variables2d& v) const {
    Projection2d proj=projection(this->space(),v);
    InterpolatedCurve::draw(c,proj);
}

} // namespace Ariadne
