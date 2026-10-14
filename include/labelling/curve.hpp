#ifndef ARIADNE_LABELLING_CURVE_HPP
#define ARIADNE_LABELLING_CURVE_HPP

#include "geometry/curve.hpp"
#include "labelling/point.hpp"

namespace Ariadne {

class LabelledInterpolatedCurve
    : public LabelledDrawable2dInterface,
      public InterpolatedCurve,
      public Labelling
{
  public:
    typedef LabelledPoint<FloatDPApproximation> PointType;

    explicit LabelledInterpolatedCurve(const PointType& pt);

    LabelledInterpolatedCurve* labelled_clone() const override;

    using InterpolatedCurve::draw;
    Void draw(CanvasInterface&, const Variables2d&) const override;
};

} // namespace Ariadne

#endif
