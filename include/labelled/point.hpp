#ifndef ARIADNE_LABELLED_POINT_HPP
#define ARIADNE_LABELLED_POINT_HPP

#include "geometry/point.hpp"
#include "symbolic/identifier.hpp"
#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "io/graphics_interface.hpp"
#include "io/figure.hpp"

namespace Ariadne {

template<class X> class LabelledPoint;

using LabelledRealPoint = LabelledPoint<Real>;
template<class IVL> class LabelledBox;
using LabelledExactBoxType = LabelledBox<ExactIntervalType>;

template<class X>
class LabelledPoint : public LabelledDrawable2dInterface, public Point<X> {
  public:
    LabelledPoint(RealSpace const& state_space, Point<X> const& pt);
    LabelledPoint(List<Assignment<RealVariable,X>> const& val);
    LabelledPoint(InitializerList<Assignment<RealVariable,X>> const& val);
    template<class T> requires Convertible<T,X>
    LabelledPoint(const LabelledPoint<T>& pt) : Point<X>(pt), _state_variables(pt.state_variables()) { }
    template<class Y, class... PRS> requires Constructible<X,Y,PRS...>
    LabelledPoint(const LabelledPoint<Y>& pt, PRS... prs) : Point<X>(pt,prs...), _state_variables(pt.state_variables()) { }

    LabelledPoint* labelled_clone() const override;
    List<Identifier> const& state_variables() const { return _state_variables; }
    RealSpace state_space() const;

    using Point<X>::draw;
    Void draw(CanvasInterface&, const Variables2d&) const override;

  private:
    List<Identifier> _state_variables;
};

template<class X> inline LabelledPoint<X> project(LabelledPoint<X> const& pt, Variables2d const& v) {
    Projection2d p = projection(pt.state_space(),v);
    return LabelledPoint<X>(pt.state_space(),project(static_cast<Point<X> const&>(pt),p));
}

} // namespace Ariadne

#endif
