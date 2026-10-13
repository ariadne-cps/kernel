#ifndef ARIADNE_LABELLED_POINT_TPL_HPP
#define ARIADNE_LABELLED_POINT_TPL_HPP

#include "labelled/point.hpp"

namespace Ariadne {

template<class X> LabelledPoint<X>::LabelledPoint(RealSpace const& state_space, Point<X> const& pt)
    : Point<X>(pt), _state_variables(state_space.variable_names()) { }

template<class X> LabelledPoint<X>::LabelledPoint(List<Assignment<RealVariable,X>> const& x)
    : Point<X>(x.size(),[&x](SizeType i){return x[i].right_hand_side();}),
      _state_variables(RealSpace(left_hand_sides(x)).variable_names()) { }

template<class X> LabelledPoint<X>::LabelledPoint(InitializerList<Assignment<RealVariable,X>> const& x)
    : LabelledPoint<X>(List<Assignment<RealVariable,X>>(x)) { }

template<class X> LabelledPoint<X>* LabelledPoint<X>::labelled_clone() const {
    return new LabelledPoint(*this);
}

template<class X> RealSpace LabelledPoint<X>::state_space() const {
    return RealSpace(this->_state_variables);
}

template<class X> Void LabelledPoint<X>::draw(CanvasInterface& canvas, const Variables2d& axes) const {
    Projection2d proj=projection(this->state_space(),axes);
    Point<X>::draw(canvas,proj);
}

} // namespace Ariadne

#endif
