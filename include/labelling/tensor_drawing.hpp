/***************************************************************************
 *            labelling/tensor_drawing.hpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file labelling/tensor_drawing.hpp
 *  \brief Drawing adapters for algebraic Tensor data.
 */

#ifndef ARIADNE_LABELLING_TENSOR_DRAWING_HPP
#define ARIADNE_LABELLING_TENSOR_DRAWING_HPP


#include "algebra/tensor.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "labelling/figure.hpp"
#include "labelling/graphics_interface.hpp"

namespace Ariadne {

template<SizeType N, class PR>
class TensorDrawable
    : public LabelledDrawable2dInterface
    , public Drawable2dInterface
{
    static_assert(N == 2u, "TensorDrawable supports rank-2 tensors only.");

    using X = Float<PR>;
    Tensor<N,X> _tensor;

  public:
    explicit TensorDrawable(Tensor<N,X> const& tensor) : _tensor(tensor) { }

    TensorDrawable* clone() const override { return new TensorDrawable(*this); }
    TensorDrawable* labelled_clone() const override { return new TensorDrawable(*this); }

    DimensionType dimension() const override { return _tensor.rank(); }

    Void draw(CanvasInterface& canvas, Projection2d const&) const override {
        this->_draw_rank2(canvas);
    }

    Void draw(CanvasInterface& canvas, Variables2d const&) const override {
        this->_draw_rank2(canvas);
    }

  private:
    Void _draw_rank2(CanvasInterface& canvas) const {
        for (SizeType frame=0; frame!=_tensor.size(1); ++frame) {
            canvas.move_to(0.0, _tensor[{0,frame}].get_d());
            for (SizeType i=1; i!=_tensor.size(0); ++i) {
                canvas.line_to(numeric_cast<double>(i), _tensor[{i,frame}].get_d());
            }
            canvas.fill();
        }
    }

};

template<SizeType N, class PR>
TensorDrawable<N,PR> tensor_drawable(Tensor<N,Float<PR>> const& tensor) {
    return TensorDrawable<N,PR>(tensor);
}

} // namespace Ariadne

#endif // ARIADNE_LABELLING_TENSOR_DRAWING_HPP
