/***************************************************************************
 *            function/taylor_function.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "function/functional.hpp"
#include <iostream>
#include <iomanip>

#include "utility/macros.hpp"
#include "utility/exceptions.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/polynomial.hpp"
#include "algebra/differential.hpp"
#include "algebra/evaluate.hpp"
#include "function/taylor_model.hpp"

#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/function_patch.hpp"
#include "function/scaled_function_patch.hpp"

#include "function/taylor_function.hpp"

#include "taylor_model.tpl.hpp"
#include "scaled_function_patch.tpl.hpp"
#include "function_mixin.tpl.hpp"

#define VOLATILE ;

namespace Ariadne {

// FunctionMixin::_call is defined out-of-line. Instantiate only those virtual
// members required by the exported patch vtables. Instantiating the whole class
// also materializes unrelated virtuals and triggers premature instantiation.
#define ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(F, ...) \
    template auto FunctionMixin<F,ValidatedTag,RealScalar(RealVector)>:: \
    _call(Vector<__VA_ARGS__> const&) const -> Scalar<__VA_ARGS__>;
#define ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(F, ...) \
    template auto FunctionMixin<F,ValidatedTag,RealVector(RealVector)>:: \
    _call(Vector<__VA_ARGS__> const&) const -> Vector<__VA_ARGS__>;

#define ARIADNE_INSTANTIATE_VALIDATED_FUNCTION_MIXIN_CALLS(M) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,ApproximateNumber) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,FloatDPApproximation) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,FloatMPApproximation) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,Differential<FloatDPApproximation>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,Differential<FloatMPApproximation>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ApproximateTag,FloatDP>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ApproximateTag,FloatMP>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,Formula<ApproximateNumber>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,ElementaryAlgebra<ApproximateNumber>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,ValidatedNumber) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,FloatDPBounds) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,FloatMPBounds) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,Differential<FloatDPBounds>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,Differential<FloatMPBounds>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatDP>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatMP>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatDPBounds>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatMPBounds>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatDPUpperInterval>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatMPUpperInterval>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,Formula<ValidatedNumber>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,ElementaryAlgebra<ValidatedNumber>) \
    ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL(ScaledFunctionPatch<M>,ValidatedScalarMultivariateFunction) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,ApproximateNumber) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,FloatDPApproximation) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,FloatMPApproximation) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,Differential<FloatDPApproximation>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,Differential<FloatMPApproximation>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ApproximateTag,FloatDP>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ApproximateTag,FloatMP>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,Formula<ApproximateNumber>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,ElementaryAlgebra<ApproximateNumber>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,ValidatedNumber) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,FloatDPBounds) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,FloatMPBounds) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,Differential<FloatDPBounds>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,Differential<FloatMPBounds>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatDP>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatMP>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatDPBounds>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatMPBounds>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatDPUpperInterval>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,TaylorModel<ValidatedTag,FloatMPUpperInterval>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,Formula<ValidatedNumber>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,ElementaryAlgebra<ValidatedNumber>) \
    ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL(VectorScaledFunctionPatch<M>,ValidatedScalarMultivariateFunction)

ARIADNE_INSTANTIATE_VALIDATED_FUNCTION_MIXIN_CALLS(ValidatedTaylorModelDP)
template class ScaledFunctionPatch<ValidatedTaylorModelDP>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelDP>;
template class ScaledFunctionPatchFactory<ValidatedTaylorModelDP>;

ARIADNE_INSTANTIATE_VALIDATED_FUNCTION_MIXIN_CALLS(ValidatedBoundsTaylorModelDP)
template class ScaledFunctionPatch<ValidatedBoundsTaylorModelDP>;
template class VectorScaledFunctionPatch<ValidatedBoundsTaylorModelDP>;
template class ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelDP>;

ARIADNE_INSTANTIATE_VALIDATED_FUNCTION_MIXIN_CALLS(ValidatedTaylorModelMP)
template class ScaledFunctionPatch<ValidatedTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelMP>;
template class ScaledFunctionPatchFactory<ValidatedTaylorModelMP>;

ARIADNE_INSTANTIATE_VALIDATED_FUNCTION_MIXIN_CALLS(ValidatedBoundsTaylorModelMP)
template class ScaledFunctionPatch<ValidatedBoundsTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedBoundsTaylorModelMP>;
template class ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelMP>;

#undef ARIADNE_INSTANTIATE_VALIDATED_FUNCTION_MIXIN_CALLS
#undef ARIADNE_INSTANTIATE_VECTOR_FUNCTION_MIXIN_CALL
#undef ARIADNE_INSTANTIATE_SCALAR_FUNCTION_MIXIN_CALL



FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory(Sweeper<FloatDP> const& sweeper) {
    return new TaylorFunctionFactory(sweeper);
}
FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory() {
    return make_taylor_function_factory(Sweeper<FloatDP>());
}
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory(Sweeper<FloatDP> const& sweeper) {
    return new TaylorFunctionFactory(sweeper);
}
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory() {
    return make_taylor_function_patch_factory(Sweeper<FloatDP>());
}

} // namespace Ariadne
