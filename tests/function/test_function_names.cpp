/***************************************************************************
 *            test_function_names.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/symbolic_function.hpp"
#include "macro/testing.hpp"
using namespace Ariadne;

// Diagnostics identify a function's paradigm, signature and implementation family.
template<class T> Void check_name(const String& expected) {
    ARIADNE_TEST_EQUALS(class_name<T>(),expected);
}
template<class P> Void check_paradigm_names(const String& prefix) {
    using N=Number<P>;
    check_name<ScalarUnivariateFunction<P>>(prefix+"ScalarUnivariateFunction");
    check_name<VectorUnivariateFunction<P>>(prefix+"VectorUnivariateFunction");
    check_name<ScalarMultivariateFunction<P>>(prefix+"ScalarMultivariateFunction");
    check_name<VectorMultivariateFunction<P>>(prefix+"VectorMultivariateFunction");
    check_name<ScalarUnivariateFormulaFunction<N>>(prefix+"ScalarUnivariateFormulaFunction");
    check_name<VectorUnivariateFormulaFunction<N>>(prefix+"VectorUnivariateFormulaFunction");
    check_name<ScalarFormulaFunction<N>>(prefix+"ScalarMultivariateFormulaFunction");
    check_name<VectorFormulaFunction<N>>(prefix+"VectorMultivariateFormulaFunction");
    check_name<ConstantFunction<N,RealScalar>>(prefix+"UnivariateConstantFunction");
    check_name<ConstantFunction<N,RealVector>>(prefix+"MultivariateConstantFunction");
    check_name<CoordinateFunction<P,RealScalar>>(prefix+"UnivariateCoordinateFunction");
    check_name<CoordinateFunction<P,RealVector>>(prefix+"MultivariateCoordinateFunction");
    check_name<UnaryFunction<P,RealScalar>>(prefix+"UnivariateUnaryFunction");
    check_name<UnaryFunction<P,RealVector>>(prefix+"MultivariateUnaryFunction");
    check_name<BinaryFunction<P,RealScalar>>(prefix+"UnivariateBinaryFunction");
    check_name<BinaryFunction<P,RealVector>>(prefix+"MultivariateBinaryFunction");
    check_name<GradedFunction<P,RealScalar>>(prefix+"UnivariateGradedFunction");
    check_name<GradedFunction<P,RealVector>>(prefix+"MultivariateGradedFunction");
    check_name<VectorOfScalarFunction<P,RealScalar>>(prefix+"VectorOfScalarUnivariateFunction");
    check_name<VectorOfScalarFunction<P,RealVector>>(prefix+"VectorOfScalarMultivariateFunction");
    check_name<ComposedFunction<P,RealScalar,RealScalar,RealScalar>>(prefix+"ScalarUnivariateComposedFunction");
    check_name<ComposedFunction<P,RealScalar,RealVector,RealScalar>>(prefix+"ScalarUnivariateComposedFunction");
    check_name<ComposedFunction<P,RealVector,RealScalar,RealScalar>>(prefix+"VectorUnivariateComposedFunction");
    check_name<ComposedFunction<P,RealVector,RealVector,RealScalar>>(prefix+"VectorUnivariateComposedFunction");
    check_name<ComposedFunction<P,RealScalar,RealScalar,RealVector>>(prefix+"ScalarMultivariateComposedFunction");
    check_name<ComposedFunction<P,RealScalar,RealVector,RealVector>>(prefix+"ScalarMultivariateComposedFunction");
    check_name<ComposedFunction<P,RealVector,RealScalar,RealVector>>(prefix+"VectorMultivariateComposedFunction");
    check_name<ComposedFunction<P,RealVector,RealVector,RealVector>>(prefix+"VectorMultivariateComposedFunction");
}
Int main() {
    ARIADNE_TEST_CALL(check_paradigm_names<ApproximateTag>("Approximate"));
    ARIADNE_TEST_CALL(check_paradigm_names<ValidatedTag>("Validated"));
    ARIADNE_TEST_CALL(check_paradigm_names<EffectiveTag>("Effective"));
    check_name<JoinedFunction<ValidatedTag,EuclideanDomain,EuclideanDomain,EuclideanDomain>>("ValidatedJoinedFunction");
    check_name<JoinedFunction<EffectiveTag,EuclideanDomain,EuclideanDomain,EuclideanDomain>>("EffectiveJoinedFunction");
    check_name<EmbeddedFunction<EffectiveTag,EuclideanDomain,EuclideanDomain,EuclideanDomain,RealDomain>>("EffectiveEmbeddedFunction");
    check_name<EmbeddedFunction<EffectiveTag,EuclideanDomain,EuclideanDomain,EuclideanDomain,EuclideanDomain>>("EffectiveEmbeddedFunction");
    check_name<LieDerivativeFunction<ValidatedTag>>("ValidatedLieDerivativeFunction");
    check_name<LieDerivativeFunction<EffectiveTag>>("EffectiveLieDerivativeFunction");
    return ARIADNE_TEST_FAILURES;
}
