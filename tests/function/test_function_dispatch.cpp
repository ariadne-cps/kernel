/***************************************************************************
 *            test_function_dispatch.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"
#include <array>
using namespace Ariadne;
using X=FloatDPApproximation;

template<class P> Void test_formula_constructors_and_dispatch() {
    using Y=Number<P>;using C=FunctionConstructors<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    auto t=Formula<Y>::coordinate(0u),y=Formula<Y>::coordinate(1u);
    U scalar(RealDomain(),sqr(t)+Y(1));
    W curve(RealDomain(),Vector<Formula<Y>>{t,t+Y(2)});
    S sum(EuclideanDomain(2u),t+y);
    V identity(EuclideanDomain(2u),Vector<Formula<Y>>{t,y});
    X parameter(3,dp);Vector<X> point({3,1},dp);
    ARIADNE_TEST_EQUALS(scalar(parameter),X(10,dp));
    ARIADNE_TEST_EQUALS(curve(parameter),Vector<X>({3,5},dp));
    ARIADNE_TEST_EQUALS(sum(point),X(4,dp));
    ARIADNE_TEST_EQUALS(identity(point),point);
    auto coordinate=C::coordinate();
    using UA=AlgebraOperations<U,Y>;using SA=AlgebraOperations<S,Y>;
    ARIADNE_TEST_EQUALS(UA::apply(BinaryElementaryOperator(Sub()),Y(5),coordinate)(parameter),X(2,dp));
    ARIADNE_TEST_EQUALS(UA::apply(GradedElementaryOperator(Pow()),coordinate,3)(parameter),X(27,dp));
    std::array<S,2> first{C::coordinate(2u,0u),C::coordinate(EuclideanDomain(2u),0u)};
    std::array<S,2> second{C::coordinate(2u,1u),C::coordinate(EuclideanDomain(2u),1u)};
    for(const auto& a:first) {
        ARIADNE_TEST_EQUALS((-a)(point),X(-3,dp));
        ARIADNE_TEST_EQUALS(SA::apply(BinaryElementaryOperator(Sub()),Y(5),a)(point),X(2,dp));
        ARIADNE_TEST_EQUALS(SA::apply(BinaryElementaryOperator(Sub()),a,Y(5))(point),X(-2,dp));
        V single{a};
        ARIADNE_TEST_EQUALS(single(point)[0],X(3,dp));
        for(const auto& b:second) {
            ARIADNE_TEST_EQUALS((a+b)(point),X(4,dp));
            ARIADNE_TEST_EQUALS((a-b)(point),X(2,dp));
            ARIADNE_TEST_EQUALS((a*b)(point),X(3,dp));
            ARIADNE_TEST_EQUALS((a/b)(point),X(3,dp));
        }
    }
    ARIADNE_TEST_FAIL(SA::apply(BinaryElementaryOperator(Add()),C::coordinate(1u,0u),first[0]));
}

Void test_model_composition() {
    using S=ValidatedScalarMultivariateFunction;using V=ValidatedVectorMultivariateFunction;
    auto t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedVectorUnivariateFunction curve{t,ValidatedNumber(2)*t};
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    V vector{x+y,x-y};auto scalar=ValidatedNumber(2)*x-y+ValidatedNumber(1);
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    S scalar_model=factory.create_constant(domain,ValidatedNumber(2));
    V vector_model=factory.create_identity(domain);
    ARIADNE_TEST_FAIL(compose(S::coordinate(3u,0u),vector_model));
    Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    ARIADNE_TEST_ASSERT(definitely(compose(t+ValidatedNumber(1),scalar_model)(point)==FloatDPBounds(3,dp)));
    ARIADNE_TEST_EQUALS(compose(curve,scalar_model)(point),Vector<FloatDPBounds>({2,4},dp));
    ARIADNE_TEST_ASSERT(definitely(compose(scalar,vector_model)(point)==FloatDPBounds(1,dp)));
    ARIADNE_TEST_EQUALS(compose(vector,vector_model)(point),Vector<FloatDPBounds>({0.75_x,-0.25_x},dp));
}

Void test_interval_ranges() {
    using S=ValidatedScalarMultivariateFunction;using V=ValidatedVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    auto scalar=ValidatedNumber(2)*x-y+ValidatedNumber(1);V vector{x+y,x-y};
    Vector<UpperIntervalType> input{UpperIntervalType(-1,1),UpperIntervalType(2,3)};
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(evaluate_range(scalar,input)),ExactIntervalType(-4,1)));
    auto ranges=evaluate_range(vector,input);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(ranges[0]),ExactIntervalType(1,4)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(ranges[1]),ExactIntervalType(-4,-1)));
    auto gradient=gradient_range(scalar,input);auto jacobian=jacobian_range(vector,input);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(gradient[0]),ExactIntervalType(2,2)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(gradient[1]),ExactIntervalType(-1,-1)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(jacobian[0][0]),ExactIntervalType(1,1)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(jacobian[1][1]),ExactIntervalType(-1,-1)));
    using D=Differential<UpperIntervalType>;
    auto variables=D::identity(1u,input);
    auto derivatives=derivative_range(vector,variables);
    ARIADNE_TEST_EQUALS(derivatives.size(),2u);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(derivatives[0].value()),ExactIntervalType(1,4)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(derivatives[1].gradient()[1]),ExactIntervalType(-1,-1)));
    ARIADNE_TEST_EQUALS(derivative_range(V::zeros(0u,2u),variables).size(),0u);
    ARIADNE_TEST_EQUALS(evaluate_range(V::zeros(0u,2u),input).size(),0u);
    Vector<UpperIntervalType> empty(0u,UpperIntervalType(0,0));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(evaluate_range(S::constant(0u,5),empty)),ExactIntervalType(5,5)));
    ARIADNE_TEST_EQUALS(gradient_range(S::constant(0u,5),empty).size(),0u);
    ARIADNE_TEST_EQUALS(jacobian_range(V::zeros(0u,0u),empty).row_size(),0u);
    Vector<D> no_variables(0u,D(0u,1u,UpperIntervalType(0,0)));
    auto constant=derivative_range(V::constant(0u,Vector<ValidatedNumber>{2,3}),no_variables);
    ARIADNE_TEST_EQUALS(constant.size(),2u);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(constant[0].value()),ExactIntervalType(2,2)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(constant[1].value()),ExactIntervalType(3,3)));
    ARIADNE_TEST_EQUALS(constant[0].degree(),1u);
    Vector<D> unused_parameters(0u,D(3u,2u,UpperIntervalType(0,0)));
    auto lifted=derivative_range(V::constant(0u,Vector<ValidatedNumber>{2,3}),unused_parameters);
    ARIADNE_TEST_EQUALS(lifted[0].argument_size(),3u);
    ARIADNE_TEST_EQUALS(lifted[0].degree(),2u);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(lifted[0].gradient()[0]),ExactIntervalType(0,0)));
    Vector<UpperIntervalType> wrong_dimension(1u,UpperIntervalType(0,0));
    Vector<D> wrong_differentials(1u,D(2u,1u,UpperIntervalType(0,0)));
    ARIADNE_TEST_FAIL(evaluate_range(scalar,wrong_dimension));
    ARIADNE_TEST_FAIL(evaluate_range(vector,wrong_dimension));
    ARIADNE_TEST_FAIL(gradient_range(scalar,wrong_dimension));
    ARIADNE_TEST_FAIL(jacobian_range(vector,wrong_dimension));
    ARIADNE_TEST_FAIL(derivative_range(vector,wrong_differentials));
}
template<class P> Void test_vector_contracts() {
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    auto s1=S::coordinate(1u,0u),s2=S::coordinate(2u,0u);
    auto vector=V::identity(2u),short_result=V::zeros(1u,2u),long_result=V::zeros(3u,2u),wrong_arguments=V::zeros(2u,3u);
    if constexpr (!Same<P,ApproximateTag>) {
        ARIADNE_TEST_FAIL(vector+short_result);
        ARIADNE_TEST_FAIL(vector+long_result);
        ARIADNE_TEST_FAIL(vector+wrong_arguments);
        ARIADNE_TEST_FAIL(vector-short_result);
        ARIADNE_TEST_FAIL(vector-long_result);
        ARIADNE_TEST_FAIL(vector-wrong_arguments);
        ARIADNE_TEST_FAIL(s1*vector);
        ARIADNE_TEST_FAIL(vector*s1);
        ARIADNE_TEST_FAIL(vector/s1);
        auto t=U::coordinate();W curve{t,t},short_curve{t},long_curve{t,t,t};
        ARIADNE_TEST_FAIL(curve+short_curve);
        ARIADNE_TEST_FAIL(curve+long_curve);
    }
    ARIADNE_TEST_FAIL(join(s1,s2));
    ARIADNE_TEST_FAIL(join(vector,s1));
    ARIADNE_TEST_FAIL(join(s1,vector));
    ARIADNE_TEST_FAIL(join(vector,wrong_arguments));
    ARIADNE_TEST_FAIL(compose(s2,long_result));
}
Void test_effective_contracts_and_derivative_failures() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto vector=V::identity(2u);Vector<EffectiveNumber> short_constant{1};
    ARIADNE_TEST_FAIL(short_constant+vector);
    ARIADNE_TEST_FAIL(vector+short_constant);
    ARIADNE_TEST_FAIL(short_constant-vector);
    ARIADNE_TEST_FAIL(vector-short_constant);
    auto s3=S::coordinate(3u,0u),s2=S::coordinate(2u,0u);
    auto non_square=V::zeros(2u,3u),wrong_g=V::zeros(2u,3u);
    ARIADNE_TEST_FAIL(lie_derivative(s3,vector));
    ARIADNE_TEST_FAIL(lie_derivative(s2,non_square));
    ARIADNE_TEST_FAIL(lie_derivative(wrong_g,vector));
    ARIADNE_TEST_FAIL(lie_derivative(vector,non_square));
    auto x=Formula<EffectiveNumber>::coordinate(0u),y=Formula<EffectiveNumber>::coordinate(1u);
    S nonsmooth(EuclideanDomain(2u),max(x,y));
    ARIADNE_TEST_FAIL(lie_derivative(nonsmooth,vector));
    ARIADNE_TEST_FAIL(lie_derivative(V{nonsmooth,s2},vector));
}
Int main() {
    ARIADNE_TEST_CALL(test_formula_constructors_and_dispatch<ApproximateTag>());
    ARIADNE_TEST_CALL(test_formula_constructors_and_dispatch<ValidatedTag>());
    ARIADNE_TEST_CALL(test_formula_constructors_and_dispatch<EffectiveTag>());
    ARIADNE_TEST_CALL(test_vector_contracts<ApproximateTag>());
    ARIADNE_TEST_CALL(test_vector_contracts<ValidatedTag>());
    ARIADNE_TEST_CALL(test_vector_contracts<EffectiveTag>());
    ARIADNE_TEST_CALL(test_effective_contracts_and_derivative_failures());
    ARIADNE_TEST_CALL(test_model_composition());
    ARIADNE_TEST_CALL(test_interval_ranges());
    return ARIADNE_TEST_FAILURES;
}
