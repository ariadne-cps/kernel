/***************************************************************************
 *            test_affine_model.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/affine_model.hpp"
#include "function/affine.hpp"
#include "function/taylor_model.hpp"
#include "function/function.hpp"
#include "function/affine_model.tpl.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

template<class F> Void test_product(typename F::PrecisionType pr) {
    using A=AffineModel<ApproximateTag,F>;using X=Approximation<F>;
    A a(X(2,pr),{X(3,pr),X(-1,pr)}),b(X(5,pr),{X(7,pr),X(2,pr)});
    auto product=a*b;
    ARIADNE_TEST_EQUALS(product.value(),X(10,pr));
    ARIADNE_TEST_EQUALS(product.gradient(0u),X(29,pr));
    ARIADNE_TEST_EQUALS(product.gradient(1u),X(-1,pr));
    auto constant=A::constant(0u,X(2,pr))*A::constant(0u,X(5,pr));
    ARIADNE_TEST_EQUALS(constant.value(),X(10,pr));
    using V=AffineModel<ValidatedTag,F>;using E=typename V::ErrorType;
    V factor(1u,F(2,pr)),uncertain(1u,F(0,pr));uncertain.set_error(E(3u,pr));
    auto enclosure=factor*uncertain;
    ARIADNE_TEST_ASSERT(decide(enclosure.error().raw()>=F(6,pr)));
    ARIADNE_TEST_ASSERT(decide(enclosure.range().lower_bound()<=F(-6,pr)));
    ARIADNE_TEST_ASSERT(decide(enclosure.range().upper_bound()>=F(6,pr)));
    ARIADNE_TEST_ASSERT(decide((uncertain*factor).error().raw()>=F(6,pr)));
}

template<class F> Void test_approximate_model(typename F::PrecisionType pr) {
    using A=AffineModel<ApproximateTag,F>;using X=Approximation<F>;
    A a(X(2,pr),{X(3,pr),X(-1,pr),X(0,pr)});
    Vector<X> point({0.25_x,0.5_x,0.0_x},pr);
    ARIADNE_TEST_EQUALS(a.evaluate(point),X(2.25_x,pr));
    ARIADNE_TEST_EQUALS(a.range().lower_bound(),X(-2,pr));
    ARIADNE_TEST_EQUALS(a.range().upper_bound(),X(6,pr));
    ARIADNE_TEST_EQUALS(a.b(),a.value());
    ARIADNE_TEST_EQUALS(a.a(),a.gradient());
    ARIADNE_TEST_EQUALS(a.precision(),pr);
    ARIADNE_TEST_EQUALS(std::get<0>(a.characteristics()),3u);
    ARIADNE_TEST_EQUALS(a.create_zero().evaluate(point),X(0,pr));
    const A& ca=a;
    ARIADNE_TEST_EQUALS(ca[0u],X(3,pr));
    a[0u]=X(4,pr);a.set_value(X(3,pr));
    a.set_gradient(0u,X(3,pr));a.set_gradient(1u,ApproximateNumber(-1));
    ARIADNE_TEST_EQUALS(a.evaluate(point),X(3.25_x,pr));
    auto b=A::coordinate(3u,1u,pr);X coefficient(2,pr);
    ARIADNE_TEST_EQUALS((-a).evaluate(point),X(-3.25_x,pr));
    ARIADNE_TEST_EQUALS((a+b).evaluate(point),X(3.75_x,pr));
    ARIADNE_TEST_EQUALS((coefficient+a).evaluate(point),X(5.25_x,pr));
    ARIADNE_TEST_EQUALS((a+coefficient).evaluate(point),X(5.25_x,pr));
    ARIADNE_TEST_EQUALS((coefficient*a).evaluate(point),X(6.5_x,pr));
    ARIADNE_TEST_EQUALS((a*coefficient).evaluate(point),X(6.5_x,pr));
    A from_covector(X(3,pr),Covector<X>{X(3,pr),X(-1,pr),X(0,pr)});
    ARIADNE_TEST_EQUALS(from_covector.evaluate(point),a.evaluate(point));
    a=coefficient;
    ARIADNE_TEST_EQUALS(a.evaluate(point),coefficient);
    A constant(0u,coefficient);constant=X(5,pr);
    Vector<X> empty(0u,pr);
    ARIADNE_TEST_EQUALS(constant.evaluate(empty),X(5,pr));
    ARIADNE_TEST_EQUALS(constant.range().lower_bound(),X(5,pr));
    ARIADNE_TEST_FAIL(a+A(2u,pr));
    ARIADNE_TEST_FAIL(a*A(2u,pr));
    std::ostringstream text;text<<from_covector<<constant;
    ARIADNE_TEST_PRINT(text.str());
    ARIADNE_TEST_ASSERT(text.str().find("*x0")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("*x1")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("*x2")==String::npos);
    BoxDomainType domain{{2,6},{-4,2}};
    auto scalings=A::scalings(domain,pr);
    ARIADNE_TEST_EQUALS(scalings[0u].evaluate(Vector<X>({-1,1},pr)),X(2,pr));
    ARIADNE_TEST_EQUALS(scalings[1u].evaluate(Vector<X>({-1,1},pr)),X(2,pr));
    ARIADNE_TEST_EQUALS(A::scalings(BoxDomainType(0u),pr).size(),0u);
}
template<class F> Void test_validated_model(typename F::PrecisionType pr) {
    using A=AffineModel<ValidatedTag,F>;using X=Bounds<F>;using E=typename A::ErrorType;
    A a(F(2,pr),{F(3,pr),F(-1,pr),F(0,pr)});a.set_error(E(F(0.5_x,pr)));
    const A& ca=a;
    Vector<X> point({0.25_x,0.5_x,0.0_x},pr);
    auto value=a.evaluate(point);
    ARIADNE_TEST_SAME(value,X(1.75_x,2.75_x,pr));
    ARIADNE_TEST_ASSERT(decide(a.range().lower_bound().raw()==F(-2.5_x,pr)));
    ARIADNE_TEST_ASSERT(decide(a.range().upper_bound().raw()==F(6.5_x,pr)));
    ARIADNE_TEST_EQUALS(ca.b(),ca.value());ARIADNE_TEST_EQUALS(ca.a(),ca.gradient());
    ARIADNE_TEST_EQUALS(ca.e().raw(),ca.error().raw());
    ARIADNE_TEST_EQUALS(ca[0u],F(3,pr));ARIADNE_TEST_EQUALS(ca.gradient(1u),F(-1,pr));
    ARIADNE_TEST_EQUALS(std::get<0>(ca.characteristics()),3u);
    ARIADNE_TEST_EQUALS(ca.precision(),pr);
    ARIADNE_TEST_ASSERT(definitely(a.create_zero().evaluate(point)==X(0,pr)));
    a.set_value(F(3,pr));a.set_gradient(0u,F(3,pr));a.set_gradient(1u,Dyadic(-1));a[2u]=F(0,pr);
    A from_covector(F(3,pr),Covector<F>{F(3,pr),F(-1,pr),F(0,pr)},E(F(0.5_x,pr)));
    ARIADNE_TEST_SAME(from_covector.evaluate(point),a.evaluate(point));
    auto b=A::coordinate(3u,1u,pr);X coefficient(2,pr);
    ARIADNE_TEST_SAME((-a).evaluate(point),X(-3.75_x,-2.75_x,pr));
    ARIADNE_TEST_SAME((a+b).evaluate(point),X(3.25_x,4.25_x,pr));
    ARIADNE_TEST_SAME((coefficient+a).evaluate(point),X(4.75_x,5.75_x,pr));
    ARIADNE_TEST_SAME((a+coefficient).evaluate(point),X(4.75_x,5.75_x,pr));
    ARIADNE_TEST_SAME((coefficient*a).evaluate(point),X(5.5_x,7.5_x,pr));
    ARIADNE_TEST_SAME((a*coefficient).evaluate(point),X(5.5_x,7.5_x,pr));
    a=F(2,pr);ARIADNE_TEST_EQUALS(a.error().raw(),F(0,pr));
    ARIADNE_TEST_ASSERT(definitely(a.evaluate(point)==X(2,pr)));
    A constant=A::constant(0u,F(2,pr));constant=F(5,pr);
    Vector<X> empty(0u,pr);
    ARIADNE_TEST_ASSERT(definitely(constant.evaluate(empty)==X(5,pr)));
    ARIADNE_TEST_EQUALS((constant*constant).value(),F(25,pr));
    ARIADNE_TEST_EQUALS((constant*constant).error().raw(),F(0,pr));
    ARIADNE_TEST_FAIL(a+A(2u,pr));ARIADNE_TEST_FAIL(a*A(2u,pr));
    std::ostringstream text;text<<from_covector<<constant;
    ARIADNE_TEST_PRINT(text.str());
    ARIADNE_TEST_ASSERT(text.str().find("*x0")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("*x1")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("*x2")==String::npos);
    BoxDomainType domain{{2,6},{-4,2}};
    auto scalings=A::scalings(domain,pr);
    ARIADNE_TEST_ASSERT(definitely(scalings[0u].evaluate(Vector<X>({-1,1},pr))==X(2,pr)));
    ARIADNE_TEST_ASSERT(definitely(scalings[1u].evaluate(Vector<X>({-1,1},pr))==X(2,pr)));
    ARIADNE_TEST_EQUALS(A::scalings(BoxDomainType(0u),pr).size(),0u);
    Affine<X> uncertain(X(1,3,pr),{X(2,4,pr),X(-2,0,pr)});
    auto old_rounding=F::get_rounding_mode();
    auto converted=affine_model(uncertain);
    ARIADNE_TEST_ASSERT(F::get_rounding_mode()==old_rounding);
    ARIADNE_TEST_EQUALS(converted.value(),F(2,pr));
    ARIADNE_TEST_EQUALS(converted.gradient(0u),F(3,pr));
    ARIADNE_TEST_EQUALS(converted.error().raw(),F(3,pr));
    Sweeper<F> sweeper=ThresholdSweeper<F>(pr,1e-12);
    Expansion<MultiIndex,F> expansion({{{0,0},2.0_x},{{1,0},3.0_x},{{0,1},-1.0_x},{{2,0},0.5_x}},pr);
    TaylorModel<ValidatedTag,F> taylor(expansion,F(0.25_x,pr),sweeper);
    auto linear=affine_model(taylor);
    ARIADNE_TEST_ASSERT(F::get_rounding_mode()==old_rounding);
    ARIADNE_TEST_EQUALS(linear.value(),F(2,pr));
    ARIADNE_TEST_EQUALS(linear.gradient(1u),F(-1,pr));
    ARIADNE_TEST_EQUALS(linear.error().raw(),F(0.75_x,pr));
    auto models=affine_models(Vector<TaylorModel<ValidatedTag,F>>{taylor});
    ARIADNE_TEST_EQUALS(models[0u].value(),F(2,pr));
}
template<class P,class F> Void test_composition(typename F::PrecisionType pr) {
    using A=AffineModel<P,F>;using C=typename A::CoefficientType;using X=typename A::NumericType;
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    Vector<A> input{A(C(2,pr),{C(1,pr),C(0,pr)}),A(C(3,pr),{C(0,pr),C(1,pr)})};
    auto scalar=compose(x+y,input);auto vector=compose(V{x+y,x-y},input);
    Vector<X> point({0.25_x,0.5_x},pr);
    ARIADNE_TEST_ASSERT(decide(scalar.evaluate(point)==X(5.75_x,pr)));
    ARIADNE_TEST_ASSERT(decide(vector[1u].evaluate(point)==X(-1.25_x,pr)));
    auto nonlinear=compose(sqr(x)+y,input);
    ARIADNE_TEST_EQUALS(nonlinear.value(),C(7,pr));
    ARIADNE_TEST_EQUALS(nonlinear.gradient(1u),C(1,pr));
}
Void test_function_conversion() {
    MultiplePrecision mp(128u);
    ARIADNE_TEST_ASSERT(decide(FloatMPApproximation(3,mp)<FloatDP(4,dp)));
    ARIADNE_TEST_ASSERT(decide(FloatMP(3,mp)<FloatDP(4,dp)));
    BoxDomainType domain{{-1,1},{-1,1}};
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    auto scalar=affine_model(domain,x+y,dp);
    auto vector=affine_models(domain,EffectiveVectorMultivariateFunction{x+y,x-y},dp);
    Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    ARIADNE_TEST_ASSERT(definitely(scalar.evaluate(point)==FloatDPBounds(0.75_x,dp)));
    ARIADNE_TEST_ASSERT(definitely(vector[1u].evaluate(point)==FloatDPBounds(-0.25_x,dp)));
}
Int main() {
    ARIADNE_TEST_CALL(test_product<FloatDP>(dp));
    ARIADNE_TEST_CALL(test_product<FloatMP>(MultiplePrecision(128u)));
    ARIADNE_TEST_CALL(test_approximate_model<FloatDP>(dp));
    ARIADNE_TEST_CALL(test_approximate_model<FloatMP>(MultiplePrecision(128u)));
    ARIADNE_TEST_CALL(test_validated_model<FloatDP>(dp));
    ARIADNE_TEST_CALL(test_validated_model<FloatMP>(MultiplePrecision(128u)));
    ARIADNE_TEST_CALL((test_composition<ApproximateTag,FloatDP>(dp)));
    ARIADNE_TEST_CALL((test_composition<ValidatedTag,FloatDP>(dp)));
    ARIADNE_TEST_CALL(test_function_conversion());
    return ARIADNE_TEST_FAILURES;
}
