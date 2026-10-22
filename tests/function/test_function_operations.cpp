/***************************************************************************
 *            test_function_operations.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/symbolic_function.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;
using X=FloatDPApproximation;

template<class F> Void check_display(const F& f) {
    std::ostringstream text; text<<f<<" "<<representation(f);
    ARIADNE_TEST_ASSERT(!text.str().empty());
}

template<class P> Void test_constructors() {
    using C=FunctionConstructors<P>;using Y=Number<P>;
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    const X scalar(3,dp);Vector<X> point({3,2},dp);
    EuclideanDomain domain(2u);RealDomain real;
    S null_scalar; V null_vector; U null_univariate; W null_curve;
    ARIADNE_TEST_ASSERT(null_scalar.raw_pointer()==nullptr);
    ARIADNE_TEST_ASSERT(null_vector.raw_pointer()==nullptr);
    ARIADNE_TEST_ASSERT(null_univariate.raw_pointer()==nullptr);
    ARIADNE_TEST_ASSERT(null_curve.raw_pointer()==nullptr);
    ARIADNE_TEST_EQUALS(C::zero()(scalar),X(0,dp));
    ARIADNE_TEST_EQUALS(C::zero(real)(scalar),X(0,dp));
    ARIADNE_TEST_EQUALS(C::constant(Y(2))(scalar),X(2,dp));
    ARIADNE_TEST_EQUALS(C::constant(real,Y(2))(scalar),X(2,dp));
    ARIADNE_TEST_EQUALS(C::coordinate()(scalar),scalar);
    ARIADNE_TEST_EQUALS(C::coordinate(real)(scalar),scalar);
    ARIADNE_TEST_EQUALS(C::coordinate(real,IndexZero())(scalar),scalar);
    ARIADNE_TEST_EQUALS(C::identity()(scalar),scalar);
    ARIADNE_TEST_EQUALS(C::identity(real)(scalar),scalar);
    ARIADNE_TEST_EQUALS(C::zero(2u)(point),X(0,dp));
    ARIADNE_TEST_EQUALS(C::zero(domain)(point),X(0,dp));
    ARIADNE_TEST_EQUALS(C::constant(2u,Y(2))(point),X(2,dp));
    ARIADNE_TEST_EQUALS(C::constant(domain,Y(2))(point),X(2,dp));
    ARIADNE_TEST_EQUALS(C::coordinate(2u,1u)(point),point[1]);
    ARIADNE_TEST_EQUALS(C::coordinate(domain,1u)(point),point[1]);
    ARIADNE_TEST_EQUALS(C::coordinates(2u)[1](point),point[1]);
    ARIADNE_TEST_EQUALS(C::coordinates(domain)[0](point),point[0]);
    ARIADNE_TEST_ASSERT(C::coordinates(0u).empty());
    ARIADNE_TEST_ASSERT(C::coordinates(EuclideanDomain(0u)).empty());
    ARIADNE_TEST_EQUALS(C::identity(2u)(point),point);
    ARIADNE_TEST_EQUALS(C::identity(domain)(point),point);
    Vector<X> zeroes(2u,dp);Vector<Y> constants{Y(2),Y(4)};
    Vector<X> expected({2,4},dp);
    ARIADNE_TEST_EQUALS(C::zeros(2u,2u)(point),zeroes);
    ARIADNE_TEST_EQUALS(C::zeros(2u,domain)(point),zeroes);
    ARIADNE_TEST_EQUALS(C::zeros(2u)(scalar),zeroes);
    ARIADNE_TEST_EQUALS(C::zeros(2u,real)(scalar),zeroes);
    ARIADNE_TEST_EQUALS(C::constant(2u,constants)(point),expected);
    ARIADNE_TEST_EQUALS(C::constant(domain,constants)(point),expected);
    ARIADNE_TEST_EQUALS(C::constant(constants)(scalar),expected);
    ARIADNE_TEST_EQUALS(C::constant(real,constants)(scalar),expected);
    ARIADNE_TEST_EQUALS(C::zeros(0u,2u)(point).size(),0u);
    ARIADNE_TEST_EQUALS(C::zeros(0u,domain)(point).size(),0u);
    ARIADNE_TEST_EQUALS(C::zeros(0u)(scalar).size(),0u);
    ARIADNE_TEST_EQUALS(C::zeros(0u,real)(scalar).size(),0u);
    Vector<X> empty_point(0u,dp);
    ARIADNE_TEST_EQUALS(C::identity(0u)(empty_point).size(),0u);
    ARIADNE_TEST_EQUALS(C::identity(EuclideanDomain(0u))(empty_point).size(),0u);
    ARIADNE_TEST_EQUALS(C::constant(2u,Vector<Y>{})(point).size(),0u);
    ARIADNE_TEST_EQUALS(C::constant(domain,Vector<Y>{})(point).size(),0u);
    ARIADNE_TEST_EQUALS(C::constant(real,Vector<Y>{})(scalar).size(),0u);
    S s(domain);U u(real);V v(2u,domain);W w(2u,real);
    ARIADNE_TEST_EQUALS(s(point),X(0,dp));
    ARIADNE_TEST_EQUALS(u(scalar),X(0,dp));
    ARIADNE_TEST_EQUALS(v(point),zeroes);
    ARIADNE_TEST_EQUALS(w(scalar),zeroes);
    S s2(SizeOne(),2u);U u2{SizeOne(),SizeOne()};V v2(2u,2u);W w2(2u,SizeOne());
    ARIADNE_TEST_EQUALS(s2(point),X(0,dp));
    ARIADNE_TEST_EQUALS(u2(scalar),X(0,dp));
    ARIADNE_TEST_EQUALS(v2(point),zeroes);
    ARIADNE_TEST_EQUALS(w2(scalar),zeroes);
    V domain_only(domain);W real_only(real);
    ARIADNE_TEST_EQUALS(domain_only.result_size(),0u);
    ARIADNE_TEST_EQUALS(real_only.result_size(),0u);
    auto x=C::coordinate(domain,0u);auto t=C::coordinate(real);
    ARIADNE_TEST_EQUALS(x.create_zero()(point),X(0,dp));
    ARIADNE_TEST_EQUALS(x.create_constant(Y(2))(point),X(2,dp));
    ARIADNE_TEST_EQUALS(x.create_coordinate(1u)(point),point[1]);
    ARIADNE_TEST_EQUALS(x.create_constant(constants)(point),expected);
    ARIADNE_TEST_EQUALS(t.create_zero()(scalar),X(0,dp));
    ARIADNE_TEST_EQUALS(t.create_constant(Y(2))(scalar),X(2,dp));
    ARIADNE_TEST_EQUALS(t.create_coordinate(IndexZero())(scalar),scalar);
    ARIADNE_TEST_EQUALS(t.create_constant(constants)(scalar),expected);
    V repeated(2u,x);W repeated_curve(2u,t);
    ARIADNE_TEST_EQUALS(repeated(point),Vector<X>({3,3},dp));
    ARIADNE_TEST_EQUALS(repeated_curve(scalar),Vector<X>({3,3},dp));
    V from_list(List<S>{x,x});W curve_list(List<U>{t,t});
    ARIADNE_TEST_EQUALS(from_list(point),repeated(point));
    ARIADNE_TEST_EQUALS(curve_list(scalar),repeated_curve(scalar));
    V from_vector(Vector<S>{x,x});W curve_vector(Vector<U>{t,t});
    ARIADNE_TEST_EQUALS(from_vector(point),repeated(point));
    ARIADNE_TEST_EQUALS(curve_vector(scalar),repeated_curve(scalar));
    s=Y(4);u=Y(4);v=constants;w=constants;
    ARIADNE_TEST_EQUALS(s(point),X(4,dp));
    ARIADNE_TEST_EQUALS(u(scalar),X(4,dp));
    ARIADNE_TEST_EQUALS(v(point),expected);
    ARIADNE_TEST_EQUALS(w(scalar),expected);
    ARIADNE_TEST_EQUALS(x.domain().dimension(),2u);
    ARIADNE_TEST_EQUALS(x.codomain().dimension(),1u);
    ARIADNE_TEST_EQUALS(v.codomain().dimension(),2u);
    check_display(x);check_display(v);check_display(t);check_display(w);
}

template<class P> Void test_composition() {
    using C=FunctionConstructors<P>;using Y=Number<P>;
    using V=VectorMultivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    auto x=C::coordinate(2u,0u),y=C::coordinate(2u,1u);auto t=C::coordinate();
    auto f=Y(2)*x-y+Y(1);auto h=sqr(t)+Y(1);
    V g{x+y,x-y};W q{t,t+Y(2)};
    X scalar(3,dp);Vector<X> point({3,2},dp);
    auto hh=compose(h,h),fq=compose(f,q);auto qh=compose(q,h),gq=compose(g,q);
    auto hf=compose(h,f),fg=compose(f,g);auto qf=compose(q,f),gg=compose(g,g);
    ARIADNE_TEST_EQUALS(hh(scalar),X(101,dp));
    ARIADNE_TEST_EQUALS(fq(scalar),X(2,dp));
    ARIADNE_TEST_EQUALS(qh(scalar),Vector<X>({10,12},dp));
    ARIADNE_TEST_EQUALS(gq(scalar),Vector<X>({8,-2},dp));
    ARIADNE_TEST_EQUALS(hf(point),X(26,dp));
    ARIADNE_TEST_EQUALS(fg(point),X(10,dp));
    ARIADNE_TEST_EQUALS(qf(point),Vector<X>({5,7},dp));
    ARIADNE_TEST_EQUALS(gg(point),Vector<X>({6,4},dp));
    ARIADNE_TEST_EQUALS(derivative(hh)(scalar),X(120,dp));
    ARIADNE_TEST_EQUALS(derivative(fq)(scalar),X(1,dp));
    ARIADNE_TEST_EQUALS(derivative(qh)(scalar),Vector<X>({6,6},dp));
    ARIADNE_TEST_EQUALS(derivative(gq)(scalar),Vector<X>({2,0},dp));
    ARIADNE_TEST_EQUALS(derivative(hf,0u)(point),X(20,dp));
    ARIADNE_TEST_EQUALS(derivative(fg,0u)(point),X(1,dp));
    ARIADNE_TEST_EQUALS(derivative(qf,0u)(point),Vector<X>({2,2},dp));
    ARIADNE_TEST_EQUALS(derivative(gg,0u)(point),Vector<X>({2,0},dp));
    ARIADNE_TEST_EQUALS(h.slope(scalar),X(6,dp));
    ARIADNE_TEST_EQUALS(q.tangent(scalar),Vector<X>({1,1},dp));
    ARIADNE_TEST_EQUALS(f.gradient(point)[1],X(-1,dp));
    ARIADNE_TEST_EQUALS(g.jacobian(point)[1][1],X(-1,dp));
    ARIADNE_TEST_EQUALS(qh.get(1u)(scalar),X(12,dp));
    ARIADNE_TEST_EQUALS(gg.get(1u)(point),X(4,dp));
    check_display(hh);check_display(fq);check_display(qh);check_display(gq);
    check_display(hf);check_display(fg);check_display(qf);check_display(gg);
    auto both=join(g,g),head=join(x,g),tail=join(g,y),pair=join(x,y);
    ARIADNE_TEST_EQUALS(both(point),Vector<X>({5,1,5,1},dp));
    ARIADNE_TEST_EQUALS(head(point),Vector<X>({3,5,1},dp));
    ARIADNE_TEST_EQUALS(tail(point),Vector<X>({5,1,2},dp));
    ARIADNE_TEST_EQUALS(pair(point),point);
    ARIADNE_TEST_EQUALS(both.get(3u)(point),X(1,dp));
    ARIADNE_TEST_EQUALS(derivative(both,1u)(point),Vector<X>({1,-1,1,-1},dp));
    check_display(both);
    auto scalar_expression=pow(x,3)+Y(2)-Y(4)*x/Y(2);
    ARIADNE_TEST_EQUALS(scalar_expression(point),X(23,dp));
    ARIADNE_TEST_EQUALS((Y(8)/x)(point),X(8,dp)/X(3,dp));
    ARIADNE_TEST_EQUALS((x/Y(2))(point),X(1.5_x,dp));
    ARIADNE_TEST_EQUALS((Y(2)-x)(point),X(-1,dp));
    ARIADNE_TEST_EQUALS((-x)(point),X(-3,dp));
    ARIADNE_TEST_EQUALS((+x)(point),X(3,dp));
    ARIADNE_TEST_EQUALS((Y(2)+x)(point),X(5,dp));
}
Void test_effective_vector_algebra() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    auto scalar=S::constant(2u,2);V vector{x+y,x-y};
    Vector<X> point({3,1},dp),expected({4,2},dp),constants({1,2},dp);
    Vector<EffectiveNumber> c{1,2};
    ARIADNE_TEST_EQUALS((+vector)(point),expected);
    ARIADNE_TEST_EQUALS((-vector)(point),-expected);
    ARIADNE_TEST_EQUALS((vector+vector)(point),Vector<X>({8,4},dp));
    ARIADNE_TEST_EQUALS((vector-vector)(point),Vector<X>({0,0},dp));
    ARIADNE_TEST_EQUALS((scalar*vector)(point),Vector<X>({8,4},dp));
    ARIADNE_TEST_EQUALS((vector*scalar)(point),Vector<X>({8,4},dp));
    ARIADNE_TEST_EQUALS((vector/scalar)(point),Vector<X>({2,1},dp));
    ARIADNE_TEST_EQUALS((c+vector)(point),expected+constants);
    ARIADNE_TEST_EQUALS((vector+c)(point),expected+constants);
    ARIADNE_TEST_EQUALS((c-vector)(point),constants-expected);
    ARIADNE_TEST_EQUALS((vector-c)(point),expected-constants);
    ARIADNE_TEST_EQUALS((EffectiveNumber(2)*vector)(point),Vector<X>({8,4},dp));
    ARIADNE_TEST_EQUALS((vector*EffectiveNumber(2))(point),Vector<X>({8,4},dp));
    ARIADNE_TEST_EQUALS((vector/EffectiveNumber(2))(point),Vector<X>({2,1},dp));
    auto t=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction curve{t,t+2};X parameter(3,dp);
    ARIADNE_TEST_EQUALS((curve+curve)(parameter),Vector<X>({6,10},dp));
    ARIADNE_TEST_EQUALS((curve*t)(parameter),Vector<X>({9,15},dp));
    auto embedded_scalar=embed(1u,x+y,1u);auto embedded_vector=embed(1u,vector,1u);
    Vector<X> extended({99,3,1,88},dp);
    ARIADNE_TEST_EQUALS(embedded_scalar(extended),X(4,dp));
    ARIADNE_TEST_EQUALS(embedded_vector(extended),expected);
    ARIADNE_TEST_EQUALS(embedded_scalar.derivative(0u)(extended),X(0,dp));
    ARIADNE_TEST_EQUALS(embedded_scalar.derivative(1u)(extended),X(1,dp));
    ARIADNE_TEST_EQUALS(embedded_scalar.derivative(3u)(extended),X(0,dp));
    ARIADNE_TEST_EQUALS(embedded_vector.derivative(0u)(extended),Vector<X>({0,0},dp));
    ARIADNE_TEST_EQUALS(embedded_vector.derivative(1u)(extended),Vector<X>({1,1},dp));
    ARIADNE_TEST_EQUALS(embedded_vector.get(1u)(extended),X(2,dp));
    auto f=sqr(x)+y;
    ARIADNE_TEST_EQUALS(derivatives(f)(point),Vector<X>({6,1},dp));
    ARIADNE_TEST_EQUALS(lie_derivative(f,vector)(point),X(26,dp));
    ARIADNE_TEST_EQUALS(lie_derivative(vector,vector)(point),Vector<X>({6,2},dp));
    ARIADNE_TEST_EQUALS(evaluate(make_formula(f),point),X(10,dp));
    ARIADNE_TEST_EQUALS(evaluate(make_formula(vector),point),expected);
    auto fx=Formula<EffectiveNumber>::coordinate(0u),fy=Formula<EffectiveNumber>::coordinate(1u);
    V formula_vector(EuclideanDomain(2u),Vector<Formula<EffectiveNumber>>{fx+fy,fx-fy});
    ARIADNE_TEST_EQUALS(evaluate(make_formula(formula_vector),point),expected);
    ARIADNE_TEST_EQUALS(evaluate(make_formula(compose(vector,vector)),point),Vector<X>({6,2},dp));
    auto zero_vector=V::zeros(0u,0u);
    ARIADNE_TEST_EQUALS(make_formula(zero_vector).size(),0u);
    check_display(embedded_scalar);check_display(embedded_vector);
}
Void test_validated_curve_algebra() {
    auto t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedVectorUnivariateFunction curve{t,t+ValidatedNumber(2)};X parameter(3,dp);
    ARIADNE_TEST_EQUALS((curve+curve)(parameter),Vector<X>({6,10},dp));
    ARIADNE_TEST_EQUALS((curve*t)(parameter),Vector<X>({9,15},dp));
}
Int main() {
    ARIADNE_TEST_CALL(test_constructors<ApproximateTag>());
    ARIADNE_TEST_CALL(test_constructors<ValidatedTag>());
    ARIADNE_TEST_CALL(test_constructors<EffectiveTag>());
    ARIADNE_TEST_CALL(test_composition<ApproximateTag>());
    ARIADNE_TEST_CALL(test_composition<ValidatedTag>());
    ARIADNE_TEST_CALL(test_composition<EffectiveTag>());
    ARIADNE_TEST_CALL(test_effective_vector_algebra());
    ARIADNE_TEST_CALL(test_validated_curve_algebra());
    return ARIADNE_TEST_FAILURES;
}
