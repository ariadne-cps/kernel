/***************************************************************************
 *            test_function_interface.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/function_mixin.tpl.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"
#include <array>
#include <memory>
#include <sstream>
using namespace Ariadne;

template<class F> concept HasVectorIndexing = requires(F& f,const F& cf) {
    cf.get(0u);f.set(0u,cf.create_zero());f[0u];cf[0u];
};
template<class P> Void check_vector_api() {
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    static_assert(!HasVectorIndexing<S> && !HasVectorIndexing<U>);
    static_assert(HasVectorIndexing<V> && HasVectorIndexing<W>);
    static_assert(!Constructible<S,List<S>> && !Constructible<S,Vector<S>>);
    static_assert(!Constructible<U,List<U>> && !Constructible<U,Vector<U>>);
    static_assert(Constructible<V,List<S>> && Constructible<W,List<U>>);
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    V vector{x,y};const V& cv=vector;
    Vector<FloatDPApproximation> point({3,4},dp);
    ARIADNE_TEST_EQUALS(cv[Range(0u,2u)](point),point);
    ARIADNE_TEST_EQUALS(cv[Range(1u,1u)](point).size(),0u);
    ARIADNE_TEST_FAIL(vector.get(2u));
    ARIADNE_TEST_FAIL(cv[Range(0u,3u)]);
    ARIADNE_TEST_FAIL(S::coordinate(2u,2u));
    ARIADNE_TEST_FAIL(vector.set(2u,x));
    ARIADNE_TEST_FAIL(vector.set(0u,S::coordinate(1u,0u)));
    vector[0u]=vector[1u];
    ARIADNE_TEST_EQUALS(vector[0u](point),point[1]);
    ARIADNE_TEST_EQUALS(vector[0u].evaluate(point),point[1]);
    ARIADNE_TEST_EQUALS(vector[0u].domain().dimension(),2u);
    std::ostringstream text;text<<vector[0u];
    ARIADNE_TEST_ASSERT(!text.str().empty());
    auto composed=compose(V::identity(2u),vector);
    ARIADNE_TEST_FAIL(composed.set(0u,x));
    auto t=U::coordinate();W curve{t,t};
    auto composed_curve=compose(curve,t);
    ARIADNE_TEST_FAIL(composed_curve.set(0u,t));
    ARIADNE_TEST_FAIL(curve.get(2u));
    ARIADNE_TEST_FAIL(curve.set(2u,t));
}

// A valid directly implemented vector interface can evaluate and differentiate
// without implementing the optional VectorOfFunctionInterface getter.
class EvaluatedVectorFunction
    : public FunctionMixin<EvaluatedVectorFunction,EffectiveTag,RealVector(RealVector)>
{
    SharedPointer<const EffectiveVectorMultivariateFunction> _function;
  public:
    explicit EvaluatedVectorFunction(SharedPointer<const EffectiveVectorMultivariateFunction> f) : _function(f) { }
    SizeType argument_size() const override { return _function->argument_size(); }
    SizeType result_size() const override { return _function->result_size(); }
    template<class X> Vector<X> operator()(const Vector<X>& x) const { return (*_function)(x); }
    friend EvaluatedVectorFunction derivative(const EvaluatedVectorFunction& f,SizeType j) {
        return EvaluatedVectorFunction(std::make_shared<const EffectiveVectorMultivariateFunction>(derivative(*f._function,j))); }
    friend OutputStream& operator<<(OutputStream& os,const EvaluatedVectorFunction& f) {
        return os<<"EvaluatedVectorFunction("<<*f._function<<")"; }
};
class EvaluatedCurveFunction
    : public FunctionMixin<EvaluatedCurveFunction,EffectiveTag,RealVector(RealScalar)>
{
    SharedPointer<const EffectiveVectorUnivariateFunction> _function;
  public:
    explicit EvaluatedCurveFunction(SharedPointer<const EffectiveVectorUnivariateFunction> f) : _function(f) { }
    SizeOne argument_size() const override { return SizeOne(); }
    SizeType result_size() const override { return _function->result_size(); }
    template<class X> Vector<X> operator()(const X& x) const { return (*_function)(x); }
    friend EvaluatedCurveFunction derivative(const EvaluatedCurveFunction& f,IndexZero j) {
        return EvaluatedCurveFunction(std::make_shared<const EffectiveVectorUnivariateFunction>(derivative(*f._function,j))); }
    friend OutputStream& operator<<(OutputStream& os,const EvaluatedCurveFunction& f) {
        return os<<"EvaluatedCurveFunction("<<*f._function<<")"; }
};
Void test_optional_getter() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    V function(new EvaluatedVectorFunction(std::make_shared<const V>(V{x+y,x-y})));
    Vector<FloatDPApproximation> point({3,1},dp);
    ARIADNE_TEST_EQUALS(function.get(0u)(point),FloatDPApproximation(4,dp));
    ARIADNE_TEST_EQUALS(function.get(1u)(point),FloatDPApproximation(2,dp));
    ARIADNE_TEST_EQUALS(function.get(1u).derivative(1u)(point),FloatDPApproximation(-1,dp));
    ARIADNE_TEST_EQUALS(function.get(0u).result_size(),1u);
    auto check_field=[&](const auto& field) {
        ARIADNE_TEST_EQUALS(field.get(0u)(point),FloatDPApproximation(4,dp));
        ARIADNE_TEST_EQUALS(field.get(1u)(point),FloatDPApproximation(2,dp));
    };
    check_field(ValidatedVectorMultivariateFunction(function));
    check_field(ApproximateVectorMultivariateFunction(function));
    auto t=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction curve(new EvaluatedCurveFunction(
        std::make_shared<const EffectiveVectorUnivariateFunction>(EffectiveVectorUnivariateFunction{t,t+2})));
    FloatDPApproximation parameter(3,dp);
    auto check_curve=[&](const auto& input_curve) {
        ARIADNE_TEST_EQUALS(input_curve.get(0u)(parameter),FloatDPApproximation(3,dp));
        ARIADNE_TEST_EQUALS(input_curve.get(1u)(parameter),FloatDPApproximation(5,dp));
        ARIADNE_TEST_EQUALS(input_curve.get(1u).derivative(IndexZero())(parameter),FloatDPApproximation(1,dp));
    };
    check_curve(curve);
    check_curve(ValidatedVectorUnivariateFunction(curve));
    check_curve(ApproximateVectorUnivariateFunction(curve));
    std::ostringstream text;text<<function.get(0u);
    ARIADNE_TEST_ASSERT(!text.str().empty());
}

class TaylorBackedFactory : public FunctionFactoryInterface<ValidatedTag> {
    TaylorFunctionFactory _factory;
  public:
    TaylorBackedFactory() : _factory(Sweeper<FloatDP>()) { }
    TaylorBackedFactory* clone() const override { return new TaylorBackedFactory(*this); }
    OutputStream& _write(OutputStream& os) const override { return os<<"TaylorBackedFactory"; }
  private:
    ValidatedScalarMultivariateFunction::Interface* _create(const BoxDomainType& domain,const ValidatedScalarMultivariateFunction::Interface& f) const override {
        auto model=_factory.create(domain,f);return new decltype(model)(model); }
    ValidatedVectorMultivariateFunction::Interface* _create(const BoxDomainType& domain,const ValidatedVectorMultivariateFunction::Interface& f) const override {
        auto model=_factory.create(domain,f);return new decltype(model)(model); }
};
Void test_function_factory() {
    TaylorBackedFactory backend;
    std::array<ValidatedFunctionFactory,3> factories{
        ValidatedFunctionFactory(backend),
        ValidatedFunctionFactory(new TaylorBackedFactory()),
        ValidatedFunctionFactory(std::make_shared<const TaylorBackedFactory>())};
    BoxDomainType domain{{-1,1},{-1,1}};
    auto source=ValidatedScalarMultivariateFunction::coordinate(2u,0u)+ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    auto identity=ValidatedVectorMultivariateFunction::identity(2u);
    Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    for(const auto& factory:factories) {
        ARIADNE_TEST_ASSERT(definitely(factory.create(domain,source.reference())(point)==FloatDPBounds(0.75_x,dp)));
        ARIADNE_TEST_EQUALS(factory.create(domain,identity.reference())(point),point);
        ARIADNE_TEST_ASSERT(definitely(factory.create_zero(domain)(point)==FloatDPBounds(0,dp)));
        ARIADNE_TEST_EQUALS(factory.create_identity(domain)(point),point);
        std::ostringstream text;text<<factory;
        ARIADNE_TEST_EQUALS(text.str(),String("TaylorBackedFactory"));
    }
    std::unique_ptr<FunctionFactoryInterface<ValidatedTag>> copy(backend._copy());
    ARIADNE_TEST_EQUALS(copy->create_identity(domain)(point),point);
    std::ostringstream text;text<<*copy;
    ARIADNE_TEST_EQUALS(text.str(),String("TaylorBackedFactory"));
}
Int main() {
    ARIADNE_TEST_CALL(test_optional_getter());
    ARIADNE_TEST_CALL(check_vector_api<ApproximateTag>());
    ARIADNE_TEST_CALL(check_vector_api<ValidatedTag>());
    ARIADNE_TEST_CALL(check_vector_api<EffectiveTag>());
    ARIADNE_TEST_CALL(test_function_factory());
    return ARIADNE_TEST_FAILURES;
}
