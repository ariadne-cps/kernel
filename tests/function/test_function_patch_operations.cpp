/***************************************************************************
 *            test_function_patch_operations.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/function_patch.hpp"
#include "macro/testing.hpp"
#include <array>
#include <sstream>
using namespace Ariadne;

Void test_patch_constant_dispatch() {
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    auto model=factory.create_constant(domain,ValidatedNumber(2));
    ValidatedScalarMultivariateFunction f=model;
    Vector<FloatDPBounds> point({0,0},dp);
    ARIADNE_TEST_ASSERT(definitely((f+ValidatedNumber(3))(point)==FloatDPBounds(5,dp)));
    using A=AlgebraOperations<ValidatedScalarMultivariateFunction,ValidatedNumber>;
    ARIADNE_TEST_ASSERT(definitely(A::apply(BinaryElementaryOperator(Sub()),ValidatedNumber(3),f)(point)==FloatDPBounds(1,dp)));
    ARIADNE_TEST_ASSERT(definitely((f-ValidatedNumber(3))(point)==FloatDPBounds(-1,dp)));
    ARIADNE_TEST_ASSERT(definitely((f*ValidatedNumber(3))(point)==FloatDPBounds(6,dp)));
    ARIADNE_TEST_ASSERT(definitely((f/ValidatedNumber(2))(point)==FloatDPBounds(1,dp)));
}
Void test_patch_vector_product() {
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    ValidatedScalarMultivariateFunction scalar=factory.create_constant(domain,ValidatedNumber(2));
    ValidatedVectorMultivariateFunction vector=factory.create_constants(domain,Vector<ValidatedNumber>{4,8});
    Vector<FloatDPBounds> point({0,0},dp);
    ARIADNE_TEST_EQUALS((scalar*vector)(point),Vector<FloatDPBounds>({8,16},dp));
}
Void test_mixed_patch_operations() {
    using S=ValidatedScalarMultivariateFunction;
    using V=ValidatedVectorMultivariateFunction;
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    std::array<S,2> scalar_a{S::constant(2u,ValidatedNumber(2)),factory.create_constant(domain,ValidatedNumber(2))};
    std::array<S,2> scalar_b{S::constant(2u,ValidatedNumber(4)),factory.create_constant(domain,ValidatedNumber(4))};
    Vector<ValidatedNumber> va{2,6},vb{4,8};
    std::array<V,2> vector_a{V::constant(2u,va),factory.create_constants(domain,va)};
    std::array<V,2> vector_b{V::constant(2u,vb),factory.create_constants(domain,vb)};
    Vector<FloatDPBounds> point({0,0},dp);
    for(SizeType i=0;i!=2u;++i) {
        const auto& a=scalar_a[i];
        ARIADNE_TEST_ASSERT(definitely((-a)(point)==FloatDPBounds(-2,dp)));
        ARIADNE_TEST_ASSERT(definitely(pow(a,2)(point)==FloatDPBounds(4,dp)));
        ARIADNE_TEST_ASSERT(definitely((ValidatedNumber(3)+a)(point)==FloatDPBounds(5,dp)));
        ARIADNE_TEST_ASSERT(definitely((ValidatedNumber(3)-a)(point)==FloatDPBounds(1,dp)));
        for(SizeType j=0;j!=2u;++j) {
            const auto& b=scalar_b[j];const auto& av=vector_a[i];const auto& bv=vector_b[j];
            ARIADNE_TEST_ASSERT(definitely((a+b)(point)==FloatDPBounds(6,dp)));
            ARIADNE_TEST_ASSERT(definitely((a-b)(point)==FloatDPBounds(-2,dp)));
            ARIADNE_TEST_ASSERT(definitely((a*b)(point)==FloatDPBounds(8,dp)));
            ARIADNE_TEST_ASSERT(definitely((a/b)(point)==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_EQUALS((av+bv)(point),Vector<FloatDPBounds>({6,14},dp));
            ARIADNE_TEST_EQUALS((av-bv)(point),Vector<FloatDPBounds>({-2,-2},dp));
            ARIADNE_TEST_EQUALS((a*bv)(point),Vector<FloatDPBounds>({8,16},dp));
            ARIADNE_TEST_EQUALS((bv*a)(point),Vector<FloatDPBounds>({8,16},dp));
            ARIADNE_TEST_EQUALS((bv/a)(point),Vector<FloatDPBounds>({2,4},dp));
            ARIADNE_TEST_EQUALS(join(a,b)(point),Vector<FloatDPBounds>({2,4},dp));
            ARIADNE_TEST_EQUALS(join(a,bv)(point),Vector<FloatDPBounds>({2,4,8},dp));
            ARIADNE_TEST_EQUALS(join(av,b)(point),Vector<FloatDPBounds>({2,6,4},dp));
            ARIADNE_TEST_EQUALS(join(av,bv)(point),Vector<FloatDPBounds>({2,6,4,8},dp));
        }
    }
    ValidatedScalarMultivariateFunctionPatch patch=factory.create_constant(domain,ValidatedNumber(2));
    ValidatedScalarMultivariateFunctionModelDP model=factory.create_constant(domain,ValidatedNumber(2));
    S scalar_model=model;
    V one_component{scalar_model};
    ARIADNE_TEST_EQUALS(one_component(point),Vector<FloatDPBounds>{FloatDPBounds(2,dp)});
    std::ostringstream text;text<<patch<<model<<factory;
    ARIADNE_TEST_ASSERT(!text.str().empty());
}
Int main() {
    ARIADNE_TEST_CALL(test_patch_constant_dispatch());
    ARIADNE_TEST_CALL(test_patch_vector_product());
    ARIADNE_TEST_CALL(test_mixed_patch_operations());
    return ARIADNE_TEST_FAILURES;
}
