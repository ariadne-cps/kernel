/***************************************************************************
 *            test_paver.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "geometry/paver.hpp"
#include "geometry/function_set.hpp"
#include "geometry/grid_paving.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

Void check_paver(const Paver& paver) {
    ARIADNE_TEST_PRINT(paver);
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction identity({x});
    ExactBoxType domain{{0.125_x,0.375_x}};
    ValidatedConstrainedImageSet set(domain,identity);
    Grid grid(1u);
    GridTreePaving approximation(grid);
    paver.adjoin_outer_approximation(approximation,set,2u);
    ARIADNE_TEST_ASSERT(approximation.superset(GridCell(grid,0u,BinaryWord("00"))));
    ARIADNE_TEST_ASSERT(approximation.superset(GridCell(grid,0u,BinaryWord("01"))));
    ARIADNE_TEST_ASSERT(definitely(subset(cast_exact_box(approximation.bounding_box()),ExactBoxType{{-0.5_x,0.75_x}})));
    GridTreePaving full(grid);
    full.adjoin(GridCell(grid,0u,BinaryWord()));
    paver.adjoin_outer_approximation(full,set,2u);
    ARIADNE_TEST_ASSERT(full.superset(GridCell(grid,0u,BinaryWord())));
    ValidatedConstrainedImageSet constrained(domain,identity,List<ValidatedConstraint>{0<=x<=1});
    GridTreePaving constrained_result(grid);
    paver.adjoin_outer_approximation(constrained_result,constrained,2u);
    ARIADNE_TEST_ASSERT(constrained_result.superset(GridCell(grid,0u,BinaryWord("00"))));
    ValidatedConstrainedImageSet below(domain,identity,List<ValidatedConstraint>{x<=-1});
    ValidatedConstrainedImageSet above(domain,identity,List<ValidatedConstraint>{x>=1});
    GridTreePaving excluded(grid);
    paver.adjoin_outer_approximation(excluded,below,2u);
    paver.adjoin_outer_approximation(excluded,above,2u);
    ARIADNE_TEST_ASSERT(excluded.is_empty());
    ValidatedConstrainedImageSet empty(ExactBoxType(1u),identity);
    paver.adjoin_outer_approximation(excluded,empty,2u);
    ARIADNE_TEST_ASSERT(excluded.is_empty());
    GridTreePaving wrong_dimension(Grid(2u));
    ARIADNE_TEST_FAIL(paver.adjoin_outer_approximation(wrong_dimension,set,2u));
    std::ostringstream repr;
    repr<<paver;
    ARIADNE_TEST_ASSERT(repr.str().find("Paver()")!=String::npos);
}

Void test_paver_strategies() {
    ARIADNE_TEST_CALL(check_paver(Paver(AffinePaver())));
    ARIADNE_TEST_CALL(check_paver(Paver(SubdivisionPaver())));
    ARIADNE_TEST_CALL(check_paver(Paver(ReducePaver())));
    ARIADNE_TEST_CALL(check_paver(Paver(ConstraintPaver())));
    ARIADNE_TEST_CALL(check_paver(Paver(OptimalConstraintPaver())));
}

Void test_nonlinear_pavers() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction quadratic({x*x});
    ValidatedConstrainedImageSet set(ExactBoxType{{-1.0_x,1.0_x}},quadratic);
    Grid grid(1u);
    const Paver strategies[]={Paver(AffinePaver()),Paver(SubdivisionPaver()),Paver(ReducePaver()),Paver(ConstraintPaver()),Paver(OptimalConstraintPaver())};
    for(const auto& paver:strategies) {
        ARIADNE_TEST_PRINT(paver);
        GridTreePaving approximation(grid);
        paver.adjoin_outer_approximation(approximation,set,1u);
        const ExactBoxType bounds=cast_exact_box(approximation.bounding_box());
        ARIADNE_TEST_ASSERT(bounds[0].lower_bound()<=FloatDP(0,dp));
        ARIADNE_TEST_ASSERT(bounds[0].upper_bound()>=FloatDP(1,dp));
        ARIADNE_TEST_ASSERT(definitely(subset(bounds,ExactBoxType{{-2.0_x,2.0_x}})));
    }
}

Void test_intrinsic_uncertainty() {
    const Paver strategies[]={Paver(AffinePaver()),Paver(SubdivisionPaver()),Paver(ReducePaver()),Paver(ConstraintPaver()),Paver(OptimalConstraintPaver())};
    for(SizeType arity=0u; arity!=2u; ++arity) {
        auto constant=ValidatedScalarMultivariateFunction::constant(arity,ValidatedNumber(FloatDPBounds(0.25_x,0.75_x,dp)));
        ValidatedVectorMultivariateFunction mapping({constant});
        ValidatedConstrainedImageSet set(ExactBoxType(arity,ExactIntervalType(-1,1)),mapping);
        for(const auto& paver:strategies) {
            ARIADNE_TEST_PRINT(arity);
            ARIADNE_TEST_PRINT(paver);
            GridTreePaving approximation(Grid(1u));
            paver.adjoin_outer_approximation(approximation,set,3u);
            const auto bounds=cast_exact_box(approximation.bounding_box());
            ARIADNE_TEST_ASSERT(bounds[0].lower_bound()<=FloatDP(0.25_x,dp));
            ARIADNE_TEST_ASSERT(bounds[0].upper_bound()>=FloatDP(0.75_x,dp));
            ARIADNE_TEST_ASSERT(definitely(subset(bounds,ExactBoxType{{0.0_x,1.0_x}})));
        }
    }
    auto identity=ValidatedVectorMultivariateFunction::identity(1u);
    ValidatedConstrainedImageSet set(ExactBoxType{{0.125_x,0.375_x}},identity);
    GridTreePaving paving(Grid(1u));
    SubdivisionPaver subdivision;
    ARIADNE_TEST_FAIL(subdivision.adjoin_outer_approximation_recursion(paving,set,2u,Vector<FloatDP>(0u,dp)));
    ARIADNE_TEST_FAIL(subdivision.adjoin_outer_approximation_recursion(paving,set,2u,Vector<FloatDP>(1u,FloatDP(0,dp))));
}

Int main() {
    ARIADNE_TEST_CALL(test_intrinsic_uncertainty());
    ARIADNE_TEST_CALL(test_paver_strategies());
    ARIADNE_TEST_CALL(test_nonlinear_pavers());
    return ARIADNE_TEST_FAILURES;
}
