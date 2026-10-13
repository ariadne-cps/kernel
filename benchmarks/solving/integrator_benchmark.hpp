#pragma once

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ariadne-kernel.hpp"

namespace IntegratorBenchmark {

using namespace Ariadne;

struct BenchmarkCase {
    std::string name;
    ValidatedVectorMultivariateFunction vector_field;
    ExactBoxType domain;
};

struct NamedIntegrator {
    std::string key;
    std::string name;
    std::unique_ptr<IntegratorInterface> integrator;
};

inline ExactBoxType box2() {
    return ExactBoxType({
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    });
}

inline ExactBoxType box4() {
    return ExactBoxType({
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    });
}

inline BenchmarkCase make_affine_2d_dense_rotating() {
    // A1: isolate geometric wrapping in a linear system.  The dense,
    // non-normal linear part rotates and shears sets without introducing
    // polynomial truncation error.  A stateful QR-preconditioned method
    // should benefit from an oriented affine/Taylor-model representation,
    // whereas an axis-aligned box is expected to wrap.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]-2*x[1],
        2*x[0]-x[1]
    });
    return {"affine_2d_dense_rotating",f,box2()};
}

inline BenchmarkCase make_affine_4d_sparse_multiscale() {
    // A2: isolate time-scale separation in an axis-aligned, sparse affine
    // system with only weak one-way coupling.  This should favour methods
    // that preserve simple scaled coordinates; QR rotation is not expected
    // to help much while the principal directions remain axis-aligned.
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        -16*x[0]+x[1],
        -4*x[1],
        -x[2]+x[3]/4,
        -x[3]/4
    });
    return {"affine_4d_sparse_multiscale",f,box4()};
}

inline BenchmarkCase make_polynomial_low_2d_sparse_weak() {
    // P1: baseline low-degree nonlinear case.  Sparsity and weak coupling
    // keep mixed dependencies limited and the spatial degree of the flow
    // grows gradually.  A graded Taylor representation should be able to
    // exploit this structure without requiring geometric preconditioning.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]+sqr(x[0]),
        -x[1]+x[0]/4
    });
    return {"polynomial_low_2d_sparse_weak",f,box2()};
}

inline BenchmarkCase make_polynomial_low_2d_dense_strong() {
    // P2: keep degree low but make both components strongly coupled.  The
    // purpose is to stress the generation of mixed monomials independently
    // of high polynomial degree.  A full multivariate Taylor model, and in
    // particular a total-degree series method, should be a natural fit.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]+x[1]+x[0]*x[1],
        x[0]-x[1]+sqr(x[0])-sqr(x[1])
    });
    return {"polynomial_low_2d_dense_strong",f,box2()};
}

inline BenchmarkCase make_polynomial_low_2d_sparse_fastgrowth() {
    // P3: low vector-field degree but fast growth of spatial degree in the
    // time-series coefficients, with sparse asymmetric coupling.  This is
    // intended to expose insufficient MaximumSpacialOrder in graded methods;
    // an oriented set representation cannot compensate for missing spatial
    // polynomial order.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[1],
        x[1]-x[0]-sqr(x[0])*x[1]
    });
    return {"polynomial_low_2d_sparse_fastgrowth",f,box2()};
}

inline BenchmarkCase make_polynomial_low_2d_rotating_anisotropic() {
    // G1: moderate polynomial complexity with a strongly non-normal linear
    // part.  The main stress is anisotropy plus rotation/shear rather than
    // rapid degree growth.  This is the case where stateful QR
    // preconditioning should have the clearest opportunity to outperform an
    // axis-aligned Taylor-model representation.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]-8*x[1]+sqr(x[0])/4,
        x[0]-x[1]
    });
    return {"polynomial_low_2d_rotating_anisotropic",f,box2()};
}

inline BenchmarkCase make_polynomial_low_4d_sparse_strong_propagation() {
    // P6: the vector field is sparse, but coupling forms a nonlinear ring so
    // successive time derivatives propagate dependencies through all four
    // variables.  This stresses flow densification rather than initial
    // density.  Sparse/graded representations may win early, while a full
    // multivariate Taylor model becomes relevant as mixed terms accumulate.
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        x[1]+x[0]*x[1],
        x[2]+x[1]*x[2],
        x[3]+x[2]*x[3],
        -x[0]+x[3]*x[0]
    });
    return {"polynomial_low_4d_sparse_strong_propagation",f,box4()};
}

inline BenchmarkCase make_analytic_2d_dense() {
    // T1: dense non-polynomial analytic dynamics.  This stresses elementary
    // function expansion and rigorous remainder bounds rather than finite
    // polynomial degree.  Picard/series-bounder approaches are expected to
    // be competitive when paired with Taylor models carrying interval
    // remainders.
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        sin(x[0])+exp(x[1])-1,
        cos(x[1])+x[0]*x[1]-1
    });
    return {"analytic_2d_dense",f,box2()};
}

inline std::vector<BenchmarkCase> make_cases() {
    std::vector<BenchmarkCase> cases;
    cases.push_back(make_affine_2d_dense_rotating());
    cases.push_back(make_affine_4d_sparse_multiscale());
    cases.push_back(make_polynomial_low_2d_sparse_weak());
    cases.push_back(make_polynomial_low_2d_dense_strong());
    cases.push_back(make_polynomial_low_2d_sparse_fastgrowth());
    cases.push_back(make_polynomial_low_2d_rotating_anisotropic());
    cases.push_back(make_polynomial_low_4d_sparse_strong_propagation());
    cases.push_back(make_analytic_2d_dense());
    return cases;
}

inline std::vector<NamedIntegrator> make_all_integrators(double step_error) {
    constexpr DegreeType order_value=8;

    std::vector<NamedIntegrator> integrators;
    integrators.push_back({
        "taylor-picard",
        "TaylorPicardIntegrator",
        std::make_unique<TaylorPicardIntegrator>(StepMaximumError(step_error))
    });
    integrators.push_back({
        "graded-taylor-picard",
        "GradedTaylorPicardIntegrator",
        std::make_unique<GradedTaylorPicardIntegrator>(
            StepMaximumError(step_error),Order(order_value))
    });
    integrators.push_back({
        "taylor-series",
        "TaylorSeriesIntegrator",
        std::make_unique<TaylorSeriesIntegrator>(
            StepMaximumError(step_error),Order(order_value))
    });
    integrators.push_back({
        "taylor-series-bounder",
        "TaylorSeriesBounderIntegrator",
        std::make_unique<TaylorSeriesBounderIntegrator>(
            StepMaximumError(step_error),Order(order_value))
    });
    integrators.push_back({
        "graded-taylor-series",
        "GradedTaylorSeriesIntegrator",
        std::make_unique<GradedTaylorSeriesIntegrator>(
            StepMaximumError(step_error))
    });
    integrators.push_back({
        "preconditioned-graded-taylor-series",
        "PreconditionedGradedTaylorSeriesIntegrator",
        std::make_unique<PreconditionedGradedTaylorSeriesIntegrator>(
            StepMaximumError(step_error))
    });
    return integrators;
}

inline std::vector<NamedIntegrator> select_integrators(
        const std::string& selection,
        double step_error)
{
    auto all=make_all_integrators(step_error);
    if(selection=="all") {
        return all;
    }

    std::vector<NamedIntegrator> selected;
    for(auto& named_integrator : all) {
        if(named_integrator.key==selection) {
            selected.push_back(std::move(named_integrator));
            return selected;
        }
    }

    std::cerr << "Unknown integrator '" << selection << "'.\n";
    std::cerr << "Available integrators:\n";
    std::cerr << "  all\n";
    for(const auto& named_integrator : all) {
        std::cerr << "  " << named_integrator.key << "\n";
    }
    std::exit(1);
}

} // namespace IntegratorBenchmark
