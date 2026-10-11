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

inline BenchmarkCase make_affine_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        1+x[0]-2*x[1],
        -1+3*x[0]+x[1]
    });
    return {"affine_2d",f,box2()};
}

inline BenchmarkCase make_polynomial_low_order_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]+x[1]+x[0]*x[1],
        x[0]-x[1]+sqr(x[0])-sqr(x[1])
    });
    return {"polynomial_low_order_2d",f,box2()};
}

inline BenchmarkCase make_polynomial_low_order_4d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        x[0]+x[1]+x[2]*x[3],
        x[0]-x[1]+x[2]+x[0]*x[3],
        x[1]-x[2]+x[3]+x[0]*x[2],
        -x[0]+x[2]-x[3]+x[1]*x[3]
    });
    return {"polynomial_low_order_4d",f,box4()};
}

inline BenchmarkCase make_polynomial_high_order_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[1]+pow(x[0],3u)-pow(x[0],5u)+pow(x[0],3u)*pow(x[1],4u),
        -x[0]+pow(x[1],3u)+pow(x[0],4u)*pow(x[1],3u)-pow(x[1],7u)
    });
    return {"polynomial_high_order_2d",f,box2()};
}

inline BenchmarkCase make_polynomial_high_order_4d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        x[1]+pow(x[0],5u)*x[2]+x[0]*x[1]*x[2]*x[3],
        -x[0]+pow(x[1],4u)*pow(x[3],2u)+pow(x[2],6u),
        x[3]+pow(x[0],2u)*pow(x[2],3u)+pow(x[1],7u),
        -x[2]+pow(x[0],3u)*pow(x[3],3u)+pow(x[2],5u)*x[1]
    });
    return {"polynomial_high_order_4d",f,box4()};
}

inline BenchmarkCase make_transcendental_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        sin(x[0])+exp(x[1])-1,
        cos(x[1])+x[0]*x[1]-1
    });
    return {"transcendental_2d",f,box2()};
}

inline BenchmarkCase make_transcendental_4d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        sin(x[0])+exp(x[1])-1,
        cos(x[1])+sin(x[2])+x[0]*x[3]-1,
        exp(x[2])-1+x[1]*x[3],
        sin(x[3])+cos(x[0])+x[2]-1
    });
    return {"transcendental_4d",f,box4()};
}

inline std::vector<BenchmarkCase> make_cases() {
    std::vector<BenchmarkCase> cases;
    cases.push_back(make_affine_2d());
    cases.push_back(make_polynomial_low_order_2d());
    cases.push_back(make_polynomial_low_order_4d());
    cases.push_back(make_polynomial_high_order_2d());
    cases.push_back(make_polynomial_high_order_4d());
    cases.push_back(make_transcendental_2d());
    cases.push_back(make_transcendental_4d());
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
