#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ariadne-kernel.hpp"

using namespace Ariadne;

namespace {

using Clock = std::chrono::steady_clock;
using Duration = std::chrono::nanoseconds;

constexpr auto minimum_total_measurement_time = std::chrono::milliseconds(1);
constexpr Nat maximum_repetitions = 1024u;

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

struct Measurement {
    Duration elapsed;
    FlowStepModelType flow;
};

ExactBoxType box2() {
    return ExactBoxType({
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    });
}

ExactBoxType box4() {
    return ExactBoxType({
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    });
}

BenchmarkCase make_affine_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        1+x[0]-2*x[1],
        -1+3*x[0]+x[1]
    });
    return {"affine_2d",f,box2()};
}

BenchmarkCase make_polynomial_low_order_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[0]+x[1]+x[0]*x[1],
        x[0]-x[1]+sqr(x[0])-sqr(x[1])
    });
    return {"polynomial_low_order_2d",f,box2()};
}

BenchmarkCase make_polynomial_high_order_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[1]+pow(x[0],3u)-pow(x[0],5u)+pow(x[0],3u)*pow(x[1],4u),
        -x[0]+pow(x[1],3u)+pow(x[0],4u)*pow(x[1],3u)-pow(x[1],7u)
    });
    return {"polynomial_high_order_2d",f,box2()};
}

BenchmarkCase make_polynomial_high_order_4d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        x[1]+pow(x[0],5u)*x[2]+x[0]*x[1]*x[2]*x[3],
        -x[0]+pow(x[1],4u)*pow(x[3],2u)+pow(x[2],6u),
        x[3]+pow(x[0],2u)*pow(x[2],3u)+pow(x[1],7u),
        -x[2]+pow(x[0],3u)*pow(x[3],3u)+pow(x[2],5u)*x[1]
    });
    return {"polynomial_high_order_4d",f,box4()};
}

BenchmarkCase make_transcendental_2d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        sin(x[0])+exp(x[1])-1,
        cos(x[1])+x[0]*x[1]-1
    });
    return {"transcendental_2d",f,box2()};
}

BenchmarkCase make_transcendental_4d() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(4);
    ValidatedVectorMultivariateFunction f({
        sin(x[0])+exp(x[1])-1,
        cos(x[1])+sin(x[2])+x[0]*x[3]-1,
        exp(x[2])-1+x[1]*x[3],
        sin(x[3])+cos(x[0])+x[2]-1
    });
    return {"transcendental_4d",f,box4()};
}

std::vector<BenchmarkCase> make_cases() {
    std::vector<BenchmarkCase> cases;
    cases.push_back(make_affine_2d());
    cases.push_back(make_polynomial_low_order_2d());
    cases.push_back(make_polynomial_high_order_2d());
    cases.push_back(make_polynomial_high_order_4d());
    cases.push_back(make_transcendental_2d());
    cases.push_back(make_transcendental_4d());
    return cases;
}

std::vector<NamedIntegrator> make_all_integrators() {
    constexpr double error=1e-8;
    constexpr DegreeType order_value=8;

    std::vector<NamedIntegrator> integrators;
    integrators.push_back({
        "taylor-picard",
        "TaylorPicardIntegrator",
        std::make_unique<TaylorPicardIntegrator>(StepMaximumError(error))
    });
    integrators.push_back({
        "graded-taylor-picard",
        "GradedTaylorPicardIntegrator",
        std::make_unique<GradedTaylorPicardIntegrator>(
            StepMaximumError(error),Order(order_value))
    });
    integrators.push_back({
        "taylor-series",
        "TaylorSeriesIntegrator",
        std::make_unique<TaylorSeriesIntegrator>(
            StepMaximumError(error),Order(order_value))
    });
    integrators.push_back({
        "taylor-series-bounder",
        "TaylorSeriesBounderIntegrator",
        std::make_unique<TaylorSeriesBounderIntegrator>(
            StepMaximumError(error),Order(order_value))
    });
    integrators.push_back({
        "graded-taylor-series",
        "GradedTaylorSeriesIntegrator",
        std::make_unique<GradedTaylorSeriesIntegrator>(
            StepMaximumError(error))
    });
    integrators.push_back({
        "preconditioned-graded-taylor-series",
        "PreconditionedGradedTaylorSeriesIntegrator",
        std::make_unique<PreconditionedGradedTaylorSeriesIntegrator>(
            StepMaximumError(error))
    });
    return integrators;
}

std::vector<NamedIntegrator> select_integrators(const std::string& selection) {
    auto all=make_all_integrators();
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

Bool accepts_step(
        const IntegratorInterface& integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step)
{
    try {
        integrator.flow_step(
            benchmark_case.vector_field,
            benchmark_case.domain,
            step);
        return true;
    } catch(const FlowTimeStepException&) {
        return false;
    } catch(const IncompleteFlowException&) {
        return false;
    } catch(const FlowBoundsException&) {
        return false;
    }
}

StepSizeType determine_step(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases)
{
    StepSizeType step=0.03125_x;

    for(Nat attempt=0u; attempt!=16u; ++attempt) {
        Bool accepted=true;

        for(const auto& named_integrator : integrators) {
            for(const auto& benchmark_case : cases) {
                if(!accepts_step(
                        *named_integrator.integrator,
                        benchmark_case,
                        step)) {
                    accepted=false;
                }
            }
        }

        if(accepted) {
            return step;
        }
        step=hlf(step);
    }

    throw std::runtime_error(
        "Unable to determine a common flow step for all benchmark cases.");
}

Measurement measure(
        const IntegratorInterface& integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step,
        Nat repetitions)
{
    FlowStepModelType flow=integrator.flow_step(
        benchmark_case.vector_field,
        benchmark_case.domain,
        step);

    const auto start=Clock::now();
    for(Nat i=0u; i!=repetitions; ++i) {
        flow=integrator.flow_step(
            benchmark_case.vector_field,
            benchmark_case.domain,
            step);
    }
    const auto stop=Clock::now();

    return {
        std::chrono::duration_cast<Duration>(stop-start),
        std::move(flow)
    };
}

Nat determine_repetitions(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step)
{
    Nat repetitions=1u;
    while(true) {
        Bool all_sufficient=true;

        for(const auto& named_integrator : integrators) {
            for(const auto& benchmark_case : cases) {
                const auto measurement=measure(
                    *named_integrator.integrator,
                    benchmark_case,
                    step,
                    repetitions);
                if(measurement.elapsed<minimum_total_measurement_time) {
                    all_sufficient=false;
                }
            }
        }

        if(all_sufficient || repetitions>=maximum_repetitions) {
            return repetitions;
        }
        repetitions*=2u;
    }
}

void print_configuration(
        const std::vector<NamedIntegrator>& integrators,
        const StepSizeType& step,
        Nat repetitions)
{
    std::cout << "# step=" << step << "\n";
    std::cout << "# repetitions=" << repetitions << "\n";
    std::cout << "# minimum_calibration_time_ms="
              << std::chrono::duration<double,std::milli>(
                     minimum_total_measurement_time).count()
              << "\n";

    for(const auto& named_integrator : integrators) {
        std::cout << "# " << named_integrator.key
                  << "=" << *named_integrator.integrator << "\n";
    }

    std::cout
        << "integrator\tfunction\trepetitions\ttime_us_per_step"
        << "\tremainder\tcomponent_remainders\n";
}

void run_benchmarks(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step,
        Nat repetitions)
{
    for(const auto& named_integrator : integrators) {
        for(const auto& benchmark_case : cases) {
            const auto measurement=measure(
                *named_integrator.integrator,
                benchmark_case,
                step,
                repetitions);

            const double time_us_per_step=
                std::chrono::duration<double,std::micro>(
                    measurement.elapsed).count()
                / static_cast<double>(repetitions);

            if(!(time_us_per_step>0.0) || !std::isfinite(time_us_per_step)) {
                std::cerr
                    << "Invalid timing for integrator "
                    << named_integrator.name
                    << " and function "
                    << benchmark_case.name
                    << ": "
                    << time_us_per_step
                    << " us\n";
                std::exit(2);
            }

            std::cout
                << named_integrator.key << "\t"
                << benchmark_case.name << "\t"
                << repetitions << "\t"
                << std::fixed << std::setprecision(6)
                << time_us_per_step << "\t"
                << measurement.flow.error() << "\t"
                << measurement.flow.errors()
                << "\n";
        }
    }
}

} // namespace

Int main(Int argc, const char* argv[]) {
    const std::string selection=(argc>1 ? argv[1] : "all");

    if(argc>2) {
        std::cerr
            << "Usage: benchmark_integrators_single_step "
            << "[all|integrator-key]\n";
        return 1;
    }

    const auto cases=make_cases();
    const auto integrators=select_integrators(selection);
    const StepSizeType step=determine_step(integrators,cases);
    const Nat repetitions=determine_repetitions(integrators,cases,step);

    print_configuration(integrators,step,repetitions);
    run_benchmarks(integrators,cases,step,repetitions);

    return 0;
}
