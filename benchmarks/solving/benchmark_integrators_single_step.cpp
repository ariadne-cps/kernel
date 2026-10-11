#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

#include "integrator_benchmark.hpp"

using namespace Ariadne;
using namespace IntegratorBenchmark;

namespace {

using Clock = std::chrono::steady_clock;
using Duration = std::chrono::nanoseconds;

constexpr auto minimum_total_measurement_time=std::chrono::milliseconds(1);
constexpr Nat maximum_repetitions=1024u;
constexpr double benchmark_step_maximum_error=1e-2;

struct Measurement {
    Duration elapsed;
    FlowStepModelType flow;
};

StepSizeType benchmark_step() {
    return 0.03125_x;
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
    Duration minimum_elapsed=Duration::max();

    for(const auto& named_integrator : integrators) {
        for(const auto& benchmark_case : cases) {
            std::cerr
                << "[repetition-calibration] integrator="
                << named_integrator.key
                << " function=" << benchmark_case.name
                << " START"
                << std::endl;

            const auto measurement=measure(
                *named_integrator.integrator,
                benchmark_case,
                step,
                1u);

            std::cerr
                << "[repetition-calibration] integrator="
                << named_integrator.key
                << " function=" << benchmark_case.name
                << " elapsed_ms="
                << std::chrono::duration<double,std::milli>(
                    measurement.elapsed).count()
                << " remainder=" << measurement.flow.error()
                << std::endl;

            if(measurement.elapsed<minimum_elapsed) {
                minimum_elapsed=measurement.elapsed;
            }
        }
    }

    if(minimum_elapsed<=Duration::zero()) {
        return maximum_repetitions;
    }

    const double required=
        std::chrono::duration<double>(minimum_total_measurement_time).count()
        / std::chrono::duration<double>(minimum_elapsed).count();

    Nat repetitions=static_cast<Nat>(std::ceil(required));
    if(repetitions<1u) {
        repetitions=1u;
    }
    if(repetitions>maximum_repetitions) {
        repetitions=maximum_repetitions;
    }
    return repetitions;
}

void run_benchmarks(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step,
        Nat repetitions)
{
    std::cout
        << "# step=" << step << "\n"
        << "# benchmark_step_maximum_error=" << benchmark_step_maximum_error << "\n"
        << "# repetitions=" << repetitions << "\n"
        << "integrator\tfunction\trepetitions\ttime_us_per_step"
        << "\tremainder\tcomponent_remainders\n";

    for(const auto& named_integrator : integrators) {
        for(const auto& benchmark_case : cases) {
            std::cerr
                << "[benchmark] integrator=" << named_integrator.key
                << " function=" << benchmark_case.name
                << " START"
                << std::endl;

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
                    << ": " << time_us_per_step << " us\n";
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

            std::cerr
                << "[benchmark] integrator=" << named_integrator.key
                << " function=" << benchmark_case.name
                << " OK time_us_per_step=" << time_us_per_step
                << " remainder=" << measurement.flow.error()
                << std::endl;
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
    const auto integrators=select_integrators(selection,benchmark_step_maximum_error);
    const StepSizeType step=benchmark_step();
    const Nat repetitions=determine_repetitions(integrators,cases,step);

    run_benchmarks(integrators,cases,step,repetitions);
    return 0;
}
