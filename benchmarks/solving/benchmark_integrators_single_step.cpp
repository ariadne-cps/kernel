#include <chrono>
#include <cmath>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "io/command_line_interface.hpp"
#include "logging/logging.hpp"
#include "integrator_benchmark.hpp"

using namespace Ariadne;
using namespace IntegratorBenchmark;

namespace {

using Clock = std::chrono::steady_clock;
using Duration = std::chrono::nanoseconds;

constexpr auto minimum_total_measurement_time=std::chrono::milliseconds(1);
constexpr Nat maximum_repetitions=1024u;
constexpr double benchmark_step_maximum_error=1e-5;

struct Measurement {
    Duration elapsed;
    FlowStepModelType flow;
};

FlowStepModelType run_integrator_step(
        const IntegratorInterface& integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step)
{
    std::optional<FlowStepModelType> flow;
    try {
        LOGGING_RUN_AT(2,
            flow=integrator.flow_step(
                benchmark_case.vector_field,
                benchmark_case.domain,
                step))
    } catch(...) {
        Logger::instance().decrease_level(2);
        throw;
    }
    return std::move(flow.value());
}

struct BenchmarkResult {
    std::string integrator;
    std::string function;
    Nat repetitions;
    double time_us_per_step;
    decltype(std::declval<FlowStepModelType>().error()) remainder;
    decltype(std::declval<FlowStepModelType>().errors()) component_remainders;
};

StepSizeType benchmark_step() {
    return 0.03125_x;
}

Bool is_integrator_selection(const std::string& value) {
    return value=="all"
        || value=="taylor-picard"
        || value=="graded-taylor-picard"
        || value=="taylor-series"
        || value=="taylor-series-bounder"
        || value=="graded-taylor-series"
        || value=="preconditioned-graded-taylor-series";
}

std::string acquire_integrator_selection(
        Int argc,
        const char* argv[],
        List<String>& cli_arguments)
{
    std::string selection="all";
    Bool selection_found=false;

    cli_arguments.append(String(argv[0]));
    for(Int i=1; i!=argc; ++i) {
        const std::string argument(argv[i]);
        if(is_integrator_selection(argument)) {
            if(selection_found) {
                throw std::runtime_error(
                    "Integrator selection specified more than once.");
            }
            selection=argument;
            selection_found=true;
        } else {
            cli_arguments.append(String(argv[i]));
        }
    }
    return selection;
}

Measurement measure(
        const IntegratorInterface& integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step,
        Nat repetitions)
{
    FlowStepModelType flow=run_integrator_step(
        integrator,
        benchmark_case,
        step);

    const auto start=Clock::now();
    for(Nat i=0u; i!=repetitions; ++i) {
        flow=run_integrator_step(
            integrator,
            benchmark_case,
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
    const double minimum_measurement_time_ms=
        std::chrono::duration<double,std::milli>(
            minimum_total_measurement_time).count();

    LOGGING_PRINTLN_AT(1,
        "Calibrating repetitions with step="<<step
        <<", minimum measurement time="
        <<minimum_measurement_time_ms
        <<" ms")

    for(const auto& named_integrator : integrators) {
        for(const auto& benchmark_case : cases) {
            LOGGING_PRINTLN_AT(2,
                "Calibration start: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name)

            const auto measurement=measure(
                *named_integrator.integrator,
                benchmark_case,
                step,
                1u);

            const double elapsed_ms=
                std::chrono::duration<double,std::milli>(
                    measurement.elapsed).count();

            LOGGING_PRINTLN_AT(1,
                "Calibration: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", elapsed="<<elapsed_ms<<" ms")

            LOGGING_PRINTLN_AT(2,
                "Calibration remainder: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", remainder="<<measurement.flow.error()
                <<", components="<<measurement.flow.errors())

            if(measurement.elapsed<minimum_elapsed) {
                minimum_elapsed=measurement.elapsed;
            }
        }
    }

    if(minimum_elapsed<=Duration::zero()) {
        LOGGING_PRINTLN_AT(1,
            "Calibration selected maximum repetitions="
            <<maximum_repetitions
            <<" because the minimum measured duration was non-positive")
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

    LOGGING_PRINTLN_AT(1,
        "Calibration selected repetitions="<<repetitions)

    return repetitions;
}

std::vector<BenchmarkResult> run_benchmarks(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const StepSizeType& step,
        Nat repetitions)
{
    std::vector<BenchmarkResult> results;
    results.reserve(integrators.size()*cases.size());

    LOGGING_PRINTLN_AT(1,
        "Running benchmark with step="<<step
        <<", step maximum error="<<benchmark_step_maximum_error
        <<", repetitions="<<repetitions)

    for(const auto& named_integrator : integrators) {
        for(const auto& benchmark_case : cases) {
            LOGGING_PRINTLN_AT(2,
                "Benchmark start: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name)

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
                LOGGING_PRINTLN(
                    "Invalid timing for integrator "<<named_integrator.name
                    <<" and function "<<benchmark_case.name
                    <<": "<<time_us_per_step<<" us")
                throw std::runtime_error("Invalid benchmark timing.");
            }

            LOGGING_PRINTLN_AT(1,
                "Benchmark result: integrator="<<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", time="<<time_us_per_step<<" us/step"
                <<", remainder="<<measurement.flow.error())

            LOGGING_PRINTLN_AT(2,
                "Benchmark component remainders: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", components="<<measurement.flow.errors())

            results.push_back({
                named_integrator.key,
                benchmark_case.name,
                repetitions,
                time_us_per_step,
                measurement.flow.error(),
                measurement.flow.errors()
            });
        }
    }

    return results;
}

void print_results(
        const StepSizeType& step,
        Nat repetitions,
        const std::vector<BenchmarkResult>& results)
{
    std::ostringstream table;
    table
        <<"Single-step integrator benchmark\n"
        <<"step="<<step
        <<", step_maximum_error="<<benchmark_step_maximum_error
        <<", repetitions="<<repetitions<<"\n"
        <<"integrator\tfunction\trepetitions\ttime_us_per_step"
        <<"\tremainder\tcomponent_remainders\n";

    table<<std::fixed<<std::setprecision(6);
    for(const auto& result : results) {
        table
            <<result.integrator<<"\t"
            <<result.function<<"\t"
            <<result.repetitions<<"\t"
            <<result.time_us_per_step<<"\t"
            <<result.remainder<<"\t"
            <<result.component_remainders
            <<"\n";
    }

    LOGGING_PRINTLN(table.str())
}

} // namespace

Int main(Int argc, const char* argv[]) {
    Logger::instance().configuration().set_verbosity(1u);

    List<String> cli_arguments;
    const std::string selection=
        acquire_integrator_selection(argc,argv,cli_arguments);

    if(not CommandLineInterface::instance().acquire(cli_arguments)) {
        return -1;
    }

    LOGGING_PRINTLN_AT(1,
        "Integrator selection="<<selection)

    const auto cases=make_cases();
    const auto integrators=select_integrators(
        selection,benchmark_step_maximum_error);
    const StepSizeType step=benchmark_step();
    const Nat repetitions=determine_repetitions(integrators,cases,step);
    const auto results=run_benchmarks(
        integrators,cases,step,repetitions);

    print_results(step,repetitions,results);
    return 0;
}
