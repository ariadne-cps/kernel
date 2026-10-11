#include <chrono>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "io/command_line_interface.hpp"
#include "logging/logging.hpp"
#include "integrator_benchmark.hpp"

using namespace Ariadne;
using namespace IntegratorBenchmark;

namespace {

using Clock=std::chrono::steady_clock;

constexpr double benchmark_step_maximum_error=1e-8;
constexpr Nat maximum_halvings=20u;
constexpr Nat refinement_iterations=6u;

struct Trial {
    Bool accepted;
    double elapsed_ms;
    std::optional<FlowStepModelType> flow;
};

struct MaximumStepResult {
    std::string display;
};

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

Trial try_step(
        const NamedIntegrator& named_integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step)
{
    LOGGING_PRINTLN_AT(2,
        "Trial start: integrator="<<named_integrator.key
        <<", function="<<benchmark_case.name
        <<", step="<<step)

    const auto start=Clock::now();
    try {
        auto flow=named_integrator.integrator->flow_step(
            benchmark_case.vector_field,
            benchmark_case.domain,
            step);

        const double elapsed_ms=
            std::chrono::duration<double,std::milli>(
                Clock::now()-start).count();
        const auto remainder=flow.error();

        LOGGING_PRINTLN_AT(2,
            "Trial accepted: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", elapsed="<<elapsed_ms<<" ms"
            <<", remainder="<<remainder)

        return {true,elapsed_ms,std::move(flow)};
    } catch(const FlowTimeStepException& exception) {
        LOGGING_PRINTLN_AT(2,
            "Trial rejected: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", reason=FlowTimeStepException"
            <<", message="<<exception.what())
    } catch(const IncompleteFlowException& exception) {
        LOGGING_PRINTLN_AT(2,
            "Trial rejected: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", reason=IncompleteFlowException"
            <<", message="<<exception.what())
    } catch(const FlowBoundsException& exception) {
        LOGGING_PRINTLN_AT(2,
            "Trial rejected: integrator="<<named_integrator.key
            <<", function="<<benchmark_case.name
            <<", step="<<step
            <<", reason=FlowBoundsException"
            <<", message="<<exception.what())
    }

    const double elapsed_ms=
        std::chrono::duration<double,std::milli>(
            Clock::now()-start).count();

    return {false,elapsed_ms,std::nullopt};
}

MaximumStepResult benchmark_pair(
        const NamedIntegrator& named_integrator,
        const BenchmarkCase& benchmark_case)
{
    LOGGING_PRINTLN_AT(1,
        "Maximum-step search: integrator="<<named_integrator.key
        <<", function="<<benchmark_case.name)

    StepSizeType lower=0.0_x;
    StepSizeType upper=0.03125_x;
    Trial upper_trial=try_step(named_integrator,benchmark_case,upper);

    if(upper_trial.accepted) {
        lower=upper;
        for(Nat i=0u; i!=maximum_halvings; ++i) {
            StepSizeType candidate=upper+upper;
            Trial trial=try_step(
                named_integrator,benchmark_case,candidate);
            if(!trial.accepted) {
                upper=candidate;
                upper_trial=std::move(trial);
                break;
            }
            lower=candidate;
            upper=candidate;
            upper_trial=std::move(trial);
        }

        if(upper_trial.accepted) {
            std::ostringstream stream;
            stream << ">=" << lower;

            LOGGING_PRINTLN_AT(1,
                "Maximum-step search completed: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", result="<<stream.str())

            return {stream.str()};
        }
    } else {
        for(Nat i=0u; i!=maximum_halvings; ++i) {
            upper=hlf(upper);
            Trial trial=try_step(
                named_integrator,benchmark_case,upper);
            if(trial.accepted) {
                lower=upper;
                upper=upper+upper;
                break;
            }
        }

        if(lower==0.0_x) {
            std::ostringstream stream;
            stream << "<" << upper;

            LOGGING_PRINTLN_AT(1,
                "Maximum-step search completed: integrator="
                <<named_integrator.key
                <<", function="<<benchmark_case.name
                <<", result="<<stream.str())

            return {stream.str()};
        }
    }

    LOGGING_PRINTLN_AT(1,
        "Maximum-step bracket: integrator="<<named_integrator.key
        <<", function="<<benchmark_case.name
        <<", lower="<<lower
        <<", upper="<<upper)

    for(Nat i=0u; i!=refinement_iterations; ++i) {
        const StepSizeType middle=hlf(lower+upper);
        Trial trial=try_step(
            named_integrator,benchmark_case,middle);
        if(trial.accepted) {
            lower=middle;
        } else {
            upper=middle;
        }
    }

    LOGGING_PRINTLN_AT(1,
        "Maximum-step search completed: integrator="
        <<named_integrator.key
        <<", function="<<benchmark_case.name
        <<", accepted="<<lower
        <<", rejected="<<upper)

    std::ostringstream stream;
    stream << lower;
    return {stream.str()};
}

void print_table(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const std::vector<std::vector<MaximumStepResult>>& results)
{
    std::ostringstream table;
    table
        <<"Maximum-step integrator benchmark\n"
        <<"step_maximum_error="<<benchmark_step_maximum_error<<"\n"
        <<"Values are largest accepted steps found; "
        <<">= means no failing upper bound was found; "
        <<"< means no accepted step was found\n"
        <<"function";

    for(const auto& named_integrator : integrators) {
        table<<"\t"<<named_integrator.key;
    }
    table<<"\n";

    for(SizeType case_index=0u; case_index!=cases.size(); ++case_index) {
        table<<cases[case_index].name;
        for(SizeType integrator_index=0u;
                integrator_index!=integrators.size();
                ++integrator_index) {
            table
                <<"\t"
                <<results[case_index][integrator_index].display;
        }
        table<<"\n";
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
        "Integrator selection="<<selection
        <<", step maximum error="<<benchmark_step_maximum_error)

    const auto cases=make_cases();
    const auto integrators=select_integrators(
        selection,benchmark_step_maximum_error);

    std::vector<std::vector<MaximumStepResult>> results(
        cases.size(),
        std::vector<MaximumStepResult>(integrators.size()));

    for(SizeType integrator_index=0u;
            integrator_index!=integrators.size();
            ++integrator_index) {
        for(SizeType case_index=0u;
                case_index!=cases.size();
                ++case_index) {
            results[case_index][integrator_index]=benchmark_pair(
                integrators[integrator_index],
                cases[case_index]);
        }
    }

    print_table(integrators,cases,results);
    return 0;
}
