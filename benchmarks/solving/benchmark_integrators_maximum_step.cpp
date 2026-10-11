#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

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

Trial try_step(
        const NamedIntegrator& named_integrator,
        const BenchmarkCase& benchmark_case,
        const StepSizeType& step)
{
    std::cerr
        << "[maximum-step] integrator=" << named_integrator.key
        << " function=" << benchmark_case.name
        << " step=" << step
        << " START"
        << std::endl;

    const auto start=Clock::now();
    try {
        auto flow=named_integrator.integrator->flow_step(
            benchmark_case.vector_field,
            benchmark_case.domain,
            step);
        const double elapsed_ms=
            std::chrono::duration<double,std::milli>(
                Clock::now()-start).count();

        std::cerr
            << "[maximum-step] integrator=" << named_integrator.key
            << " function=" << benchmark_case.name
            << " step=" << step
            << " OK elapsed_ms=" << elapsed_ms
            << " remainder=" << flow.error()
            << std::endl;

        return {true,elapsed_ms,std::move(flow)};
    } catch(const FlowTimeStepException& exception) {
        std::cerr
            << "[maximum-step] integrator=" << named_integrator.key
            << " function=" << benchmark_case.name
            << " step=" << step
            << " FAIL FlowTimeStepException: " << exception.what()
            << std::endl;
    } catch(const IncompleteFlowException& exception) {
        std::cerr
            << "[maximum-step] integrator=" << named_integrator.key
            << " function=" << benchmark_case.name
            << " step=" << step
            << " FAIL IncompleteFlowException: " << exception.what()
            << std::endl;
    } catch(const FlowBoundsException& exception) {
        std::cerr
            << "[maximum-step] integrator=" << named_integrator.key
            << " function=" << benchmark_case.name
            << " step=" << step
            << " FAIL FlowBoundsException: " << exception.what()
            << std::endl;
    }

    return {
        false,
        std::chrono::duration<double,std::milli>(Clock::now()-start).count(),
        std::nullopt
    };
}

MaximumStepResult benchmark_pair(
        const NamedIntegrator& named_integrator,
        const BenchmarkCase& benchmark_case)
{
    StepSizeType lower=0.0_x;
    StepSizeType upper=0.03125_x;
    Trial upper_trial=try_step(named_integrator,benchmark_case,upper);

    if(upper_trial.accepted) {
        lower=upper;
        for(Nat i=0u; i!=maximum_halvings; ++i) {
            StepSizeType candidate=upper+upper;
            Trial trial=try_step(named_integrator,benchmark_case,candidate);
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
            return {stream.str()};
        }
    } else {
        for(Nat i=0u; i!=maximum_halvings; ++i) {
            upper=hlf(upper);
            Trial trial=try_step(named_integrator,benchmark_case,upper);
            if(trial.accepted) {
                lower=upper;
                upper=upper+upper;
                break;
            }
        }

        if(lower==0.0_x) {
            std::ostringstream stream;
            stream << "<" << upper;
            return {stream.str()};
        }
    }

    for(Nat i=0u; i!=refinement_iterations; ++i) {
        const StepSizeType middle=hlf(lower+upper);
        Trial trial=try_step(named_integrator,benchmark_case,middle);
        if(trial.accepted) {
            lower=middle;
        } else {
            upper=middle;
        }
    }

    std::cerr
        << "[maximum-step] integrator=" << named_integrator.key
        << " function=" << benchmark_case.name
        << " bracket=[" << lower << "," << upper << ")"
        << std::endl;

    std::ostringstream stream;
    stream << lower;
    return {stream.str()};
}

void print_table(
        const std::vector<NamedIntegrator>& integrators,
        const std::vector<BenchmarkCase>& cases,
        const std::vector<std::vector<MaximumStepResult>>& results)
{
    std::cout
        << "# benchmark_step_maximum_error="
        << benchmark_step_maximum_error << "\n"
        << "# table values are largest accepted steps found; "
        << ">= means no failing upper bound was found; "
        << "< means no accepted step was found\n";

    std::cout << "function";
    for(const auto& named_integrator : integrators) {
        std::cout << "\t" << named_integrator.key;
    }
    std::cout << "\n";

    for(SizeType case_index=0u; case_index!=cases.size(); ++case_index) {
        std::cout << cases[case_index].name;
        for(SizeType integrator_index=0u;
                integrator_index!=integrators.size();
                ++integrator_index) {
            std::cout
                << "\t"
                << results[case_index][integrator_index].display;
        }
        std::cout << "\n";
    }
}

} // namespace

Int main(Int argc, const char* argv[]) {
    const std::string selection=(argc>1 ? argv[1] : "all");
    if(argc>2) {
        std::cerr
            << "Usage: benchmark_integrators_maximum_step "
            << "[all|integrator-key]\n";
        return 1;
    }

    const auto cases=make_cases();
    const auto integrators=select_integrators(
        selection,benchmark_step_maximum_error);

    std::vector<std::vector<MaximumStepResult>> results(
        cases.size(),
        std::vector<MaximumStepResult>(integrators.size()));

    for(SizeType integrator_index=0u;
            integrator_index!=integrators.size();
            ++integrator_index) {
        for(SizeType case_index=0u; case_index!=cases.size(); ++case_index) {
            results[case_index][integrator_index]=benchmark_pair(
                integrators[integrator_index],
                cases[case_index]);
        }
    }

    print_table(integrators,cases,results);
    return 0;
}
