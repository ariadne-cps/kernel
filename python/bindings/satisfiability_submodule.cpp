/***************************************************************************
 *            satisfiability_submodule.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

#include "pybind11.hpp"
#include "kernel-utilities.hpp"

#include "satisfiability/smt_solver.hpp"

using namespace Ariadne;

Void satisfiability_submodule(pybind11::module& module)
{
    pybind11::enum_<SmtResultStatus>(module,"SmtResultStatus")
        .value("UNSAT",SmtResultStatus::UNSAT)
        .value("EPSILON_SAT",SmtResultStatus::EPSILON_SAT)
        .value("UNKNOWN",SmtResultStatus::UNKNOWN);

    pybind11::enum_<SmtUnknownReason>(module,"SmtUnknownReason")
        .value("NONE",SmtUnknownReason::NONE)
        .value("RESOURCE_EXHAUSTED",SmtUnknownReason::RESOURCE_EXHAUSTED)
        .value("DP_RESOLUTION_EXHAUSTED",SmtUnknownReason::DP_RESOLUTION_EXHAUSTED)
        .value("MIXED",SmtUnknownReason::MIXED);

    pybind11::class_<SmtSolverConfiguration> configuration_class(
        module,"SmtSolverConfiguration");
    configuration_class.def(pybind11::init<ExactDouble>());
    configuration_class.def("epsilon",&SmtSolverConfiguration::epsilon);
    configuration_class.def(
        "theory_minimization_budget",
        &SmtSolverConfiguration::theory_minimization_budget);
    configuration_class.def(
        "learned_clause_limit",
        &SmtSolverConfiguration::learned_clause_limit);
    configuration_class.def(
        "box_processing_limit",
        &SmtSolverConfiguration::box_processing_limit);

    pybind11::class_<SmtSearchStatistics> statistics_class(
        module,"SmtSearchStatistics");
    statistics_class.def_readonly(
        "boxes_processed",&SmtSearchStatistics::boxes_processed);
    statistics_class.def_readonly(
        "boxes_pruned",&SmtSearchStatistics::boxes_pruned);
    statistics_class.def_readonly(
        "boxes_split",&SmtSearchStatistics::boxes_split);
    statistics_class.def_readonly(
        "boolean_decisions",&SmtSearchStatistics::boolean_decisions);
    statistics_class.def_readonly(
        "theory_checks",&SmtSearchStatistics::theory_checks);

    pybind11::class_<SmtResult> result_class(module,"SmtResult");
    result_class.def("status",&SmtResult::status);
    result_class.def("is_unsat",&SmtResult::is_unsat);
    result_class.def("is_epsilon_sat",&SmtResult::is_epsilon_sat);
    result_class.def("is_unknown",&SmtResult::is_unknown);
    result_class.def("unknown_reason",&SmtResult::unknown_reason);
    result_class.def("has_witness",&SmtResult::has_witness);
    result_class.def(
        "witness",&SmtResult::witness,
        pybind11::return_value_policy::reference_internal);
    result_class.def(
        "statistics",&SmtResult::statistics,
        pybind11::return_value_policy::reference_internal);

    pybind11::class_<SmtSolver> solver_class(module,"SmtSolver");
    solver_class.def(pybind11::init<SmtSolverConfiguration>());
    solver_class.def(
        "solve",
        pybind11::overload_cast<
            RealSpace const&,
            ExactBoxType const&,
            ContinuousPredicate const&
        >(&SmtSolver::solve,pybind11::const_));
    solver_class.def(
        "solve_parallel",
        pybind11::overload_cast<
            RealSpace const&,
            ExactBoxType const&,
            ContinuousPredicate const&
        >(&SmtSolver::solve_parallel,pybind11::const_));
}
