/***************************************************************************
 *            solver_submodule.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
#include "kernel-utilities.hpp"

#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "solving/equation_solver_interface.hpp"
#include "solving/equation_solver.hpp"
#include "solving/interval_newton_method.hpp"
#include "solving/krawczyk_method.hpp"
#include "function/taylor_function.hpp"


using namespace Ariadne;

namespace Ariadne {

template<class T> Nat __hash__(const T&);
template<> Nat __hash__<FloatDPBoundsVector>(const FloatDPBoundsVector&) {
    return 0;
}

class EquationSolverWrapper
  : public pybind11::wrapper< EquationSolverInterface >
{
    typedef EquationSolverInterface::ValidatedNumericType ValidatedNumericType;
    typedef EquationSolverInterface::ApproximateNumericType ApproximateNumericType;
  public:
    EquationSolverInterface* clone() const { return this->get_override("clone")(); }
    Void set_maximum_error(ApproximateDouble me) { this->get_override("set_maximum_error")(me); }
    ExactDouble maximum_error() const { return this->get_override("maximum_error")(); }
    Void set_maximum_number_of_steps(Nat ns) { this->get_override("set_maximum_number_of_steps")(ns); }
    Nat maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    Vector<ValidatedNumericType> zero(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("zero")(f,bx); }
    Vector<ValidatedNumericType> fixed_point(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("fixed_point")(f,bx); }
    Vector<ValidatedNumericType> solve(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& pt) const {
        return this->get_override("solve")(f,pt); }
    Vector<ValidatedNumericType> solve(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve")(f,bx); }
    ValidatedVectorMultivariateFunctionPatch implicit(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& pd, const ExactBoxType& bx) const {
        return this->get_override("implicit")(f,pd,bx); }
    ValidatedScalarMultivariateFunctionPatch implicit(const ValidatedScalarMultivariateFunction& f, const ExactBoxType& pd, const ExactIntervalType& ivl) const {
        return this->get_override("implicit")(f,pd,ivl); }
    ValidatedVectorMultivariateFunctionPatch continuation(const ValidatedVectorMultivariateFunction& f, const Vector< ApproximateNumericType>& a, const ExactBoxType& X,  const ExactBoxType& A) const {
        return this->get_override("continuation")(f,a,X,A); }
    List< Vector<ValidatedNumericType> > list_solve_all(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        Set< Vector<ValidatedNumericType> > set = this->get_override("solve_all")(f,bx);
        return make_list(set);
    }
    Void _write(OutputStream& os) const { this->get_override("_write")(os); }
};




} // namespace Ariadne


Void export_solvers(pybind11::module& module)
{
    typedef EquationSolverInterface::ValidatedNumericType ValidatedNumericType;
    pybind11::class_<EquationSolverInterface,EquationSolverWrapper> solver_interface_class(module,"EquationSolverInterface");
    solver_interface_class.def("maximum_error", &EquationSolverInterface::maximum_error);
    solver_interface_class.def("set_maximum_error", &EquationSolverInterface::set_maximum_error);
    solver_interface_class.def("maximum_number_of_steps", &EquationSolverInterface::maximum_number_of_steps);
    solver_interface_class.def("set_maximum_number_of_steps", &EquationSolverInterface::set_maximum_number_of_steps);
    solver_interface_class.def("zero", &EquationSolverInterface::zero);
    solver_interface_class.def("fixed_point", &EquationSolverInterface::fixed_point);
    solver_interface_class.def("solve", (Vector<ValidatedNumericType>(EquationSolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &EquationSolverInterface::solve);
    solver_interface_class.def("implicit",(ValidatedVectorMultivariateFunctionPatch(EquationSolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&,const ExactBoxType&)const) &EquationSolverInterface::implicit);
    solver_interface_class.def("implicit",(ValidatedScalarMultivariateFunctionPatch(EquationSolverInterface::*)(const ValidatedScalarMultivariateFunction&,const ExactBoxType&,const ExactIntervalType&)const) &EquationSolverInterface::implicit);
    solver_interface_class.def("solve_all",[](EquationSolverInterface const& s, const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx){ return static_cast< List<Vector<ValidatedNumericType>> >( s.solve_all(f,bx) );});
    solver_interface_class.def("__str__",&__cstr__<EquationSolverInterface>);

    pybind11::class_<IntervalNewtonSolver, EquationSolverInterface> interval_newton_solver_class(module,"IntervalNewtonSolver");
    interval_newton_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

    pybind11::class_<KrawczykSolver, EquationSolverInterface> krawczyk_solver_class(module,"KrawczykSolver");
    krawczyk_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

    pybind11::class_<FactoredKrawczykSolver, EquationSolverInterface> factored_krawczyk_solver_class(module,"FactoredKrawczykSolver");
    factored_krawczyk_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

}





Void export_solving_optimization(pybind11::module& module);

Void solving_submodule(pybind11::module& module)
{
    export_solvers(module);
    export_solving_optimization(module);
}

/***************************************************************************
 *            optimization_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
#include "kernel-utilities.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "geometry/box.hpp"

#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "solving/linear_programming.hpp"
#include "solving/interior_point_method.hpp"
#include "solving/simplex_algorithm.hpp"
#include "solving/nonlinear_programming.hpp"
#include "solving/constraint_solver.hpp"

using namespace Ariadne;




Void export_slackness(pybind11::module& module)
{
    pybind11::enum_<Slackness> variable_enum(module,"Slackness");
    variable_enum.value("BASIS", Slackness::BASIS);
    variable_enum.value("LOWER", Slackness::LOWER);
    variable_enum.value("UPPER", Slackness::UPPER);
}

Void export_constraint(pybind11::module& module)
{
    pybind11::class_<EffectiveConstraint> effective_nonlinear_constraint_class(module,"EffectiveConstraint");
    effective_nonlinear_constraint_class.def(pybind11::init<EffectiveConstraint>());
    effective_nonlinear_constraint_class.def(pybind11::init<EffectiveNumber,EffectiveScalarMultivariateFunction,EffectiveNumber>());
    effective_nonlinear_constraint_class.def(pybind11::init([](EffectiveScalarMultivariateFunction const& f, EffectiveNumber const& u){return EffectiveConstraint(EffectiveNumber(-infty),f,u);}));
    effective_nonlinear_constraint_class.def(pybind11::init([](EffectiveNumber const& l, EffectiveScalarMultivariateFunction const& f){return EffectiveConstraint(l,f,EffectiveNumber(infty));}));
    effective_nonlinear_constraint_class.def(pybind11::init<EffectiveScalarMultivariateFunction,EffectiveNumber>());

    effective_nonlinear_constraint_class.def("__str__",&__cstr__<EffectiveConstraint>);

    pybind11::class_<ValidatedConstraint> validated_nonlinear_constraint_class(module,"ValidatedConstraint");
    validated_nonlinear_constraint_class.def(pybind11::init<ValidatedNumber,ValidatedScalarMultivariateFunction,ValidatedNumber
    >());
    validated_nonlinear_constraint_class.def(pybind11::init<ValidatedConstraint>());
    validated_nonlinear_constraint_class.def(pybind11::init<EffectiveConstraint>());
    validated_nonlinear_constraint_class.def("lower_bound", &ValidatedConstraint::lower_bound);
    validated_nonlinear_constraint_class.def("upper_bound", &ValidatedConstraint::upper_bound);
    validated_nonlinear_constraint_class.def("function", (const ValidatedScalarMultivariateFunction&(ValidatedConstraint::*)()const) &ValidatedConstraint::function);
    validated_nonlinear_constraint_class.def("__str__",&__cstr__<ValidatedConstraint>);
}

Void export_optimiser_interface(pybind11::module& module)
{
    typedef OptimiserInterface::ValidatedNumericType ValidatedNumericType;
    typedef OptimiserInterface::ExactVectorType ExactVectorType;

    pybind11::class_<OptimiserInterface> optimiser_interface_class(module,"OptimiserInterface");
    optimiser_interface_class.def("minimise", (Vector<ValidatedNumericType>(OptimiserInterface::*)(ValidatedScalarMultivariateFunction, ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType)const) &OptimiserInterface::minimise);
    optimiser_interface_class.def("minimise", (Vector<ValidatedNumericType>(OptimiserInterface::*)(ValidatedScalarMultivariateFunction, ExactBoxType, ValidatedVectorMultivariateFunction, ValidatedVectorMultivariateFunction)const) &OptimiserInterface::minimise);
    optimiser_interface_class.def("feasible", (ValidatedKleenean(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType)const) &OptimiserInterface::feasible);
    optimiser_interface_class.def("feasible", (ValidatedKleenean(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ValidatedVectorMultivariateFunction)const) &OptimiserInterface::feasible);
    optimiser_interface_class.def("almost_feasible_point", (Bool(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,FloatDPApproximationVector, FloatDPApproximation)const) &OptimiserInterface::almost_feasible_point);
    optimiser_interface_class.def("is_feasible_point", (Bool(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,FloatDPVector)const) &OptimiserInterface::is_feasible_point);
    optimiser_interface_class.def("validate_feasibility", (Bool(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,ExactVectorType)const) &OptimiserInterface::validate_feasibility);
    optimiser_interface_class.def("validate_feasibility", (Bool(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,ExactVectorType,ExactVectorType)const) &OptimiserInterface::validate_feasibility);
    optimiser_interface_class.def("validate_infeasibility", (Bool(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,ExactVectorType,ExactVectorType)const) &OptimiserInterface::validate_infeasibility);
    optimiser_interface_class.def("contains_feasible_point", (ValidatedKleenean(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,FloatDPBoundsVector)const) &OptimiserInterface::contains_feasible_point);
    optimiser_interface_class.def("is_infeasibility_certificate", (Bool(OptimiserInterface::*)(ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType,ExactVectorType)const) &OptimiserInterface::is_infeasibility_certificate);
    //NOTE: Not in C++ API
    //optimiser_interface_class.def("__str__", &__cstr__<OptimiserInterface>);
}

Void export_interior_point_solvers(pybind11::module& module)
{
    pybind11::class_<NonlinearInfeasibleInteriorPointOptimiser,OptimiserInterface> nonlinear_infeasible_interior_point_solver_class(module,"NonlinearInfeasibleInteriorPointOptimiser");
    nonlinear_infeasible_interior_point_solver_class.def(pybind11::init<>());

    pybind11::class_<NonlinearInteriorPointOptimiser,OptimiserInterface> nonlinear_interior_point_solver_class(module,"NonlinearInteriorPointOptimiser");
    nonlinear_interior_point_solver_class.def(pybind11::init<>());
}

Void export_interior_point_solver(pybind11::module& module)
{
    pybind11::class_<InteriorPointSolver> interior_point_solver_class(module,"InteriorPointSolver");
    interior_point_solver_class.def(pybind11::init<>());
    interior_point_solver_class.def("minimise", &InteriorPointSolver::minimise);
    interior_point_solver_class.def("feasible", (ValidatedKleenean(InteriorPointSolver::*)(const Vector<FloatDP>&,const Vector<FloatDP>&, const Matrix<FloatDP>&,const Vector<FloatDP>&)const) &InteriorPointSolver::feasible);
    interior_point_solver_class.def("validate_feasibility", &InteriorPointSolver::validate_feasibility);
}


Void export_constraint_solver(pybind11::module& module)
{
    pybind11::class_<ConstraintSolver> constraint_solver_class(module,"ConstraintSolver");
    constraint_solver_class.def(pybind11::init<>());
    constraint_solver_class.def("feasible", (Pair<ValidatedKleenean,ConstraintSolver::ExactPointType>(ConstraintSolver::*)(const ExactBoxType&,const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &ConstraintSolver::feasible);
    constraint_solver_class.def("check_feasibility", &ConstraintSolver::check_feasibility);
    constraint_solver_class.def("reduce_box", [](ConstraintSolver const& solver, ExactBoxType const& domain, ValidatedVectorMultivariateFunction const& function, ExactBoxType const& codomain) {
        UpperBoxType reduced(domain);
        Bool empty = solver.reduce(reduced,function,codomain);
        return pybind11::make_tuple(empty,reduced);
    });
    constraint_solver_class.def("hull_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&)const) &ConstraintSolver::hull_reduce);
    constraint_solver_class.def("box_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&,SizeType)const) &ConstraintSolver::box_reduce);
    constraint_solver_class.def("monotone_reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedScalarMultivariateFunction&,const ExactIntervalType&,SizeType)const) &ConstraintSolver::monotone_reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const List<ValidatedConstraint>&)const) &ConstraintSolver::reduce);
    constraint_solver_class.def("reduce", (Bool(ConstraintSolver::*)(UpperBoxType&,const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &ConstraintSolver::reduce);
}


template<class X>
Void export_simplex_solver(pybind11::module& module)
{
    typedef Array<SizeType> SizeArray;

    typedef RigorousNumericType<X> XX;

    pybind11::class_< SimplexSolver<X> > simplex_solver_class(module,(class_name<X>()+"SimplexSolver").c_str());
    simplex_solver_class.def(pybind11::init<>());
    simplex_solver_class.def("lpstep",(Bool(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,Array<Slackness>& ,SizeArray&,Matrix<XX>&,Vector<XX>&)const) &SimplexSolver<X>::lpstep);


    simplex_solver_class.def("minimise",(Vector<XX>(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&)const) &SimplexSolver<X>::minimise);
    simplex_solver_class.def("feasible",(ValidatedKleenean(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&)const) &SimplexSolver<X>::feasible);

    simplex_solver_class.def("verify_feasibility",(ValidatedKleenean(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,const Array<Slackness>&)const) &SimplexSolver<X>::verify_feasibility);

    simplex_solver_class.def("compute_basis",(Pair< SizeArray, Matrix<XX> >(SimplexSolver<X>::*)(const Matrix<X>&)const) &SimplexSolver<X>::compute_basis);

}


Void export_solving_optimization(pybind11::module& module) {
    export_slackness(module);
    export_constraint(module);
    export_optimiser_interface(module);
    export_interior_point_solvers(module);
    export_simplex_solver<FloatDPApproximation>(module);
    export_simplex_solver<FloatDP>(module);
    export_simplex_solver<Rational>(module);
    export_interior_point_solver(module);
    export_constraint_solver(module);
}
