/***************************************************************************
 *            solving/interval_newton_method.hpp
 *
 *  Copyright  2006-20  Pieter Collins
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

/*! \file solving/interval_newton_method.hpp
 *  \brief Interval Newton method for validated equation solving.
 */

#ifndef ARIADNE_INTERVAL_NEWTON_METHOD_HPP
#define ARIADNE_INTERVAL_NEWTON_METHOD_HPP

#include "solving/equation_solver.hpp"

namespace Ariadne {

/*! \ingroup AlgebraicEquationSubModule
 *  \brief ExactIntervalType Newton solver. Uses the contractor \f$[x']=x_0-Df^{-1}([x])f(x_0)\f$.
 */
class IntervalNewtonSolver
    : public EquationSolverBase
{
  public:
    /*! \brief Constructor. */
    IntervalNewtonSolver(ApproximateDouble max_error, Nat max_steps) : EquationSolverBase(max_error,max_steps) { }
    IntervalNewtonSolver(MaximumError max_error, MaximumNumericTypeOfSteps max_steps) : EquationSolverBase(max_error,max_steps) { }
    /*! \brief Cloning operator. */
    virtual IntervalNewtonSolver* clone() const { return new IntervalNewtonSolver(*this); }
    /*! \brief Write to an output stream. */
    virtual Void _write(OutputStream& os) const;

    using EquationSolverBase::implicit;
  public:
    virtual ValidatedVectorMultivariateFunctionPatch implicit_step(const ValidatedVectorMultivariateFunction& f, const ValidatedVectorMultivariateFunctionPatch& p, const ValidatedVectorMultivariateFunctionPatch& x) const;

    virtual Vector<ValidatedNumericType> step(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& pt) const;
};



} // namespace Ariadne

#endif /* ARIADNE_INTERVAL_NEWTON_METHOD_HPP */
