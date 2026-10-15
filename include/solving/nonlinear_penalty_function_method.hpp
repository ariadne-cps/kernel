/***************************************************************************
 *            solving/nonlinear_penalty_function_method.hpp
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

/*! \file solving/nonlinear_penalty_function_method.hpp
 *  \brief Penalty-function methods for nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_PENALTY_FUNCTION_METHOD_HPP
#define ARIADNE_NONLINEAR_PENALTY_FUNCTION_METHOD_HPP

#include "solving/nonlinear_programming.hpp"

namespace Ariadne {

//! \ingroup OptimisationSubModule
//! \brief Solver for feasibility problems based on a penalty-function approach
//!
//! For the feasibility problem \f$x\in D,\ g(x)\in C,\ h(x)=0\f$ where \f$D\f$ is bounded and \f$D,C\f$ have nonempty interiors.
//! Introduce slack variable \f$w=g(x)\f$, and minimise \f[ \sum_{j} (g_j(x)-w_j)^2 + \sum_k h_k(x)^2 \f] with \f$x\in D\f$ and \f$w\in C\f$.
//! Since the minimiser is not unique, we need to add penalty terms \f$-\mu/2(\log(x_u-x)+log(x-x_l))\f$ and \f$-\nu/2(\log(w_u-w)+log(w-w_l))\f$.
//! It suffices to add a penalty in \f$x\f$.
class NonlinearPenaltyFunctionSolver
    : public NonlinearProgrammingSolverBase
{
  public:
    using NonlinearProgrammingSolverBase::minimise;
    using NonlinearProgrammingSolverBase::feasible;
  public:
    virtual NonlinearPenaltyFunctionSolver* clone() const;
    virtual ValidatedKleenean check_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ExactVectorType x, ExactVectorType y) const;
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;
    virtual Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                                  FloatDPApproximationVector& x, FloatDPApproximationVector& w, FloatDPApproximation& mu) const;
    virtual Void feasibility_step(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                                  FloatDPBoundsVector& x, FloatDPBoundsVector& w) const;
    virtual Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                                  FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& z) const;
};





} // namespace Ariadne

#endif
