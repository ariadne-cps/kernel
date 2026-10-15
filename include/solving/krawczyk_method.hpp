/***************************************************************************
 *            solving/krawczyk_method.hpp
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

/*! \file solving/krawczyk_method.hpp
 *  \brief Krawczyk methods for validated equation solving.
 */

#ifndef ARIADNE_KRAWCZYK_METHOD_HPP
#define ARIADNE_KRAWCZYK_METHOD_HPP

#include "solving/equation_solver.hpp"

namespace Ariadne {

/*! \ingroup AlgebraicEquationSubModule
 *  \brief Krawczyk solver. Uses the contractor \f$[x']=x_0-Mf(x_0)+(I-MDf([x])([x]-x_0)\f$
 *  where \f$M\f$ is typically \f$Df^{-1}(x_0)\f$.
 */
class KrawczykSolver
    : public EquationSolverBase
{
  public:
    /*! \brief Constructor. */
    KrawczykSolver(ApproximateDouble max_error, Nat max_steps) : EquationSolverBase(max_error,max_steps) { }
    KrawczykSolver(MaximumError max_error, MaximumNumericTypeOfSteps max_steps) : EquationSolverBase(max_error,max_steps) { }
    /*! \brief Cloning operator. */
    virtual KrawczykSolver* clone() const { return new KrawczykSolver(*this); }
    /*! \brief Write to an output stream. */
    virtual Void _write(OutputStream& os) const;

    /*! \brief Solve \f$f(a,x)=0\f$ for a in \a par, looking for solutions with x in \a ix. */
    virtual ValidatedVectorMultivariateFunctionPatch implicit_step(const ValidatedVectorMultivariateFunction& f, const ValidatedVectorMultivariateFunctionPatch& p, const ValidatedVectorMultivariateFunctionPatch& x) const;
  public:
    /*! \brief A single step of the Krawczyk contractor. */
    virtual Vector<ValidatedNumericType>
    step(const ValidatedVectorMultivariateFunction& f,
          const Vector<ValidatedNumericType>& pt) const;
};


/*! \ingroup Solvers
 *  \brief Modified Krawczyk solver. Uses the contractor \f$[x']=x_0-J^{-1}(f(x_0)+(J-Df([x])([x]-x_0)\f$
 *  where \f$J\f$ is typically \f$Df(x_0)\f$.
 */
class FactoredKrawczykSolver
    : public KrawczykSolver
{
  public:
    /*! \brief Constructor. */
    FactoredKrawczykSolver(ApproximateDouble max_error, Nat max_steps) : KrawczykSolver(max_error,max_steps) { }
    FactoredKrawczykSolver(MaximumError max_error, MaximumNumericTypeOfSteps max_steps) : KrawczykSolver(max_error,max_steps) { }
    /*! \brief Cloning operator. */
    virtual FactoredKrawczykSolver* clone() const { return new FactoredKrawczykSolver(*this); }
    /*! \brief Write to an output stream. */
    virtual Void _write(OutputStream& os) const;
  public:
    /*! \brief A single step of the modified Krawczyk contractor. */
    virtual Vector<ValidatedNumericType>
    step(const ValidatedVectorMultivariateFunction& f,
          const Vector<ValidatedNumericType>& pt) const;
};


} // namespace Ariadne

#endif /* ARIADNE_KRAWCZYK_METHOD_HPP */
