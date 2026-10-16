/***************************************************************************
 *            solving/nonlinear_programming.hpp
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

/*! \file solving/nonlinear_programming.hpp
 *  \brief Nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_PROGRAMMING_HPP
#define ARIADNE_NONLINEAR_PROGRAMMING_HPP

#include "numeric/declarations.hpp"
#include "algebra/declarations.hpp"
#include "function/declarations.hpp"
#include "geometry/declarations.hpp"

#include "logging/logging.hpp"
#include "numeric/numeric.hpp"
#include "utility/tuple.hpp"


namespace Ariadne {

template<class X, class R> class Constraint;

class InfeasibleProblemException : public std::runtime_error {
  public: InfeasibleProblemException() : std::runtime_error("InfeasibleProblemException") { }
};
class IndeterminateFeasibilityException : public std::runtime_error {
  public: IndeterminateFeasibilityException() : std::runtime_error("IndeterminateFeasibilityException") { }
};
class DegenerateNonlinearFeasibilityProblemException : public std::runtime_error {
  public: DegenerateNonlinearFeasibilityProblemException() : std::runtime_error("DegenerateNonlinearFeasibilityProblemException") { }
};
class NearBoundaryOfFeasibleDomainException : public std::runtime_error {
  public: NearBoundaryOfFeasibleDomainException() : std::runtime_error("NearBoundaryOfFeasibleDomainException") { }
};

//! \ingroup OptimisationSubModule EvaluationModule
//! Interface for nonlinear programming solvers.
class NonlinearProgrammingSolverInterface {
    using FLT=FloatDP;
  public:
    typedef Bounds<FLT> ValidatedNumericType;
    typedef Vector<FLT> ExactVectorType;
    typedef Vector<Bounds<FLT>> ValidatedVectorType;
    typedef Vector<Approximation<FLT>> ApproximateVectorType;
  public:
    //! \brief Virtual destructor.
    virtual ~NonlinearProgrammingSolverInterface() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual NonlinearProgrammingSolverInterface* clone() const = 0;

    //! \brief Solve the general nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    //! \brief Solve the standard nonlinear programming problem \f$\min f(x) \text{ such that } x\in D,\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const = 0;

    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    //! \brief Tests is the standard nonlinear feasibility problem \f$x\in D,\ g(x)\leq 0 \text{ and } h(x) = 0\f$ is feasible. Assumes \f$D\f$ is bounded with nonempty interior.
    //! \internal This is one of the simplest nonlinear programming problems, and is a good test case for new algorithms.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const = 0;

    //! \brief Tests if the point \a x is feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual Bool almost_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                       FloatDPApproximationVector x, FloatDPApproximation eps) const = 0;
    //! \brief Tests if the point \a x is feasible.
    virtual Bool is_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                   FloatDPVector x) const = 0;
    //! \brief Tests if the point \a x is near feasible.
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                      ExactVectorType x) const = 0;
    //! \brief Tests if the point \a x is near feasible, using approximate multipliers \a y to guide the search.
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                      ExactVectorType x, ExactVectorType y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a x.
    virtual Bool validate_infeasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                      ExactVectorType x, ExactVectorType y) const = 0;
    //! \brief Tests if the box \a X definitely containss a feasible point.
    virtual ValidatedKleenean contains_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                            FloatDPBoundsVector X) const = 0;
    //! \brief Tests if the Lagrange multipliers \a y are a certificate of infeasiblity.
    virtual Bool is_infeasibility_certificate(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                              ExactVectorType y) const = 0;
};

//! \ingroup OptimisationSubModule
//! Common routines for nonlinear minimisation
class NonlinearProgrammingSolverBase
    : public NonlinearProgrammingSolverInterface
{
  protected:
    static const FloatDP zero;
    static const FloatDP one;
  public:
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const;

    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const;

    virtual Bool almost_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                       ApproximateVectorType x, FloatDPApproximation error) const;
    virtual Bool is_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                   ExactVectorType x) const;
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                      ExactVectorType x) const;
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                      ExactVectorType x, ExactVectorType y) const;
    virtual Bool validate_infeasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                        ExactVectorType x, ExactVectorType y) const;
    virtual ValidatedKleenean contains_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                            ValidatedVectorType X) const;
    virtual Bool is_infeasibility_certificate(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                              ExactVectorType lambda) const;
};

} // namespace Ariadne

#endif
