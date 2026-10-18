/***************************************************************************
 *            solving/nonlinear_interior_point_method.hpp
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

/*! \file solving/nonlinear_interior_point_method.hpp
 *  \brief Interior-point methods for nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_INTERIOR_POINT_METHOD_HPP
#define ARIADNE_NONLINEAR_INTERIOR_POINT_METHOD_HPP

#include "solving/nonlinear_programming.hpp"

namespace Ariadne {

//! \ingroup OptimisationSubModule
//! \brief Solver for nonlinear programming problems using invalid interior point methods.
//! \details Introduces variables \f$w\f$ and attempts to find \f$x\in D\f$ and \f$w\in C\f$ such that \f$g(x)=w\f$.
//!   The dual variables \f$y\f$ are unconstrained Lagrange multipliers for \f$y\cdot(g(x)-w)=0\f$.
class NonlinearInfeasibleInteriorPointSolver
    : public NonlinearProgrammingSolverBase
{
  public:
    virtual NonlinearInfeasibleInteriorPointSolver* clone() const { return new NonlinearInfeasibleInteriorPointSolver(*this); }

    using NonlinearProgrammingSolverBase::minimise;
    using NonlinearProgrammingSolverBase::feasible;

    struct PrimalDualData;
    struct StepData;

    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;

    //! \brief Runs the nonlinear feasibility iteration and returns both its validated status and final numerical candidate.
    //! \details The candidate is heuristic data only. It is guaranteed to lie in \f$D\f$,
    //! but a caller must validate it independently before using it as a certificate.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_candidate(
        ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall primal and dual variables.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, const PrimalDualData& wxy0) const;

    Void setup_feasibility(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                           StepData& stp) const;
    Void step(const ApproximateScalarMultivariateFunction& f, const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
              StepData& stp) const;
//    FloatDPApproximation compute_mu(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
//                     FloatDPApproximationVector& w, FloatDPApproximationVector& x, FloatDPApproximationVector& y) const;
};


//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using interior point methods.
class NonlinearInteriorPointSolver
    : public NonlinearProgrammingSolverBase
{
  public:
    virtual NonlinearInteriorPointSolver* clone() const { return new NonlinearInteriorPointSolver(*this); }

    using NonlinearProgrammingSolverBase::minimise;
    using NonlinearProgrammingSolverBase::feasible;

    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal \f$x_0\f$ and dual variables \f$y_0\f$ overall constraint violation \f$t_0\f$.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                                  const FloatDPApproximationVector& x0, const FloatDPApproximationVector& y0, const FloatDPApproximation& t0) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal variables \f$x_0\f$ and dual variables \f$y_0\f$.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                                                  const FloatDPApproximationVector& x0, const FloatDPApproximationVector& y0) const;

    Void minimisation_step(const ApproximateScalarMultivariateFunction& f, const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C, const ApproximateVectorMultivariateFunction& h,
                           FloatDPApproximationVector& x, FloatDPApproximationVector& w, FloatDPApproximationVector& kappa, FloatDPApproximationVector& lambda, const FloatDPApproximation& mu) const;
    Void setup_feasibility(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                           FloatDPApproximationVector& x, FloatDPApproximationVector& lambda) const;
    Void setup_feasibility(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                           FloatDPApproximationVector& x, FloatDPApproximationVector& lambda, FloatDPApproximation& t) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                          FloatDPApproximationVector& x, FloatDPApproximationVector& lambda) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                          FloatDPApproximationVector& x, FloatDPApproximationVector& lambda, FloatDPApproximation& violation) const;
    Void initialise_lagrange_multipliers(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                                         const FloatDPApproximationVector& x, FloatDPApproximationVector& lambda) const;
    FloatDPApproximation compute_mu(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                     const FloatDPApproximationVector& x, const FloatDPApproximationVector& lambda) const;

  private:
    FloatDPApproximation compute_mu(const ApproximateScalarMultivariateFunction& f, const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                     const FloatDPApproximationVector& x, const FloatDPApproximationVector& y) const;
    Void compute_violation(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                           FloatDPApproximationVector& x, FloatDPApproximation& t) const;
};




class IntervalNonlinearInteriorPointSolver
    : public NonlinearInteriorPointSolver
{
    virtual IntervalNonlinearInteriorPointSolver* clone() const { return new IntervalNonlinearInteriorPointSolver(*this); }
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    Void feasibility_step(const FloatDPVector& xl, const FloatDPVector& xu, const ValidatedVectorMultivariateFunction& h,
                          FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const;
};


class ApproximateNonlinearInteriorPointSolver
    : public NonlinearInteriorPointSolver
{
    virtual ApproximateNonlinearInteriorPointSolver* clone() const { return new ApproximateNonlinearInteriorPointSolver(*this); }
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                          FloatDPApproximationVector& X, FloatDPApproximationVector& Lambda) const;
};



} // namespace Ariadne

#endif
