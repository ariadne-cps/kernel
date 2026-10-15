/***************************************************************************
 *            solving/linear_interior_point.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*! \file solving/linear_interior_point.hpp
 *  \brief Interior-point algorithms for linear programming.
 */

#ifndef ARIADNE_LINEAR_INTERIOR_POINT_HPP
#define ARIADNE_LINEAR_INTERIOR_POINT_HPP

#include "logging/logging.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "utility/tuple.hpp"
#include "solving/linear_programming.hpp"

namespace Ariadne {

enum class LinearProgramStatus : std::uint8_t { INDETERMINATE_FEASIBILITY=0, PRIMAL_FEASIBLE=1, DUAL_FEASIBLE=2, PRIMAL_DUAL_FEASIBLE=3, DEGENERATE_FEASIBILITY=4};
OutputStream& operator<<(OutputStream& os, LinearProgramStatus lps);

//! \ingroup OptimisationSubModule
//! Solver for linear programming problems using interior point methods.
class InteriorPointSolver
{
    typedef FloatDP X;
    typedef RigorousNumericType<X> VX;
    typedef ApproximateNumericType<X> AX;
  public:
    //! \brief Find approximate optimal solution of \f$\min c^T x \text{ s.t. } Ax=b; x\geq0\f$.
    //! Returns the tuple (v,x,y) where v is the optimal value, x is the optimal point,
    //! and y the corresponding dual feasible point.
    Tuple< AX, Vector<AX>, Vector<AX> >
    minimise(const Vector<X>& c, const Vector<X>& xl, const Vector<X>& xu, const Matrix<X>& A, const Vector<X>& b) const;

    //! \brief Find approximate optimal solution of \f$\min c^T x \text{ s.t. } Ax=b; x\geq0\f$.
    //! Returns the triple (v,x,y) where v is the optimal value, x is the optimal point, and y the corresponding dual feasible point.
    Tuple< AX, Vector<AX>, Vector<AX> >
    hotstarted_minimise(const Vector<X>& c, const Vector<X>& xl, const Vector<X>& xu, const Matrix<X>& A, const Vector<X>& b,
                        Vector<AX>& x, Vector<AX>& y, Vector<AX>& zl, Vector<AX>& zu) const;

    //! \brief Test feasibility of the problem \f$Ax=b; x_l\leq x\leq x_u\f$.
    //! Returns the pair (r,x) where r is the result, and x the (potential) feasible point.
    ValidatedKleenean
    feasible(const Vector<X>& xl, const Vector<X>& xu, const Matrix<X>& A, const Vector<X>& b) const;


    //! \brief Validate that \a x is primal feasible and \a y is dual feasible.
    //! Returns the interval of possible optimal values.
    ValidatedKleenean validate_feasibility(const Vector<X>& xl, const Vector<X>& xu,
                                           const Matrix<X>& A, const Vector<X>& b,
                                           const Vector<AX>& x, const Vector<AX>& y) const;
  public:
    //! \brief Perform a step of the optimization of \f$\min c^T x \text{ s.t. } Ax=b; x_l \leq x\leq x_u\f$.
    //! Returns true if a full Newton step (alpha=1) is taken. In this case, the problem is feasible (up to roundoff error).
    LinearProgramStatus
    _minimisation_step(const Vector<X>& c,
                       const Vector<X>& xl, const Vector<X>& xu,
                       const Matrix<X>& A, const Vector<X>& b,
                       Vector<AX>& x, Vector<AX>& y, Vector<AX>& zl, Vector<AX>& zu) const;

    //! \brief Perform a step of the feasibility problem \f$Ax=b,\ x_l \leq x \leq x_u\f$.
    LinearProgramStatus
    _feasibility_step(const Vector<X>& xl, const Vector<X>& xu,
                      const Matrix<X>& A, const Vector<X>& b,
                      Vector<AX>& x, Vector<AX>& y, Vector<AX>& zl, Vector<AX>& zu) const;


};

} // namespace Ariadne

#endif
