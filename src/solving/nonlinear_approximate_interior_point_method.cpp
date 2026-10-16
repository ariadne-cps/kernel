/***************************************************************************
 *            solving/nonlinear_approximate_interior_point_method.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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

// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
#include "solving/nonlinear_interior_point_method.hpp"
#include "solving/nonlinear_programming.hpp"

namespace Ariadne {

using namespace nonlinear_programming_detail;

ValidatedKleenean ApproximateNonlinearInteriorPointSolver::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("D="<<D<<", h="<<h);
    FloatDPApproximationVector x=midpoint(D);
    FloatDPApproximationVector y(h.result_size(),zero);

    for(SizeType i=0; i!=8; ++i) {
        this->feasibility_step(D,h,x,y);
    }

    if( decide(norm(h(x))<1e-10) ) { return true; }

    if(!possibly(contains(UpperIntervalType(dot(UpperIntervalVectorType(exact_vector(y)),apply(h,D))),zero))) { return false; }

    return indeterminate;
}

Void ApproximateNonlinearInteriorPointSolver::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                 FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("x="<<x<<" y="<<y);
    static const double SCALE_FACTOR = 0.75;
    const SizeType n=x.size();
    const SizeType m=y.size();
    // Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));
    FloatDPApproximationMatrix A = ddhx.jacobian();
    LOGGING_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPApproximationMatrix H(n,n,dp);
    for(SizeType i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(SizeType j=0; j!=n; ++j) {
        H[j][j] += rec(sqr(x[j]-D[j].lower_bound()));
        H[j][j] += rec(sqr(D[j].upper_bound()-x[j]));
    }

    FloatDPApproximationVector rx = transpose(A) * y;
    for(SizeType j=0; j!=n; ++j) {
        rx[j] -= rec(x[j]-D[j].lower_bound());
        rx[j] += rec(D[j].upper_bound()-x[j]);
    }
    FloatDPApproximationVector ry = ddhx.value();
    LOGGING_PRINTLN("rx="<<rx<<" ry="<<ry);

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPApproximationMatrix Hinv=inverse(H);
    LOGGING_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPApproximationMatrix S=A*Hinv*transpose(A);
    FloatDPApproximationMatrix Sinv=inverse(S);
    LOGGING_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    FloatDPApproximationVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatDPApproximationVector dx = Hinv * ( rx - transpose(A) * dy);
    LOGGING_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPApproximation ax = one;
    FloatDPApproximationVector nx = x-ax*dx;
    while(!contains(D,exact_vector(nx))) {
        ax*=SCALE_FACTOR;
        nx = x - ax * dx;
    }
    FloatDPApproximationVector ny = y-ax*dy;
    LOGGING_PRINTLN("nx="<<nx<<" ax="<<ax<<" ny="<<ny);
    LOGGING_PRINTLN_AT(1,"h(x)="<<h(nx));

    x=nx; y=ny;
}



} // namespace Ariadne
