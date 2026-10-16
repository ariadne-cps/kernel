/***************************************************************************
 *            solving/nonlinear_programming.cpp
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
#include "solving/nonlinear_programming.hpp"

namespace Ariadne {

using namespace Detail;

const FloatDP NonlinearProgrammingSolverBase::zero = FloatDP(0,dp);
const FloatDP NonlinearProgrammingSolverBase::one = FloatDP(1,dp);

Bool NonlinearProgrammingSolverBase::
almost_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ApproximateVectorType ax, FloatDPApproximation error) const
{
    ExactVectorType ex=exact_vector(ax);
    if(!contains(D,ex)) { return false; }
    ApproximateVectorType gx=g(ax);
    return probably(contains(cast_exact_widen(C,cast_exact(error)),gx));
}


Bool NonlinearProgrammingSolverBase::
is_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ExactVectorType x) const
{
    if(!contains(D,x)) { return false; }
    Vector<FloatDPBounds> gx=g(x);
    return definitely(contains(C,gx));
}


ValidatedKleenean NonlinearProgrammingSolverBase::
contains_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ValidatedVectorType X) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", X="<<X);

    // Now test if the (reduced) box X satisfies other constraints
    if(definitely(disjoint(Box<UpperIntervalType>(X),D))) { return false; }
    if(definitely(not subset(Box<UpperIntervalType>(X),D))) { return indeterminate; }

    // Test inequality constraints
    ValidatedKleenean result = true;
    Vector<FloatDPBounds> gx=g(X);
    LOGGING_PRINTLN("g(X)="<<gx);
    for(SizeType i=0; i!=C.size(); ++i) {
        if(definitely(disjoint(UpperIntervalType(gx[i]),C[i]))) {
            return false;
        }
        if(!C[i].is_singleton()) {
            if(definitely(not subset(UpperIntervalType(gx[i]),C[i]))) { result = indeterminate; }
        }
    }

    // Break if some inequality constraints indefinite
    if(!definitely(result)) { return result; }

    // Extract the equality constraints
    List<SizeType> equality_constraints;
    equality_constraints.reserve(C.size());
    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].is_singleton()) { equality_constraints.append(i); }
    }

    // Construct the function g_e(x) = g_{e_i}(x)
    ARIADNE_ASSERT(g.result_size()>0);
    ValidatedVectorMultivariateFunction ge(equality_constraints.size(),g.domain());
    ExactBoxType ce(equality_constraints.size());
    for(SizeType i=0; i!=ge.result_size(); ++i) {
        ge[i]=g[equality_constraints[i]];
        ce[i]=C[equality_constraints[i]];
    }

    LOGGING_PRINTLN("ge="<<ge<<", ce="<<ce);

    // FIXME: Carefully change this code!
    FloatDPBoundsMatrix ivlA=jacobian(ge,X);
    LOGGING_PRINTLN("ivlA="<<ivlA);
    FloatDPApproximationVector fltD(X.size(),dp);
    for(SizeType i=0; i!=X.size(); ++i) { fltD[i]=rec(cast_positive(+sqr(X[i].error()))); }
    FloatDPApproximationMatrix fltA=midpoint(ivlA);
    LOGGING_PRINTLN("A="<<fltA);
    LOGGING_PRINTLN("D="<<fltD);
    FloatDPApproximationMatrix fltL = FloatDPApproximationDiagonalMatrix(fltD.array())*transpose(fltA);
    LOGGING_PRINTLN("L="<<fltL);

    FloatDPBoundsMatrix ivlS = ivlA * cast_exact(fltL);
    LOGGING_PRINTLN("ivlS="<<ivlS);

    FloatDPBoundsMatrix ivlR = inverse(ivlS);
    try {
        ivlR=inverse(ivlS);
    }
    catch (SingularMatrixException const&) {
        return indeterminate;
    }

    LOGGING_PRINTLN("ivlR="<<ivlR);
    FloatDPBoundsMatrix& valR=reinterpret_cast<FloatDPBoundsMatrix&>(ivlR);

    // Projected interval Newton step. For h:R^n->R^m; Dh mxn, take L nxm.
    // ExactIntervalType Newton update X' = x - L * (Dh(X)*L)^{-1} * h(x)
    // Choose L = rad(X)^2 Dh(x)^T where rad(X) is the diagonal matrix of radii of X
    Vector<FloatDPBounds> x=midpoint(X);
    Vector<FloatDPBounds> new_X = x - cast_exact(fltL) * (valR * (ge(x)-cast_singleton(ce)) );
    LOGGING_PRINTLN("old_X="<<X);
    LOGGING_PRINTLN("new_X="<<new_X);
    Vector<FloatDPBounds> reduced_X = refinement(X,new_X);
    LOGGING_PRINTLN("reduced_X="<<reduced_X);

    if(refines(new_X,X)) { return true; }
    else { return indeterminate; }
}




Bool NonlinearProgrammingSolverBase::
validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     ExactVectorType x0, ExactVectorType) const
{
    return this->validate_feasibility(D,g,C,x0);
}

Bool NonlinearProgrammingSolverBase::
validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     ExactVectorType x0) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x0.size()==D.size());
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    LOGGING_PRINTLN("x0="<<x0);

    Vector<FloatDPBounds> x(x0);
    LOGGING_PRINTLN("x="<<x);

    Vector<FloatDPBounds> gx=g(x);
    LOGGING_PRINTLN("gx="<<gx);

    List<SizeType> equalities, inequalities;
    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].lower_bound()==C[i].upper_bound()) {
            equalities.append(i);
        } else {
            inequalities.append(i);
            if(!definitely(contains(C[i],gx[i]))) {
                LOGGING_PRINTLN_AT(1,"g["<<i<<"](x)="<<gx[i]<<", C["<<i<<"]="<<C[i]);
                return false; }
        }
    }

    if(equalities.empty()) { LOGGING_PRINTLN("Feasible"); return true; }

    SizeType k=equalities.size();
    SizeType n=D.size();
    ValidatedVectorMultivariateFunction h(equalities.size(),g.domain());
    ExactFloatDPVectorType c(equalities.size(),dp);
    for(SizeType i=0; i!=equalities.size(); ++i) {
        h[i] = g[equalities[i]];
        c[i] = C[equalities[i]].lower_bound();
    }
    LOGGING_PRINTLN("h="<<h<<" c="<<c<<" h(x)-c="<<(h(x0)-c));

    // Attempt to solve h(x0+AT*w)=0
    // TODO: Change to use validated numbers
    Matrix<FloatDPBounds> AT = transpose(h.jacobian(x0));
    LOGGING_PRINTLN("A="<<transpose(AT));
    Vector<FloatDPBounds> w0(k,FloatDPBounds(0,dp));

    Bool found_solution=false;
    Bool validated_solution=false;

    Vector<FloatDPBounds> w(k,dp), mw(k,dp), nw(k,dp);
    Vector<FloatDPBounds> mx(n,dp);

    for(SizeType ii=0; ii!=12; ++ii) {
        mw=midpoint(w);
        x=x0+AT*w;
        mx=x0+AT*mw;
        nw = mw - solve(h.jacobian(x)*AT,Vector<FloatDPBounds>(h(mx)-c));
        LOGGING_PRINTLN_AT(1,"w="<<w<<", h(x0+AT*w)="<<h(x)<<", nw="<<nw<<", refines="<<refines(nw,w));

        if(!found_solution) {
            if(refines(nw,w)) {
                found_solution=true;
                w=w+FloatDPBounds(0,1,dp)*(w0-w);
            } else {
                w=w+FloatDPBounds(0,1,dp)*(FloatDPBounds(2,dp)*nw-w);
            }
        } else {
            if(refines(nw,w)) {
                validated_solution=true;
            } else if(validated_solution) {
                // Not a contraction, so presumably accurate enough
                break;
            }
            w=refinement(nw,w);
        }

    }
    LOGGING_PRINTLN("w="<<w<<", validated="<<validated_solution);

    if(!validated_solution) { return false; }

    // Compute x value
    x=x0+AT*w;
    LOGGING_PRINTLN("x="<<x);
    gx=g(x);
    LOGGING_PRINTLN("g(x)="<<gx);

    // Check that equality constraints are plausible
    ARIADNE_DEBUG_ASSERT(models(h(x)-c,ExactFloatDPVectorType(k,dp)));

    // Check inequality constraints once more
    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].lower_bound()==C[i].upper_bound()) {
            ARIADNE_DEBUG_ASSERT(models(gx[i],C[i].midpoint()));
        } else {
            if(!definitely(element(gx[i],C[i]))) {
                return false;
            }
        }
    }
    return true;
}


Bool NonlinearProgrammingSolverBase::
validate_infeasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                       ExactVectorType x, ExactVectorType y) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());
    LOGGING_SCOPE_CREATE;
    // Compute first-order approximation to g(D) centred at x.
    // For feasibilty, have yg(D) cap yC nonempty.
    // Estimate y g(X) = y g(x) + y Dg(X).(X-x)

    // Compute y.C
    UpperIntervalType yC = dot(y,cast_vector(C));

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(SizeType j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j]; }
    UpperIntervalType tygD = apply(tyg,D);

    UpperIntervalMatrixType dgD = jacobian_range(g,cast_vector(D));
    UpperIntervalVectorType ydgD = transpose(dgD) * y;

    UpperIntervalType ygx = dot(y,apply(g,UpperIntervalVectorType(x)));

    UpperIntervalType ygD = ygx;
    for(SizeType i=0; i!=x.size(); ++i) {
        ygD += ydgD[i] * (D[i]-x[i]);
    }

    LOGGING_PRINTLN("yC="<<yC<<" tygD="<<tygD<<" ygD="<<ygD);

    if(definitely(intersection(yC,ygD).is_empty())) { LOGGING_PRINTLN("Infeasible"); return true; }
    else { return false; }
}

// FIXME: Look at this code again, especially relating to generalised Lagrange multipliers
Bool NonlinearProgrammingSolverBase::
is_infeasibility_certificate(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ExactFloatDPVectorType y) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", y="<<y);

    if(y.size()==0) { return D.is_empty(); }

    // Try to prove lambda.(g(y)-c) != 0
    const SizeType n=C.size();

    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(SizeType i=0; i!=n; ++i) {
        tyg+=y[i]*ValidatedScalarMultivariateTaylorFunctionModelDP(D,g[i],default_sweeper());
    }
    ValidatedNumericType iygx = tyg(cast_singleton(D));

    UpperIntervalType iyC(0,0);
    for(SizeType i=0; i!=n; ++i) {
        iyC+=y[i]*C[i];
    }

    if(definitely(disjoint(iyC,UpperIntervalType(iygx)))) {
        return true;
    } else {
        return false;
    }
}





ValidatedVectorType NonlinearProgrammingSolverBase::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const
{
    LOGGING_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->minimise(f,D,gh,C);
}



ValidatedKleenean NonlinearProgrammingSolverBase::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const
{
    LOGGING_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->feasible(D,gh,C);
}



} // namespace Ariadne
