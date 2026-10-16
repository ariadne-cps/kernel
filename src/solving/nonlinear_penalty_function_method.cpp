/***************************************************************************
 *            solving/nonlinear_penalty_function_method.cpp
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
#include "solving/nonlinear_penalty_function_method.hpp"
#include "solving/nonlinear_programming.hpp"

namespace Ariadne {

using namespace nonlinear_programming_detail;

//------- NonlinearPenaltyFunctionSolver ------------------------------------------//

NonlinearPenaltyFunctionSolver* NonlinearPenaltyFunctionSolver::
clone() const
{
    return new NonlinearPenaltyFunctionSolver(*this);
}

ValidatedVectorType NonlinearPenaltyFunctionSolver::
minimise(ValidatedScalarMultivariateFunction, ExactBoxType, ValidatedVectorMultivariateFunction, ExactBoxType) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedKleenean NonlinearPenaltyFunctionSolver::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("D="<<D<<" g="<<g<<" C="<<C<<" ");

    FloatDPApproximationVector x=midpoint(D);

    FloatDPApproximationVector w=midpoint(C);
    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].upper_bound()==+infty) { w[i]=C[i].lower_bound()+one; }
        else if(C[i].lower_bound()==-infty) { w[i]=C[i].upper_bound()-one; }
    }

    FloatDPApproximationVector y(C.size(),zero);

    LOGGING_PRINTLN("x="<<x<<" w="<<w<<" y="<<y);

    for(SizeType i=0; i!=10; ++i) {
        this->feasibility_step(D,g,C,x,y,w);
    }
    return this->check_feasibility(D,g,C,exact_vector(x),exact_vector(y));
}

Void NonlinearPenaltyFunctionSolver::
feasibility_step(const ExactBoxType& X, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& W,
                 FloatDPApproximationVector& x, FloatDPApproximationVector& w, FloatDPApproximation& mu) const
{
    LOGGING_SCOPE_CREATE;
    ApproximateVectorMultivariateFunction h(0u,X.dimension());
    const SizeType n=X.size();
    const SizeType m=W.size();
    const SizeType l=h.result_size();

    LOGGING_PRINTLN("x="<<x);
    LOGGING_PRINTLN("w="<<w);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));

    mu *= 0.5;
    LOGGING_PRINTLN("mu="<<mu);

    // G is the constraint value vector
    FloatDPApproximationVector gx = ddgx.value();
    FloatDPApproximationVector hx = ddhx.value();
    LOGGING_PRINTLN("g(x)="<<gx);
    LOGGING_PRINTLN("h(x)="<<hx);

    // A is the transpose derivative matrix aij=dgi/dxj
    FloatDPApproximationMatrix A = transpose(ddgx.jacobian());
    LOGGING_PRINTLN("A=Dg(x)="<<A);
    FloatDPApproximationMatrix B = transpose(ddhx.jacobian());
    // FIXME: Due to problems with zero-element differential, need to resize matrix if no h
    if(l==0) { B.resize(n,0); }
    LOGGING_PRINTLN("B=Dh(x)="<<B);

    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatDPApproximationMatrix H(n,n,dp);
    for(SizeType j=0; j!=m; ++j) { H += (gx[j]-w[j]) * ddgx[j].hessian(); }
    for(SizeType k=0; k!=l; ++k) { H += (hx[k]) * ddhx[k].hessian(); }
    LOGGING_PRINTLN("H="<<H);

    FloatDPApproximationDiagonalMatrix D(n,dp);
    FloatDPApproximationDiagonalMatrix E(m,dp);
    for(SizeType i=0; i!=n; ++i) { D[i] = rec(sqr(x[i]-X[i].lower_bound())) + rec(sqr(X[i].upper_bound()-x[i])); }
    for(SizeType j=0; j!=m; ++j) { E[j] = rec(sqr(w[j]-W[j].lower_bound())) + rec(sqr(W[j].upper_bound()-w[j])); }
    LOGGING_PRINTLN("D="<<D);
    LOGGING_PRINTLN("E="<<E);

    FloatDPApproximationMatrix S = H + B * transpose(B);
    S += D;
    LOGGING_PRINTLN("S="<<S);

    FloatDPApproximationMatrix R=inverse(S);
    LOGGING_PRINTLN("inverse(S)="<<R);

    // Compute residuals
    FloatDPApproximationVector rx = A*gx + B * hx ; // + 1/(x.upper_bound()-x) + 1/x.lower_bound()-x if no regularisation
    FloatDPApproximationVector rw = w-gx;

    LOGGING_PRINTLN("rx="<<rx);
    LOGGING_PRINTLN("rw="<<rw);

    FloatDPApproximationVector dx = R * (rx + A * rw);
    FloatDPApproximationVector dw = rw + transpose(A)*dx;
    LOGGING_PRINTLN("dx="<<dx);
    LOGGING_PRINTLN("dw="<<dw);


    FloatDPApproximationVector newx(n,dp);
    FloatDPApproximationVector neww(m,dp);

    static const FloatDPApproximation ALPHA_SCALE_FACTOR = 0.75_approx;

    FloatDPApproximation alpha = 1.0_approx;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        alpha *= ALPHA_SCALE_FACTOR;
    } while ( !contains(X,exact_vector(newx)) || !contains(W,exact_vector(neww)) );
    alpha /= ALPHA_SCALE_FACTOR;

    LOGGING_PRINTLN("alpha="<<alpha);

    LOGGING_PRINTLN("newx="<<newx);
    LOGGING_PRINTLN("neww="<<neww<<"\n");

    x=newx;
    w=neww;

    return;
}


Void NonlinearPenaltyFunctionSolver::
feasibility_step(const ExactBoxType&, const ValidatedVectorMultivariateFunction&, const ExactBoxType&,
                 FloatDPBoundsVector&, FloatDPBoundsVector&) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


// Use a penalty approach without multipliers on the constraint functions
// Solve g(x)=w, x in D, w in C; Lagrangian y.(g(x)-w)
Void NonlinearPenaltyFunctionSolver::
feasibility_step(ExactBoxType const& D, ApproximateVectorMultivariateFunction const& g, ExactBoxType const& C,
                 FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& w) const
{
    LOGGING_SCOPE_CREATE;
    auto m=y.size(); auto n=x.size();

    FloatDPApproximationVector cl=lower_bounds(C);
    FloatDPApproximationVector cu=upper_bounds(C);
    FloatDPApproximationVector dl=lower_bounds(D);
    FloatDPApproximationVector du=upper_bounds(D);

    LOGGING_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    LOGGING_PRINTLN("dl="<<dl<<", du="<<du);
    LOGGING_PRINTLN("cl="<<cl<<", cu="<<cu);
    LOGGING_PRINTLN("w="<<w<<", x="<<x<<", y="<<y);

    ARIADNE_ASSERT_MSG(g.argument_size()==D.size(),"D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_ASSERT_MSG(g.result_size()==C.size(),  "D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    // Solve the problem
    //   minimise Sum -log(w-cl)-log(cu-w)-log(x-dl)-log(du-x)
    //   subject to g(x)-w=0

    // Lagrangian
    //   -log(w-cl)-log(cu-w)-log(x-dl)-log(du-x) - y.(g(x)-w)

    // Conditions for a constrained minimum
    // 1/(cu-w)-1/(w-cu) + y       = 0
    // 1/(du-x)-1/(x-dl) - y.Dg(x) = 0
    //          w - g(x)           = 0

    // Second-order conditions
    // (1/(w-cl)^2 + 1/(cu-w)^2) dw                 +   dy = - ( 1/(cu-w) - 1/(w-cl) + y       )
    //   (1/(x-dl)^2 + 1/(du-x)^2 - y.D^2x) dx - Dg'(x) dy = - ( 1/(du-x) - 1/(x-dl) - y.Dg(x) )
    //                           dw - Dg(x) dx             = - ( w - g(x) )

    Vector<Differential<ApproximateNumericType>> ddgx=g.evaluate(Differential<ApproximateNumericType>::variables(2,x));
    LOGGING_PRINTLN("ddgx="<<ddgx);

    Vector<ApproximateNumericType> gx = ddgx.value();
    LOGGING_PRINTLN("g(x)="<<gx);
    Matrix<ApproximateNumericType> A = ddgx.jacobian();
    LOGGING_PRINTLN("Dg(x)="<<A);

    Vector<ApproximateNumericType> yA=transpose(A)*y;

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<ApproximateNumericType> YH(x.size(),x.size(),dp);
    for(SizeType i=0; i!=y.size(); ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    LOGGING_PRINTLN("Y.D2g(x)="<<YH);

    Vector<ApproximateNumericType> recwu=cu-w; recwu=erec(recwu);
    Vector<ApproximateNumericType> recwl=w-cl; recwl=erec(recwl);
    Vector<ApproximateNumericType> recxu=du-x; recxu=erec(recxu);
    Vector<ApproximateNumericType> recxl=x-dl; recxl=erec(recxl);

    Vector<ApproximateNumericType> diagDw=esqr(recwu)+esqr(recwl);
    Matrix<ApproximateNumericType> Dw(m,m,dp); for(SizeType i=0; i!=m; ++i) { Dw[i][i]=diagDw[i]; }
    DiagonalMatrix<ApproximateNumericType> Dx(esqr(recxu)+esqr(recxl));


    for(SizeType i=0; i!=n; ++i) { YH[i][i]-=Dx[i]; }

    Matrix<ApproximateNumericType> AT=transpose(A);
    Matrix<ApproximateNumericType> Znm(n,m,dp);
    Matrix<ApproximateNumericType> Zmn(m,n,dp);
    Matrix<ApproximateNumericType> Zmm(m,m,dp);
    Matrix<ApproximateNumericType> Im=Matrix<ApproximateNumericType>::identity(m,dp);


    Matrix<ApproximateNumericType> S=cojoin(join(Dw,Zmn,Im),join(Znm,-YH,-AT),join(Im,-A,Zmm));
    Vector<ApproximateNumericType> r=join(recwu-recwl+y,recxu-recxl-yA,w-gx);

    for(SizeType j=0; j!=m; ++j) {
        if(C[j].lower_bound()==C[j].upper_bound()) {
            S[j][j]=1;
            S[j][m+n+j]=0;
            S[m+n+j][j]=0;
            r[j]=0;
        }
    }

    Vector<ApproximateNumericType> swxy = -solve(S,r);

    Vector<ApproximateNumericType> sw(m,dp),sx(n,dp),sy(m,dp);
    sw = project(swxy,range(0,m));
    sx = project(swxy,range(m,m+n));
    sy = project(swxy,range(m+n,m+n+m));

    ApproximateNumericType al=one;
    ApproximateVectorType nw=w+al*sw;
    ApproximateVectorType nx=x+al*sx;
    ApproximateVectorType ny(m,dp);
    LOGGING_PRINTLN("sx="<<sx);
    LOGGING_PRINTLN("sw="<<sw);
    while(!contains(C,exact_vector(nw)) || !contains(D,exact_vector(nx))) {
        al*=0.75;
        nw=w+al*sw;
        nx=x+al*sx;
    }
    LOGGING_PRINTLN("al="<<sw);
    ny=y+al*sy;

    w=nw; x=nx; y=ny;
}



ValidatedKleenean NonlinearPenaltyFunctionSolver::
check_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     ExactFloatDPVectorType x, ExactFloatDPVectorType y) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());

    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1,"D="<<D<<" C="<<C);

    FloatDPBoundsVector gx=g(x);
    LOGGING_PRINTLN_AT(1,"x="<<x<<" y="<<y<<" g(x)="<<gx);

    ValidatedKleenean result = true;

    List<SizeType> equalities;
    for(SizeType j=0; j!=C.size(); ++j) {
        if( definitely(gx[j].upper()<C[j].lower_bound() || gx[j].lower()>C[j].upper_bound()) ) {
            return false;
        }
        if(C[j].lower_bound()==C[j].upper_bound()) {
            equalities.append(j);
        } else {
            if(!definitely(contains(C[j],gx[j]))) { result = indeterminate; }
        }
    }

    if(definitely(result)) {
        if(equalities.empty()) { LOGGING_PRINTLN("Feasible"); return true; }

        ValidatedVectorMultivariateFunction h(equalities.size(),g.domain());
        FloatDPBoundsVector c(equalities.size(),dp);
        for(SizeType i=0; i!=equalities.size(); ++i) {
            h[i] = g[equalities[i]];
            c[i] = C[equalities[i]].lower_bound();
        }
        LOGGING_PRINTLN_AT(1,"g="<<g);
        LOGGING_PRINTLN_AT(1,"h="<<h<<" c="<<c<<" h(x)-c="<<(h(x)-c));

        static const ExactDouble SEARCH_TOLERANCE=1e-8_pr;
        FloatDPBoundsVector W(h.result_size(),FloatDPBounds(-1,+1,dp)*SEARCH_TOLERANCE);
        FloatDPBoundsMatrix AT = transpose(midpoint(h.jacobian(x)));
        FloatDPBoundsVector B = x+AT*W;
        FloatDPBoundsMatrix IA = h.jacobian(B);
        LOGGING_PRINTLN_AT(1,"AT="<<AT<<" IA="<<IA);
        LOGGING_PRINTLN_AT(1,"B="<<B);

        // Perform an interval Newton step to try to attain feasibility
        FloatDPBoundsVector nW = inverse(IA*AT) * FloatDPBoundsVector(h(x)-cast_exact(c));
        LOGGING_PRINTLN_AT(1,"W="<<W<<"\nnew_W="<<nW);
        if(definitely(subset(UpperBoxType(B),D)) && refines(nW,W)) { LOGGING_PRINTLN("Feasible"); return true; }
        else { result=indeterminate; }
    }

    // Compute first-order approximation to g(D) centred at x.
    // For feasibilty, have yg(D) cap yC nonempty.
    // Estimate y g(X) = y g(x) + y Dg(X).(X-x)

    // Compute y.C
    UpperIntervalVectorType iy(y);
    UpperIntervalType yC = dot(iy,C);

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(SizeType j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j]; }
    UpperIntervalType tygD = UpperIntervalType(tyg(cast_singleton(D)));

    UpperIntervalMatrixType dgD = jacobian_range(g,cast_vector(D));
    UpperIntervalVectorType ydgD = transpose(dgD) * y;

    FloatDPBounds ygx = dot(y,gx);

    UpperIntervalType ygD = UpperIntervalType(ygx);
    for(SizeType i=0; i!=x.size(); ++i) {
        ygD += ydgD[i] * (D[i]-UpperIntervalType(x[i]));
    }

    LOGGING_PRINTLN_AT(1,"yC="<<yC<<" tygD="<<tygD<<" ygD="<<ygD);

    if(definitely(is_empty(intersection(yC,ygD)))) { LOGGING_PRINTLN("Infeasible"); return false; }
    else { return indeterminate; }
}






/*

// Solve max log(x-xl) + log(xu-x) + log(zu-z) + log(z-zl) such that g(x)=z
//   if zl[i]=zu[i] then z=zc is hard constraint
//   alternatively, use the penalty (z[j]-zc[j])^2/2 instead

// KKT conditions
//     1/(x-xl) - 1/(xu-x) + y Dg(x) = 0
//     1/(z-zl) - 1/(zu-z) - y = 0
//     g(x) - z = 0
//   If zu[j]=inf, then 1/(z[j-zl[j]) - y[j] = 0, so y[j]>=0
//   If zl[j]=zu[j], then use the equation z[j]=zc[j] instead
//
// Re-write KKT conditions for x,z as
//     (xu-xl) + y g(x) (x-xl)(xu-x) = 0
//     (zu-zl) - y(z-zl)(zu-z) = 0
//  Or 1 - y(z-zl)(zu-z)/(zu-zl) = 0
//   If zu[j] = inf, then 1 - y(z-zl) = 0
//      zl[j]=zu[j]=zc[j], then y(z-zc)^2 = 0
//
// Derivative matrix
//     - 1/(x-xl)^2 - 1/(xu-x)^2 + y D^2g = 0
//
// PROBLEM:
//   Dg can be singular, even at intermediate points
//
// FJ conditions
//     mu/(x-xl) - mu/(xu-x) + y D g(x) = 0
//     mu/(z-zl) - mu/(zu-z) - y = 0
//     g(x) - z = 0
//     sum y^2 - mu = 0
//   If zu[j]=inf, then mu/(z[j-zl[j]) - y[j] = 0, so y[j]>=0
//   If zl[j]=zu[j], then -(z-zc)/mu - y = 0 instead
//
// Re-write FJ conditions for x,z as
// Derivative matrix
//     - mu/(x-xl)^2 dx - mu/(xu-x)^2 dx + y D^2g dx + Dg^T dy + (1/(x-xl) - (1/xu-x)) dmu
//     - mu/(z-zl)^2 dz - mu/(zu-z)^2 dz + I dy + (1/(z-zl) - (1/zu-z)) dmu
//     Dg dx - I dz
//    2y . dy - dmu
Void NonlinearPenaltyFunctionSolver::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                 RawFloatDPVector& x, RawFloatDPVector& y, RawFloatDPVector& z) const
{
    LOGGING_PRINTLN_AT(2,"feasibility_step");
    RawFloatDPVector xl=lower_bounds(D); RawFloatDPVector xu=upper_bounds(D);
    RawFloatDPVector zl=lower_bounds(C); RawFloatDPVector zu=upper_bounds(C);

    const SizeType n=x.size();
    const SizeType m=y.size();

    LOGGING_PRINTLN_AT(4,"x="<<x<<" y="<<y<<" z="<<z);
    Vector<FloatDPDifferential> ddx = FloatDPDifferential::variables(2,x);
    Vector<FloatDPDifferential> ddgx = g.evaluate(ddx);

    FloatDPMatrix A=ddgx.jacobian();
    LOGGING_PRINTLN_AT(6,"A="<<A);
    RawFloatDPVector v = join(join(x,z),y);

    RawFloatDPVector r(n+2*m,n+2*m);
    project(r,range(0,n)) = y * A;
    for(SizeType i=0; i!=n; ++i) {
        r[i] += ( rec(x[i]-xl[i]) - rec(xu[i]-x[i]) );
    }
    for(SizeType j=0; j!=m; ++j) {
        if(zl[j]==zu[j]) { assert(zu[j]==zl[j]); r[n+j] = z[j]-zl[j]; }
        else { r[n+j] = ( rec(z[j]-zl[j]) - rec(zu[j]-z[j]) - y[j] ); }
    }
    project(r,range(n+m,n+2*m)) = ddgx.value() - z;
    r[n+2*m]=0.0;
    LOGGING_PRINTLN_AT(5,"r="<<r);

    FloatDPMatrix S(n+2*m+1,n+2*m+1);
    for(SizeType j=0; j!=m; ++j) {
        FloatDPMatrix H=ddgx[j].hessian();
        for(SizeType i1=0; i1!=n; ++i1) {
            for(SizeType i2=0; i2!=n; ++i2) {
                S[i1][i2]+=y[j]*H[i1][i2];
            }
        }
    }
    for(SizeType j=0; j!=m; ++j) {
        for(SizeType i=0; i!=n; ++i) {
            S[i][j+m+n]=A[j][i];
            S[j+m+n][i]=A[j][i];
        }
    }
    for(SizeType j=0; j!=m; ++j) {
        S[n+j][n+m+j] = -1.0;
        S[n+m+j][n+j] = -1.0;
        //if(zl[j]==zu[j]) { S[n+j][n+j] = -1.0; S[n+j][n+m+j] = 0.0; }
        if(zl[j]==zu[j]) { S[n+j][n+j] = +inf; }
        else { S[n+j][n+j] = - rec(sqr(z[j]-zl[j])) - rec(sqr(zu[j]-z[j])); }
    }
    for(SizeType i=0; i!=n; ++i) {
        S[i][i]-= rec(sqr(xu[i]-x[i]));
        S[i][i]-= rec(sqr(x[i]-xl[i]));
    }

    for(SizeType i=0; i!=n; ++i) {
        S[i][n+2*m] -= rec(xu[i]-x[i]);
        S[i][n+2*m] += rec(x[i]-xl[i]);
    }

    for(SizeType j=0; j!=n; ++j) {
        //S[n+m+j][n+m+j] = -1.0/1024;
    }

    LOGGING_PRINTLN_AT(5,"S="<<S);
    //LOGGING_PRINTLN_AT(5,"S="<<std::fixed<<pretty(S));

    FloatDPMatrix Sinv = inverse(S);
    //LOGGING_PRINTLN_AT(9,"Sinv="<<Sinv<<"\n);
    //LOGGING_PRINTLN_AT(5,"Sinv="<<std::fixed<<pretty(Sinv));

    RawFloatDPVector dv = Sinv * r;
    LOGGING_PRINTLN_AT(5,"dv="<<dv);

    FloatDP alpha = 1.0;
    RawFloatDPVector nv = v-dv;
    while(!contains(D,RawFloatDPVector(project(nv,range(0,n)))) || !contains(C,RawFloatDPVector(project(nv,range(n,n+m)))) ) {
        alpha *= 0.75;
        nv = v-alpha*dv;
    }

    LOGGING_PRINTLN_AT(4,"nv="<<nv<<" a="<<alpha);

    x=project(nv,range(0,n));
    z=project(nv,range(n,n+m));
    y=project(nv,range(n+m,n+2*m));
    LOGGING_PRINTLN_AT(4,"g(x)-z="<<g(x)-z);

}
*/


} // namespace Ariadne
