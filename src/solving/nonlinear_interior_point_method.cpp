/***************************************************************************
 *            solving/nonlinear_interior_point_method.cpp
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
#include "solving/nonlinear_programming.tpl.hpp"

namespace Ariadne {

using namespace Detail;

//------- NonlinearInteriorPointSolver -----------------------------------//

ValidatedVectorType NonlinearInteriorPointSolver::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("f="<<f<<" D="<<D<<" g="<<g<<" C="<<C);

    // Reject inconsistent public inputs before evaluating g over D or
    // constructing the primal-dual Newton system.
    ARIADNE_PRECONDITION(f.argument_size()==D.size());
    ARIADNE_PRECONDITION(g.argument_size()==D.size());
    ARIADNE_PRECONDITION(g.result_size()==C.size());

    ValidatedVectorMultivariateFunction h(0,D.dimension());

    UpperBoxType gD = apply(g,D);
    if(definitely(disjoint(gD,C))) { throw InfeasibleProblemException(); }

    FloatDPApproximationVector x = midpoint(D);
    FloatDPApproximationVector w = midpoint(intersection(UpperBoxType(gD),C));

    FloatDPApproximationVector kappa(g.result_size(),zero);
    FloatDPApproximationVector lambda(h.result_size(),zero);
    FloatDPApproximation mu = one;


    for(SizeType i=0; i!=12; ++i) {
        this->minimisation_step(f,D,g,C,h, x,w, kappa,lambda, mu);
        // For the twelve iterations i=0,...,11, i%3==0 only at 0,3,6,9.
        // An additional i<=10 guard would therefore always be true.
        if(i%3==0) { mu *= 0.25_exact; }
    }

    const ExactFloatDPVectorType candidate=exact_vector(x);
    // The returned box is a singleton. A small numerical change, or even
    // a certificate for a nearby root, would not prove that this point
    // satisfies the constraints. Validate the actual returned point.
    if(!this->is_feasible_point(D,g,C,candidate)) {
        LOGGING_PRINTLN("Final minimisation point is not certified feasible");
        throw IndeterminateFeasibilityException();
    }
    return ValidatedVectorType(candidate);
}



// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


// min f(x) | x\in D & w\in C | g(x) = w & h(x) = 0
// Lagrange multipliers kappa d(g(x)-w); lambda dh(x)
Void NonlinearInteriorPointSolver::
minimisation_step(const ApproximateScalarMultivariateFunction& f, const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c, const ApproximateVectorMultivariateFunction& h,
                  FloatDPApproximationVector& x, FloatDPApproximationVector& w,
                  FloatDPApproximationVector& kappa, FloatDPApproximationVector& lambda, const FloatDPApproximation& mu) const
{
    const SizeType n=x.size();
    const SizeType m=kappa.size();
    const SizeType l=lambda.size();

    ARIADNE_DEBUG_PRECONDITION(w.size()==kappa.size());
    ARIADNE_DEBUG_PRECONDITION(f.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(g.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(h.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(g.result_size()==m);
    ARIADNE_DEBUG_PRECONDITION(h.result_size()==l);
    ARIADNE_DEBUG_PRECONDITION(contains(d,exact_vector(x)));
    ARIADNE_DEBUG_PRECONDITION(contains(c,exact_vector(w)));
    ARIADNE_DEBUG_PRECONDITION(mu.raw()>0);

    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("x="<<x);
    LOGGING_PRINTLN("w="<<w);
    LOGGING_PRINTLN_AT(1,"kappa="<<kappa);
    LOGGING_PRINTLN_AT(1,"lambda="<<lambda);
    LOGGING_PRINTLN_AT(1,"mu="<<mu);

    FloatDPApproximationVector slack(2*n,dp);
    FloatDPApproximationVectorRange slackl(slack,range(0,n));
    FloatDPApproximationVectorRange slacku(slack,range(n,2*n));

    FloatDPApproximationDifferential ddfx=f.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));

    // G is the constraint value vector
    FloatDPApproximation fx = ddfx.value();
    FloatDPApproximationVector gx = ddgx.value();
    FloatDPApproximationVector hx = ddhx.value();
    LOGGING_PRINTLN("f(x)="<<fx);
    LOGGING_PRINTLN("g(x)="<<gx);
    LOGGING_PRINTLN("h(x)="<<hx);
    LOGGING_PRINTLN_AT(1,"g(x)-w="<<(gx-w));

    // A, B are the derivative matrices aij=dgi/dxj
    // HACK: Need to explicitly set size of Jacobian if g or h have result_size of zero
    FloatDPApproximationVector df = transpose(ddfx.gradient());
    LOGGING_PRINTLN_AT(1,"df(x)="<<df);
    FloatDPApproximationMatrix A = ddgx.jacobian();
    if(m==0) { A=FloatDPApproximationMatrix(m,n,dp); }
    LOGGING_PRINTLN("A="<<A);
    FloatDPApproximationMatrix B = ddhx.jacobian();
    if(l==0) { B=FloatDPApproximationMatrix(l,n,dp); }
    LOGGING_PRINTLN("B="<<B);



    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j]kappa[j]*dg[j]/dx[i1]dx[i2] + Sum[k]lambda[k]*dh[k]/dx[i1]dx[i2]
    FloatDPApproximationMatrix H = ddfx.hessian();
    for(SizeType j=0; j!=m; ++j) { H += kappa[j] * ddgx[j].hessian(); }
    for(SizeType k=0; k!=l; ++k) { H += lambda[k] * ddhx[k].hessian(); }
    LOGGING_PRINTLN("H="<<H);

    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    // Compute the residuals and contributions from slack in x and w
    //   rx = df/dx[i] + Sum[j] dg[j]/dx[i] * kappa[j] + Sum[k] dh[k]/dx[i] * lambda[j] + mu *( 1/(xu[i]-x[i]) - 1/(x[i]-xl[i]) )
    FloatDPApproximationVector rx = df + transpose(A) * kappa + transpose(B)* lambda;
    FloatDPApproximationDiagonalMatrix D(n,dp);
    for(SizeType i=0; i!=n; ++i) {
        FloatDPApproximation nuu = rec(d[i].upper_bound()-x[i]);
        FloatDPApproximation nul = rec(x[i]-d[i].lower_bound());
        rx[i] += mu * ( nuu - nul );
        D[i] = mu * ( nuu*nuu + nul*nul );
    }

    //   rw = - kappa[j] + mu *( 1/(wu[i]-w[i]) - 1/(w[i]-wl[i]) )
    FloatDPApproximationVector rw = -kappa;
    FloatDPApproximationDiagonalMatrix C(m,dp);
    for(SizeType j=0; j!=m; ++j) {
        FloatDPApproximation nuu = rec(c[j].upper_bound()-w[j]);
        FloatDPApproximation nul = rec(w[j]-c[j].lower_bound());
        rw[j] += (mu*EQUALITY_RELAXATION_MULTIPLIER) * ( nuu - nul );
        C[j] = (mu*EQUALITY_RELAXATION_MULTIPLIER) * ( nuu*nuu + nul*nul );
    }

    //   rkappa = g(x) - w
    FloatDPApproximationVector rkappa = gx - w;

    //   rlambda = h(x)
    FloatDPApproximationVector const& rlambda = hx;

    LOGGING_PRINTLN("rx="<<rx);
    LOGGING_PRINTLN("rw="<<rw);
    LOGGING_PRINTLN("rkappa="<<rkappa);
    LOGGING_PRINTLN("rlambda="<<rlambda);

    // Solve the equations
    //   H+D dx        + AT dk + BT dl = rx
    //            C dw -  I dk         = rw
    //    A  dx - I dw                 = rk
    //    B  dx                        = rl

    // Eliminate dw, dk without fill-in to obtain
    //   (H+D+ATCA) dx + BT dl = rx + AT rw + ATC rk
    //         B    dx         = rl

    // Set S=(H+D+ATCA); invert, and eliminate dx
    //   dx = Sinv * (rx + AT rw + ATC rk - BT dl)
    //   (B * Sinv * BT) dl = B * Sinv * (rx + AT rw + ATC rk) - rl
    FloatDPApproximationMatrix& S=H;
    S+=D;
    S+=FloatDPApproximationMatrix(transpose(A))*C*A;
    LOGGING_PRINTLN("S="<<S);

    FloatDPApproximationMatrix Sinv=inverse(S);
    LOGGING_PRINTLN("R=Sinv="<<Sinv);

    FloatDPApproximationMatrix BSinvBT = (B*Sinv)*FloatDPApproximationMatrix( transpose(B) );
    LOGGING_PRINTLN("B*inverse(S)*BT="<<BSinvBT);
    LOGGING_PRINTLN("inverse(B*inverse(S)*BT)="<<inverse(BSinvBT));

    FloatDPApproximationVector rr = Sinv * (rx + transpose(A) * (rkappa * C + rw));
    FloatDPApproximationVector dlambda = inverse(BSinvBT) * (B * rr - rlambda);
    FloatDPApproximationVector dx = rr - transpose(B*Sinv) * dlambda;
    FloatDPApproximationVector dw = A * dx - rkappa;
    FloatDPApproximationVector dkappa = rw - C * dw;

    static const FloatDPApproximation ALPHA_SCALE_FACTOR = 0.75_approx;
    static const FloatDPApproximation MINIMUM_ALPHA = 1e-16_approx;

    // Compute distance to move variables preserving feasibility
    // FIXME: Current implementation might fail due to getting too close to boundary!
    FloatDPApproximationVector newx(n,dp);
    FloatDPApproximationVector neww(m,dp);
    FloatDPApproximation alpha = 1.0_approx;
    Bool success = false;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        if (contains(d,exact_vector(newx)) && contains(c,exact_vector(neww))) { success = true; }
        else { alpha *= ALPHA_SCALE_FACTOR; }
        if(probably(alpha<MINIMUM_ALPHA)) { throw NearBoundaryOfFeasibleDomainException(); }
    } while(!success);
    LOGGING_PRINTLN("alpha="<<alpha);

    FloatDPApproximationVector newlambda = lambda - alpha * dlambda;
    FloatDPApproximationVector newkappa = kappa - alpha * dkappa;

    LOGGING_PRINTLN("newx="<<newx);
    LOGGING_PRINTLN("neww="<<neww);
    LOGGING_PRINTLN("newkappa="<<newkappa);
    LOGGING_PRINTLN("newlambda="<<newlambda);

    x=newx; w=neww; kappa=newkappa; lambda=newlambda;
}



namespace {

// Positive slack means that the corresponding facet is satisfied with
// margin t. Infinite constraint bounds are inactive, with unit slack.
FloatDPApproximationVector feasibility_slack(
    const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    const FloatDPApproximationVector& y, const FloatDPApproximation& t)
{
    const SizeType m=d.size(), n=c.size();
    FloatDPApproximationVector z(2*(m+n),dp);
    const auto gy=g(y);
    for(SizeType j=0; j!=n; ++j) {
        z[j]=is_inf(c[j].upper_bound()) ? FloatDPApproximation(1,dp) : c[j].upper_bound()-gy[j]-t;
        z[n+j]=is_inf(c[j].lower_bound()) ? FloatDPApproximation(1,dp) : gy[j]-c[j].lower_bound()-t;
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=d[i].upper_bound()-y[i]-t;
        z[2*n+m+i]=y[i]-d[i].lower_bound()-t;
    }
    return z;
}

} // namespace

ValidatedKleenean NonlinearInteriorPointSolver::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("d="<<d<<", g="<<g<<", c="<<c);
    ARIADNE_PRECONDITION(g.argument_size()==d.size());
    ARIADNE_PRECONDITION(g.result_size()==c.size());
    if(d.is_empty() || c.is_empty()) { return false; }

    FloatDPApproximation t(dp);
    FloatDPApproximationVector x,y;
    this->setup_feasibility(d,g,c,x,y,t);
    if(this->is_feasible_point(d,g,c,exact_vector(y))) { return true; }

    for(SizeType i=0; i!=12; ++i) {
        try {
            this->feasibility_step(d,g,c,x,y,t);
        } catch(const SingularMatrixException&) {
            break;
        } catch(const NearBoundaryOfFeasibleDomainException&) {
            break;
        }
        LOGGING_PRINTLN_AT(1,"t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x);
        if(this->is_feasible_point(d,g,c,exact_vector(y))) { return true; }
    }
    // x contains weights for upper/lower constraint facets followed by
    // domain facets. Only their signed constraint differences multiply g.
    FloatDPApproximationVector lambda(c.size(),dp);
    for(SizeType j=0; j!=c.size(); ++j) { lambda[j]=x[j]-x[c.size()+j]; }
    if(this->is_infeasibility_certificate(d,g,c,exact_vector(lambda))) { return false; }
    return indeterminate;
}

Void NonlinearInteriorPointSolver::feasibility_step(
    const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    ARIADNE_PRECONDITION(g.argument_size()==d.size());
    ARIADNE_PRECONDITION(g.result_size()==c.size());
    ARIADNE_PRECONDITION(x.size()==2*(d.size()+c.size()));
    ARIADNE_PRECONDITION(y.size()==d.size());
    // A hot start without a margin is made strictly interior before the
    // same primal-dual step used by the persistent-margin driver.
    FloatDPApproximation t(0,dp);
    const auto slack=feasibility_slack(d,g,c,y,t);
    for(SizeType k=0; k!=slack.size(); ++k) { t=min(t,slack[k]); }
    // An absolute unit offset alone can round away for large violations.
    t-=max(FloatDPApproximation(1,dp),abs(t));
    this->feasibility_step(d,g,c,x,y,t);
}

Void NonlinearInteriorPointSolver::feasibility_step(
    const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const
{
    LOGGING_SCOPE_CREATE;
    const SizeType m=d.size(), n=c.size(), l=2*(m+n);
    ARIADNE_PRECONDITION(g.argument_size()==m);
    ARIADNE_PRECONDITION(g.result_size()==n);
    ARIADNE_PRECONDITION(x.size()==l);
    ARIADNE_PRECONDITION(y.size()==m);

    // Maximise the common margin t, with z=s(y)-t>0, x>0 and sum(x)=1.
    // Columns of A are gradients of -s, followed by the t derivative.
    const auto ddgy=g.evaluate(FloatDPApproximationDifferential::variables(2,y));
    FloatDPApproximationMatrix A(m+1,l,dp), H(m+1,m+1,dp);
    List<SizeType> active;
    for(SizeType j=0; j!=n; ++j) {
        if(!is_inf(c[j].upper_bound())) {
            active.append(j); A[m][j]=1;
            for(SizeType i=0; i!=m; ++i) { A[i][j]=ddgy[j].gradient()[i]; }
        }
        if(!is_inf(c[j].lower_bound())) {
            active.append(n+j); A[m][n+j]=1;
            for(SizeType i=0; i!=m; ++i) { A[i][n+j]=-ddgy[j].gradient()[i]; }
        }
    }
    for(SizeType i=0; i!=m; ++i) {
        ARIADNE_PRECONDITION(!is_inf(d[i].lower_bound()));
        ARIADNE_PRECONDITION(!is_inf(d[i].upper_bound()));
        active.append(2*n+i); active.append(2*n+m+i);
        A[i][2*n+i]=1; A[i][2*n+m+i]=-1;
        A[m][2*n+i]=1; A[m][2*n+m+i]=1;
    }
    if(active.empty()) { t=0; return; }

    const auto z=feasibility_slack(d,g,c,y,t);
    FloatDPApproximationVector weights(l,dp);
    for(const auto k : active) {
        ARIADNE_PRECONDITION(decide(x[k]>0));
        ARIADNE_PRECONDITION(decide(z[k]>0));
        weights[k]=x[k];
    }
    for(SizeType j=0; j!=n; ++j) {
        const auto hj=ddgy[j].hessian();
        for(SizeType i=0; i!=m; ++i) {
            for(SizeType k=0; k!=m; ++k) { H[i][k]+=(weights[j]-weights[n+j])*hj[i][k]; }
        }
    }
    const auto mu=dot(weights,z)/active.size();
    auto rx=emul(weights,z);
    for(const auto k : active) { rx[k]-=mu*0.125_x; }
    FloatDPApproximationVector residual=-(A*weights); residual[m]+=1;
    auto S=H;
    adat(S,A,ediv(weights,z));
    const auto dyt=inverse(S)*(residual+A*ediv(rx,z));
    const FloatDPApproximationVector dz=-(transpose(A)*dyt);
    const auto dx=-ediv(FloatDPApproximationVector(rx+emul(weights,dz)),z);
    LOGGING_PRINTLN("mu="<<mu<<", dx="<<dx<<", dyt="<<dyt);

    // Re-evaluate nonlinear slacks; a positive linear prediction alone
    // cannot preserve the barrier. Bound backtracking even on failure.
    FloatDPApproximation alpha(1,dp);
    while(decide(alpha>0x1p-52_x)) {
        const FloatDPApproximationVector nx=weights+alpha*dx;
        const FloatDPApproximationVector ny=y+alpha*FloatDPApproximationVector(project(dyt,range(0,m)));
        const FloatDPApproximation nt=t+alpha*dyt[m];
        const auto nz=feasibility_slack(d,g,c,ny,nt);
        Bool acceptable=true;
        for(const auto k : active) {
            if(!decide(nx[k]>0 && nz[k]>0 && nx[k]*nz[k]>mu*0.0009765625_x)) {
                acceptable=false; break;
            }
        }
        if(acceptable) { x=nx; y=ny; t=nt; return; }
        alpha*=0.75_x;
    }
    throw NearBoundaryOfFeasibleDomainException();
}

FloatDPApproximation NonlinearInteriorPointSolver::
compute_mu(const ExactBoxType&, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
           const FloatDPApproximationVector& x, const FloatDPApproximationVector& lambda) const
{
    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions
    // With no constraints there are no multiplier contributions to
    // average. Avoid the otherwise undefined division by C.size()==0.
    if(C.size()==0u) { return FloatDPApproximation(0.0_x,dp); }

    FloatDPApproximation mu=zero;
    FloatDPApproximationVector gx = g(x);

    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].lower_bound()==C[i].upper_bound()) { }
        else if(is_inf(C[i].lower_bound())) { mu += lambda[i] * (gx[i] - C[i].upper_bound()); }
        else if(is_inf(C[i].upper_bound())) { mu += lambda[i] * (gx[i] - C[i].lower_bound()); }
        else { // std::cerr<<"FIXME: Compute mu for singleton constraint\n";
            if ( decide(lambda[i] <=0.0) ) { mu += lambda[i] * (gx[i] - C[i].upper_bound()); }
            else { mu += lambda[i] * (gx[i] - C[i].lower_bound()); }
        }
    }
    mu /= C.size();
    return mu;
}


Void NonlinearInteriorPointSolver::
setup_feasibility(const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
                  FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    ARIADNE_PRECONDITION(g.argument_size()==d.size());
    ARIADNE_PRECONDITION(g.result_size()==c.size());
    ARIADNE_PRECONDITION(!d.is_empty());
    for(SizeType i=0; i!=d.size(); ++i) {
        ARIADNE_PRECONDITION(!is_inf(d[i].lower_bound()));
        ARIADNE_PRECONDITION(!is_inf(d[i].upper_bound()));
    }
    const SizeType m=d.size(), n=c.size(), l=2*(m+n);
    y=midpoint(d);
    x=FloatDPApproximationVector(l,dp);
    SizeType count=2*m;
    for(SizeType j=0; j!=n; ++j) {
        if(!is_inf(c[j].upper_bound())) { x[j]=1; ++count; }
        if(!is_inf(c[j].lower_bound())) { x[n+j]=1; ++count; }
    }
    for(SizeType k=2*n; k!=l; ++k) { x[k]=1; }
    if(count!=0u) { for(SizeType k=0; k!=l; ++k) { x[k]/=count; } }
}

Void NonlinearInteriorPointSolver::
setup_feasibility(const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
                  FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const
{
    this->setup_feasibility(d,g,c,x,y);
    t=0;
    const auto slack=feasibility_slack(d,g,c,y,t);
    for(SizeType k=0; k!=slack.size(); ++k) { t=min(t,slack[k]); }
    // An absolute unit offset alone can round away for large violations.
    t-=max(FloatDPApproximation(1,dp),abs(t));
}

} // namespace Ariadne
