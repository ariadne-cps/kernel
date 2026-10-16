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
        if(i%3==0 && i<=10) { mu *= 0.25_exact; }
    }

    return ValidatedVectorType(exact_vector(x));
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



ValidatedKleenean NonlinearInteriorPointSolver::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("d="<<d<<", g="<<g<<", c="<<c);

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());
    FloatDPApproximation t(dp);
    FloatDPApproximationVector x,y,z;

    this->setup_feasibility(d,g,c,x,y);

    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        LOGGING_PRINTLN_AT(1,"t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z);
        this->feasibility_step(d,g,c,x,y);
        if(probably(LogicalValue(t>0))) {
            LOGGING_PRINTLN_AT(1,"y="<<y<<", g(y)="<<g(y));
            if(this->is_feasible_point(d,g,c,exact_vector(y))) {
                return true;
            }
        }
    }
    LOGGING_PRINTLN("t="<<t<<", y="<<y<<", g(y)="<<g(y));
    if(this->is_infeasibility_certificate(d,g,c,exact_vector(x))) {
        return false;
    }
    return indeterminate;
}


Void
NonlinearInteriorPointSolver::feasibility_step(
    const ExactBoxType&, const ApproximateVectorMultivariateFunction&, const ExactBoxType&,
    FloatDPApproximationVector&, FloatDPApproximationVector&) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Void
NonlinearInteriorPointSolver::feasibility_step(
    const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const
{
    LOGGING_SCOPE_CREATE;

    static const ExactDouble GAMMA=0.0009765625_x; // 1.0/1024;
    static const ExactDouble SIGMA=0.125_x;
    static const ExactDouble SCALE=0.75_x;

    const SizeType m=d.size();
    const SizeType n=c.size();

    FloatDPApproximationVector z(n,dp);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==m);
    ARIADNE_ASSERT(y.size()==n);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    LOGGING_PRINTLN("ddgx="<<ddgx);

    Vector<FloatDPApproximation> gx = ddgx.value();
    LOGGING_PRINTLN("g(x)="<<gx<<" ");
    Matrix<FloatDPApproximation> A = transpose(ddgx.jacobian());
    LOGGING_PRINTLN("A="<<A<<" ");

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) \lambda_k$
    Matrix<FloatDPApproximation> H(m,m,dp);
    for(SizeType i=0; i!=m; ++i) {
        H+=y[i]*ddgx[i].hessian();
    }
    LOGGING_PRINTLN("H="<<H<<" ");


    // Add correction for singleton domain to diagonal elements of Hessian
    for(SizeType i=0; i!=m; ++i) {
    }

    // Compute diagonal entries of KKT Hessian
    Vector<FloatDPApproximation> D(n,dp);
    for(SizeType j=0; j!=n; ++j) {
        if(c[j].lower_bound()==c[j].upper_bound()) {
        } else if(c[j].upper_bound().raw()==+inf) {
        } else if(c[j].lower_bound().raw()==-inf) {
        } else {
            ARIADNE_DEBUG_ASSERT(definitely(-infty<c[j].lower_bound() && c[j].lower_bound()<c[j].upper_bound() && c[j].upper_bound()<+infty));
        }
    }

    FloatDPApproximation sigma(SIGMA,dp);
    FloatDPApproximation mu=dot(x,z)/m;
    if(!egtr(emul(x,z),GAMMA*mu)) {
        LOGGING_PRINTLN("WARNING: near-degeneracy in Lyapunov multipliers in interior-point solver:\n  x="<<x<<", y="<<y<<", z="<<z);
        x=(1-sigma)*x+FloatDPApproximationVector(x.size(),sigma/x.size());
        mu=dot(x,z)/m;
    }

    FloatDPApproximationVector yt=join(y,t);
    LOGGING_PRINTLN("m="<<m<<" n="<<n);
    LOGGING_PRINTLN("x="<<x<<" yt="<<yt<<" z="<<z);


    // Construct diagonal matrices
    FloatDPApproximationVector DE=ediv(x,z);
    LOGGING_PRINTLN("D="<<DE);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    FloatDPApproximationVector gye(2*(m+n),dp);
    //for(SizeType j=0; j!=n; ++j) { gxe[j]=gy[j]-c[j].upper_bound()+t; gye[n+j]=c[j].lower_bound()-gy[j]+t; }
    //for(SizeType i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper_bound()+t; gye[2*n+m+i]=d[i].lower_bound()-y[i]+t; }
    LOGGING_PRINTLN("GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    FloatDPApproximationMatrix AE(m+1,2*(m+n),dp);
    //for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    //for(SizeType i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    //for(SizeType k=0; k!=o; ++k) { AE[m][k]=1; }
    FloatDPApproximationMatrix AET=transpose(AE);

    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //LOGGING_PRINTLN("S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //LOGGING_PRINTLN("S="<<S);
    FloatDPApproximationMatrix S=feasibility_adat(H,A,DE);
    LOGGING_PRINTLN("S="<<S);
    FloatDPApproximationMatrix Sinv=inverse(S);
    LOGGING_PRINTLN("Sinv="<<Sinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    FloatDPApproximationVector rx=esub(emul(x,z),mu*sigma);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    FloatDPApproximationVector ryt=-feasibility_mul(A,x); ryt[m]+=1; // FIXME: Need hessian
    FloatDPApproximationVector rz=gye+z;
    LOGGING_PRINTLN("rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;
    FloatDPApproximationVector rr=ryt + AE*ediv(FloatDPApproximationVector(rx-emul(x,rz)),z) - ryt;


    // Compute the differences
    FloatDPApproximationVector dyt=Sinv*rr;
    //RawFloatDPVector dz=-rz-prod(AET,dyt);
    FloatDPApproximationVector dz=-rz-feasibility_trmul(A,dyt);
    FloatDPApproximationVector dx=-ediv(FloatDPApproximationVector(rx+emul(x,dz)),z);
    LOGGING_PRINTLN("dx="<<dx<<" dyt="<<dyt<<" dz="<<dz);

    FloatDPApproximationVector nx,ny,nyt,nz; FloatDPApproximation nt(dp);

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    Bool allpositive=false;
    FloatDPApproximation alpha=1/FloatDPApproximation(SCALE,dp);
    if(!egtr(emul(x,z) , GAMMA*mu/16)) {
        LOGGING_PRINTLN("WARNING: x="<<x<<", z="<<z<< ", x.z="<<emul(x,z)<<"<"<<GAMMA*mu / 16);
        throw NearBoundaryOfFeasibleDomainException();
    }
    while(!allpositive) {
        alpha=alpha*SCALE;
        nx=x+alpha*dx;
        nyt=yt+alpha*dyt;
        ny=project(nyt,range(0,m));
        nt=nyt[m];
        //NonlinearInteriorPointSolver::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),GAMMA*mu);
    }
    LOGGING_PRINTLN("alpha="<<alpha);
    LOGGING_PRINTLN("nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<emul(nx,nz));

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];
}

FloatDPApproximation NonlinearInteriorPointSolver::
compute_mu(const ExactBoxType&, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
           const FloatDPApproximationVector& x, const FloatDPApproximationVector& lambda) const
{
    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions
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
setup_feasibility(const ExactBoxType& d, const ApproximateVectorMultivariateFunction&, const ExactBoxType& c,
                  FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    const SizeType l=2*(d.size()+c.size());
    y=midpoint(d);
    x=FloatDPApproximationVector(l,one/l);
    //compute_tz(d,g,c,y,t,z);
}





} // namespace Ariadne
