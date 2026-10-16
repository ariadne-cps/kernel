/***************************************************************************
 *            solving/nonlinear_infeasible_interior_point_method.cpp
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

//------- NonlinearInfeasibleInteriorPointSolver -------------------------//

struct NonlinearInfeasibleInteriorPointSolver::PrimalDualData {
    PrimalDualData() : PrimalDualData(0u,0u,dp) { }
    PrimalDualData(SizeType m, SizeType n, DP pr) : w(m,pr), x(n,pr), y(m,pr) { }
    FloatDPApproximationVector w,x,y;
};

struct NonlinearInfeasibleInteriorPointSolver::StepData : public PrimalDualData {
    StepData() : StepData(0u,0u,dp) { }
    StepData(SizeType m, SizeType n, DP pr)
        : PrimalDualData(m,n,pr), vl(m,pr), wl(m,pr), xl(n,pr), zl(n,pr), vu(m,pr), wu(m,pr), xu(n,pr), zu(n,pr), mu(pr) { }
    FloatDPApproximationVector vl,wl,xl,zl,vu,wu,xu,zu; FloatDPApproximation mu;
};

ValidatedVectorType NonlinearInfeasibleInteriorPointSolver::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("f="<<f<<", D="<<D<<", g="<<g<<", C="<<C);

    static const double VALUE_TOLERANCE=1e-8;
    static const double STATE_TOLERANCE=1e-8;
    static const CounterType MAXIMUM_STEPS=24;

    ARIADNE_ASSERT(f.argument_size()==D.size());
    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());
    StepData v;
    FloatDPApproximationVector& x=cast_approximate(v.x);
    FloatDPApproximationVector& y=cast_approximate(v.y);
    C=intersection(C,cast_exact_box(apply(g,D)+UpperIntervalVectorType(C.size(),UpperIntervalType(-1,+1))));
    this->setup_feasibility(D,g,C,v);
    FloatDPApproximationVector oldx=x;

    static const ExactDouble MU_MIN = 1e-12_pr;

    // FIXME: Allow more steps
    for(SizeType i=0; i!=MAXIMUM_STEPS; ++i) {
        LOGGING_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        oldx=x;
        FloatDPApproximation oldfx=f(oldx);
        this->step(f,D,g,C,v);
        if(this->is_infeasibility_certificate(D,g,C,exact_vector(y))) {
            LOGGING_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            LOGGING_PRINTLN_AT(1,"Infeasible");
            throw InfeasibleProblemException();
        }
        FloatDPApproximation fx=f(x);
        if(probably(mag(fx-oldfx)<VALUE_TOLERANCE) && probably(norm(oldx-x)<STATE_TOLERANCE)) {
            break;
        }
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    LOGGING_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));

    if(this->validate_feasibility(D,g,C,exact_vector(x))) {
        LOGGING_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        return exact_vector(x);
    }
    LOGGING_PRINTLN("indeterminate_feasibility");
    throw IndeterminateFeasibilityException();
}

ValidatedKleenean NonlinearInfeasibleInteriorPointSolver::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    return this->feasible_candidate(D,g,C).first;
}

Pair<ValidatedKleenean,FloatDPApproximationVector>
NonlinearInfeasibleInteriorPointSolver::
feasible_candidate(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    LOGGING_SCOPE_CREATE
    LOGGING_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);

    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());

    StepData v;
    FloatDPApproximationVector& x=cast_approximate(v.x);
    FloatDPApproximationVector& y=cast_approximate(v.y);

    auto safe_candidate=[&D](FloatDPApproximationVector const& candidate) -> FloatDPApproximationVector {
        if(contains(D,exact_vector(candidate))) {
            return candidate;
        }
        return midpoint(D);
    };

    ApproximateScalarMultivariateFunction f(EuclideanDomain(D.dimension()));
    ExactBoxType R=intersection(cast_exact_box(widen(apply(g,D),1)),C);
    this->setup_feasibility(D,g,R,v);

    static const ExactDouble MU_MIN = 1e-12_pr;

    for(SizeType i=0; i!=12; ++i) {
        LOGGING_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        try {
            this->step(f,D,g,R,v);
        } catch(const NearBoundaryOfFeasibleDomainException&) {
            LOGGING_PRINTLN("Near boundary of feasible domain; returning indeterminate candidate");
            return {indeterminate,safe_candidate(x)};
        }
        if(this->validate_feasibility(D,g,C,exact_vector(x))) {
            LOGGING_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            LOGGING_PRINTLN("Feasible");
            return {true,x};
        }
        if(this->is_infeasibility_certificate(D,g,C,exact_vector(y))) {
            LOGGING_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            LOGGING_PRINTLN("Infeasible");
            return {false,safe_candidate(x)};
        }
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    LOGGING_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    LOGGING_PRINTLN("Indeterminate");
    return {indeterminate,safe_candidate(x)};
}

Void NonlinearInfeasibleInteriorPointSolver::
setup_feasibility(const ExactBoxType& D, const ApproximateVectorMultivariateFunction&, const ExactBoxType& C,
                  StepData& v) const
{
    ExactIntervalType I(-1,+1);
    SizeType m=C.size(); SizeType n=D.size();

    v.x=midpoint(D);
    v.y=FloatDPApproximationVector(m,zero);
    v.w=midpoint(C);

    //stp.xl=lower_bound(D)-x;
    v.wl=Vector(m,-one);
    v.wu=Vector(m,+one);
    v.xl=lower_bounds(D)-v.x;
    v.xu=upper_bounds(D)-v.x;
    v.vl=Vector(m,-one);
    v.vu=Vector(m,+one);
    v.zl=Vector(n,-one);
    v.zu=Vector(n,+one);
    // FIXME: What should relaxation parameter be?
    v.mu=1.0_x;
}


Void
NonlinearInfeasibleInteriorPointSolver::step(
    const ApproximateScalarMultivariateFunction& f, const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    StepData& v) const
{
    FloatDPApproximationVector& w=v.w; FloatDPApproximationVector& x=v.x; FloatDPApproximationVector& y=v.y;
    FloatDPApproximation& mu=v.mu;
    FloatDPApproximationVector& wl=v.wl; FloatDPApproximationVector& wu=v.wu;
    FloatDPApproximationVector& xl=v.xl; FloatDPApproximationVector& xu=v.xu;
    FloatDPApproximationVector& vl=v.vl; FloatDPApproximationVector& vu=v.vu;
    FloatDPApproximationVector& zl=v.zl; FloatDPApproximationVector& zu=v.zu;
    FloatDPApproximationVector cl=lower_bounds(c); FloatDPApproximationVector cu=upper_bounds(c);
    FloatDPApproximationVector dl=lower_bounds(d); FloatDPApproximationVector du=upper_bounds(d);

    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    LOGGING_PRINTLN("w ="<<w<<",  x ="<<x<<", y ="<<y<<" mu="<<mu);
    LOGGING_PRINTLN("wl="<<wl<<", wu="<<wu<<", xl="<<xl<<", xu="<<xu);
    LOGGING_PRINTLN("vl="<<vl<<", vu="<<vu<<", zl="<<zl<<", zu="<<zu);
    LOGGING_PRINTLN("cl-wl="<<cl-wl<<", dl-xl="<<dl-xl);
    LOGGING_PRINTLN("  w  ="<<w<<",   x  ="<<x);
    LOGGING_PRINTLN("cu-wu="<<cu-wu<<", du-xu="<<du-xu);
    static const ExactDouble gamma=0.0009765625_x;
    static const ExactDouble sigma=0.125_x;
    static const ExactDouble scale=0.75_x;

    const SizeType n=d.size();
    const SizeType m=c.size();

    ARIADNE_ASSERT_MSG(f.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.result_size()==c.size(),  "f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    mu = mu * sigma;

    FloatDPApproximationVector ax(x);
    FloatDPApproximationDifferential ddfx=f.differential(ax,2u);
    LOGGING_PRINTLN("ddfx="<<ddfx);
    Vector<FloatDPApproximationDifferential> ddgx=g.differential(ax,2u);
    LOGGING_PRINTLN("ddgx="<<ddgx);

    FloatDPApproximation fx = ddfx.value();
    Vector<FloatDPApproximation> gx = ddgx.value();
    LOGGING_PRINTLN("f(x)="<<fx);
    LOGGING_PRINTLN("g(x)="<<gx);
    Vector<FloatDPApproximation> Jfx = transpose(ddfx.gradient());
    Matrix<FloatDPApproximation> A = ddgx.jacobian();
    Matrix<FloatDPApproximation>& Jgx = A;
    LOGGING_PRINTLN("Df(x)="<<Jfx);
    LOGGING_PRINTLN("Dg(x)="<<Jgx);

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<FloatDPApproximation> YH = ddfx.hessian();
    for(SizeType i=0; i!=m; ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    LOGGING_PRINTLN("D2f(x)="<<ddfx.hessian());
    LOGGING_PRINTLN("D2f(x)+Y.D2g(x)="<<YH);

    // Set up the system of equations
    // (A^TDA + E - Y.H) dx = A^T(r_w-Dr_y)+r_x
    // dw = A \delta x + r_y
    // dy = r_w - D dw

    FloatDPApproximationDiagonalMatrix const& Vl=diagonal_matrix(vl);
    FloatDPApproximationDiagonalMatrix const& Vu=diagonal_matrix(vu);
    FloatDPApproximationDiagonalMatrix const& Wl=diagonal_matrix(wl);
    FloatDPApproximationDiagonalMatrix const& Wu=diagonal_matrix(wu);
    FloatDPApproximationDiagonalMatrix const& Xl=diagonal_matrix(xl);
    FloatDPApproximationDiagonalMatrix const& Xu=diagonal_matrix(xu);
    FloatDPApproximationDiagonalMatrix const& Zl=diagonal_matrix(zl);
    FloatDPApproximationDiagonalMatrix const& Zu=diagonal_matrix(zu);

    // Compute the diagonal matrices
    //   D=XL/ZL+XU/ZU  E=WL/VL+WU/VU
    FloatDPApproximationDiagonalMatrix Dl=Vl/Wl;
    FloatDPApproximationDiagonalMatrix Du=Vu/Wu;
    FloatDPApproximationDiagonalMatrix D=Dl+Du;
    LOGGING_PRINTLN("D="<<D);
    FloatDPApproximationDiagonalMatrix El=Zl/Xl;
    FloatDPApproximationDiagonalMatrix Eu=Zu/Xu;
    FloatDPApproximationDiagonalMatrix E=El+Eu;
    LOGGING_PRINTLN("E="<<E);

    // Construct the residuals
    // The residual for the slack variable xl is given by the duality condition xl.zl=mu as mu/xl-zl
    // The residual for the dual variable zl is given by the slackness condition x-xl-cl
    // The residual for the auxiliary variable w is given by y-(vu-vl)
    // The residual for the dual variable y is given by g(x)-w
    FloatDPApproximationVector ew=(vl+vu)-y;
    FloatDPApproximationVector ex=Jfx+transpose(Jgx)*y+(zl+zu);
    FloatDPApproximationVector ey=gx-w;
    FloatDPApproximationVector ewl=esub(vl,ediv(mu,wl));
    FloatDPApproximationVector ewu=esub(vu,ediv(mu,wu));
    FloatDPApproximationVector exl=esub(zl,ediv(mu,xl));
    FloatDPApproximationVector exu=esub(zu,ediv(mu,xu));
    FloatDPApproximationVector evl=w+wl-cl;
    FloatDPApproximationVector evu=w+wu-cu;
    FloatDPApproximationVector ezl=x+xl-dl;
    FloatDPApproximationVector ezu=x+xu-du;

    LOGGING_PRINTLN("ew="<<ew<<", ex="<<ex<<", ey="<<ey);
    LOGGING_PRINTLN("ewl="<<ewl<<", ewu="<<ewu<<", exl="<<exl<<" exu="<<exu);
    LOGGING_PRINTLN("evl="<<evl<<", evu="<<evu<<", ezl="<<ezl<<" ezu="<<ezu);

    FloatDPApproximationVector rw = ew - (ewl+ewu) + Dl*evl + Du*evu;
    FloatDPApproximationVector rx = ex - (exl+exu) + El*ezl + Eu*ezu;
    FloatDPApproximationVector& ry = ey;

    // Solve linear system
    // ( D   0  -I ) (dw)   (rw)
    // ( 0  H+E A^T) (dx) = (rx)
    // (-I   A   0 ) (dy) = (ry)

    // Eliminate dw=A*dx-ry and dy=D*dw-rw
    // Reduce to (H+E)*dx+A'*(D*(A*dx-ry)-rw)=rx
    // Simplify to (H+E+A'*D*A)*dx=rx+A'*(D*ry+rw)

    // normal equation matrix
    FloatDPApproximationMatrix S=YH;
    atda(S,A,D);
    S+=E;

    LOGGING_PRINTLN("S="<<S);
    ARIADNE_DEBUG_ASSERT(decide(norm(FloatDPApproximationMatrix(S-(YH+E+transpose(A)*(D*A))))/norm(S)<1e-8));
    LOGGING_PRINTLN("Sinv="<<inverse(S));

    FloatDPApproximationVector r = transpose(A)*(rw+D*ry)+rx;
    LOGGING_PRINTLN("rw="<<rw<<" rx="<<rx<<" ry="<<ry);
    LOGGING_PRINTLN("r="<<r);

    // Compute the differences
    FloatDPApproximationVector dx = solve(S,r);
    // Apply correction to improve accuracy
    dx = dx + solve(S,r-S*dx);
    LOGGING_PRINTLN("S*dx="<<S*dx<<" r="<<r);
    LOGGING_PRINTLN("S*inverse(S)-I="<<S*inverse(S)-FloatDPApproximationMatrix::identity(n,dp));
    ARIADNE_DEBUG_ASSERT(decide(norm(S*dx - r)/max(1.0_x,norm(r))<1e-4));

    FloatDPApproximationVector dw = A*dx-ry;
    FloatDPApproximationVector dy = D*dw-rw;
    LOGGING_PRINTLN("dw="<<dw<<" dx="<<dx<<" dy="<<dy);

    LOGGING_PRINTLN("YH*dx+E*dx+dy*A="<<(YH*dx+E*dx+transpose(A)*dy)<<", rx="<<rx);

    // Check solution of linear system for residuals
    ARIADNE_DEBUG_ASSERT(decide(norm(D*dw-dy-rw)/max(one,norm(rw))<1e-4));
    ARIADNE_DEBUG_ASSERT(decide(norm(YH*dx+E*dx+transpose(A)*dy-rx)/max(one,norm(rx))<1e-2));
    ARIADNE_DEBUG_ASSERT(decide(norm(-dw+A*dx-ry)/max(one,norm(ry))<1e-4));


    FloatDPApproximationVector dwl = evl-dw;
    FloatDPApproximationVector dwu = evu-dw;
    FloatDPApproximationVector dxl = ezl-dx;
    FloatDPApproximationVector dxu = ezu-dx;
    FloatDPApproximationVector dvl = ewl-Dl*dwl;
    FloatDPApproximationVector dvu = ewu-Du*dwu;
    FloatDPApproximationVector dzl = exl-El*dxl;
    FloatDPApproximationVector dzu = exu-Eu*dxu;

    LOGGING_PRINTLN("dwl="<<dwl<<", dwu="<<dwu<<", dxl="<<dxl<<" dxu="<<dxu);
    LOGGING_PRINTLN("dvl="<<dvl<<", dvu="<<dvu<<", dzl="<<dzl<<" dzu="<<dzu);

    LOGGING_PRINTLN("YH*dx+dy*A+dzl+dzu="<<(YH*dx+transpose(A)*dy+dzl+dzu)<<", ex="<<ex);
    // Check solution of linear system
/*
    ARIADNE_DEBUG_ASSERT(norm(-dy+dvl+dvu - ew)/max(1.0,norm(ew))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(YH*dx+dy*A+dzl+dzu - ex)/max(1.0,norm(ex))<1e-2);
    ARIADNE_DEBUG_ASSERT(norm(-dw+A*dx - ey)/max(1.0,norm(ey))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(Dl*dwl+dvl - ewl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(Du*dwu+dvu - ewu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(El*dxl+dzl - exl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(Eu*dxu+dzu - exu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dw+dwl - evl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dw+dwu - evu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dx+dxl - ezl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dx+dxu - ezu)<1e-12);
*/

    FloatDPApproximationVector nw; FloatDPApproximationVector nx; FloatDPApproximationVector ny;
    FloatDPApproximationVector nwl; FloatDPApproximationVector nwu; FloatDPApproximationVector nxl; FloatDPApproximationVector nxu;
    FloatDPApproximationVector nvl; FloatDPApproximationVector nvu; FloatDPApproximationVector nzl; FloatDPApproximationVector nzu;


    FloatDPApproximation alpha=one;
    nx = x-alpha*dx;
    // Pick an update value which minimises the objective function
    FloatDPApproximation fxmin=f(nx);
    FloatDPApproximation alphamin=one;
    static const CounterType REDUCTION_STEPS=4;
    for(SizeType i=0; i!=REDUCTION_STEPS; ++i) {
        alpha*=scale;
        nx = x-alpha*dx;
        FloatDPApproximation fnx=f(nx);
        if(decide(fnx<fxmin*scale)) {
            fxmin=fnx;
            alphamin=alpha;
        }
    }
    //alpha=alphamin;
    alpha=one;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly.
    static const double MINIMUM_ALPHA=1e-16;
    Bool allfeasible=false;
    while(decide(alpha>MINIMUM_ALPHA) && !allfeasible) {
        nwl=wl-alpha*dwl;
        nwu=wu-alpha*dwu;
        nxl=xl-alpha*dxl;
        nxu=xu-alpha*dxu;
        nvl=vl-alpha*dvl;
        nvu=vu-alpha*dvu;
        nzl=zl-alpha*dzl;
        nzu=zu-alpha*dzu;
        allfeasible = elss(nwl,mu*gamma) && egtr(nwu,mu*gamma) && elss(nxl,mu*gamma) && egtr(nxu,mu*gamma)
                           && elss(nvl,mu*gamma)&& egtr(nvu,mu*gamma) && elss(nzl,mu*gamma) && egtr(nzu,mu*gamma);
        //allfeasible = eneg(nwl) && epos(nwu) && eneg(nxl) && epos(nxu) && elss(nvl,mu*gamma) && egtr(nvu,mu*gamma) && elss(nzl,mu*gamma) && egtr(nzu,mu*gamma);
        if(!allfeasible) { alpha*=scale; }
    }
    nw=w-alpha*dw;
    nx=x-alpha*dx;
    ny=y-alpha*dy;
    if(decide(alpha<=MINIMUM_ALPHA)) {
        LOGGING_PRINTLN_AT(1,"w="<<w<<"  x="<<x<<"  y="<<y);
        LOGGING_PRINTLN_AT(1,"nw="<<nw<<"  nx="<<nx<<"  ny="<<ny);
        throw NearBoundaryOfFeasibleDomainException(); }
    LOGGING_PRINTLN("alpha="<<alpha);
    LOGGING_PRINTLN("nw="<<nw<<" nx="<<nx<<" ny="<<ny);
    LOGGING_PRINTLN("nwl="<<nwl<<", nwu="<<nwu<<", nxl="<<nxl<<" nxu="<<nxu);
    LOGGING_PRINTLN("nvl="<<nvl<<", nvu="<<nvu<<", nzl="<<nzl<<" nzu="<<nzu);

    w=nw; x=nx; y=ny;
    wl=nwl; wu=nwu; xl=nxl; xu=nxu;
    vl=nvl; vu=nvu; zl=nzl; zu=nzu;

    FloatDPApproximation nmu = zero;
    for(SizeType i=0; i!=m; ++i) {
        nmu = nmu + wl[i]*vl[i] + wu[i]*vu[i];
    }
    for(SizeType j=0; j!=n; ++j) {
        nmu = nmu + xl[j]*zl[j] + xu[j]*zu[j];
    }
    nmu /= (2*(m+n));
    mu = nmu;

    LOGGING_PRINTLN("nmu="<<nmu);

}






} // namespace Ariadne
