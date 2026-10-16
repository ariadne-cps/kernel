/***************************************************************************
 *            solving/nonlinear_interval_interior_point_method.cpp
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

// Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
ValidatedKleenean IntervalNonlinearInteriorPointSolver::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("D="<<D<<", h="<<h);

    const SizeType n=D.size();

    FloatDPBoundsVector zl(n,dp), zu(n,dp);
    ExactFloatDPVectorType xl = Ariadne::lower_bounds(D);
    ExactFloatDPVectorType xu = Ariadne::upper_bounds(D);

    FloatDPBoundsVector x=cast_singleton(D);
    FloatDPBoundsVector y(h.result_size(),FloatDPBounds(-1,+1,dp));
    FloatDPBounds mu(0,1,dp);

    for(SizeType i=0; i!=8; ++i) {
        this->feasibility_step(xl,xu,h,x,y,zl,zu,mu);
    }

    return indeterminate;
}

Void IntervalNonlinearInteriorPointSolver::
feasibility_step(const ExactFloatDPVectorType& xl, const ExactFloatDPVectorType& xu, const ValidatedVectorMultivariateFunction& h,
                 FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("[x]="<<x<<" [lambda]="<<y<<", [zl]="<<zl<<", [zu]="<<zu<<" [mu]="<<mu);

    const SizeType n=x.size();
    const SizeType m=y.size();

    FloatDPBoundsVector mx=midpoint(x);
    FloatDPBoundsVector my=midpoint(y);
    FloatDPBoundsVector mzl=midpoint(zl);
    FloatDPBoundsVector mzu=midpoint(zu);
    FloatDPBounds mmu(midpoint(mu));
    LOGGING_PRINTLN_AT(1,"x~"<<x<<" lambda~="<<y<<", mu~"<<mu);

    // Solve equations y Dh(x) - zl + zu = 0; h(x) = 0; (x-xl).zl - mu = 0;  (xu-x).zu - mu = 0; Sum_j y_j^2 - mu = 0
    Vector<FloatDPBoundsDifferential> ddhx=h.evaluate(FloatDPBoundsDifferential::variables(2,x));
    Vector<FloatDPBoundsDifferential> dhmx=h.evaluate(FloatDPBoundsDifferential::variables(1,mx));
    FloatDPBoundsMatrix A = ddhx.jacobian();
    FloatDPBoundsMatrix mA = dhmx.jacobian();
    LOGGING_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPBoundsVector rx = transpose(mA) * my;
    for(SizeType j=0; j!=n; ++j) {
        rx[j] -= mmu*rec(mx[j]-xl[j]);
        rx[j] += mmu*rec(xu[j]-mx[j]);
    }
    FloatDPBoundsVector ry = dhmx.value();
    FloatDPBoundsVector rzl = esub(emul(FloatDPBoundsVector(mx-xl),mzl),mmu);
    FloatDPBoundsVector rzu = esub(emul(FloatDPBoundsVector(xu-mx),mzu),mmu);
    LOGGING_PRINTLN("rx="<<rx<<" ry="<<ry<<" rzl="<<rzl<<" rzu="<<rzu);

    FloatDPBoundsMatrix H(n,n,dp);
    for(SizeType i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(SizeType j=0; j!=n; ++j) {
        H[j][j] += mu*rec(sqr(x[j]-xl[j]));
        H[j][j] += mu*rec(sqr(xu[j]-x[j]));
    }

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPBoundsMatrix Hinv=inverse(H);
    LOGGING_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPBoundsMatrix S=A*Hinv*transpose(A);
    FloatDPBoundsMatrix Sinv=inverse(S);
    LOGGING_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    FloatDPBoundsVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatDPBoundsVector dx = Hinv * ( rx - transpose(A) * dy);
    LOGGING_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPBoundsVector nx = x-dx;
    FloatDPBoundsVector ny = y-dy;
    LOGGING_PRINTLN("nx="<<nx<<" ny="<<ny);
    LOGGING_PRINTLN_AT(1,"h(x)="<<h(nx));

    x = refinement(x,nx); y=refinement(y,ny);
    FloatDPBounds nmu = zero;
    for(SizeType i=0; i!=m; ++i) { nmu += sqr(y[i]); }
    mu=refinement(mu,nmu);
}

/*

struct KuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<KuhnTuckerFunctionBody,ExactIntervalType>
{
    ValidatedScalarMultivariateFunction f;
    Array<ValidatedScalarMultivariateFunction> g;
    Array<ValidatedScalarMultivariateFunction> df;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    KuhnTuckerFunctionBody(ValidatedScalarMultivariateFunction _f, ValidatedVectorMultivariateFunction _g) {
        ARIADNE_ASSERT(_f.argument_size()==_g.argument_size());
        const SizeType m=_g.argument_size();
        const SizeType n=_g.result_size();
        g.resize(n); df.resize(m); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        f=_f;
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; }
        for(SizeType i=0; i!=m; ++i) { df[i]=f.derivative(i); }
        for(SizeType j=0; j!=n; ++j) { for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return g.size()*2+f.argument_size(); }
    SizeType argument_size() const { return g.size()*2+f.argument_size(); }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const SizeType m=f.argument_size();
        const SizeType n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        Vector<X> rx(m), rz(n), rs(n);
        for(SizeType i=0; i!=m; ++i) { rx[i]=df[i].evaluate(y); for(SizeType j=0; j!=n; ++j) { rx[i]=rx[i]-x[j]*dg[j][i].evaluate(y); } }
        for(SizeType j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + z[j]; }
        for(SizeType j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
    }
};

struct FeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    Array<ValidatedScalarMultivariateFunction> g;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    FeasibilityKuhnTuckerFunctionBody(ValidatedVectorMultivariateFunction _g) {
        const SizeType m=_g.argument_size();
        const SizeType n=_g.result_size();
        g.resize(n); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return g.size()*2+g[0].argument_size()+1; }
    SizeType argument_size() const { return g.size()*2+g[0].argument_size()+1; }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const SizeType m=g[0].argument_size();
        const SizeType n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        X t(arg[n+m+n]);
        Vector<X> rx(m), rz(n), rs(n); X rt;
        for(SizeType i=0; i!=m; ++i) { rx[i]=x[0]*dg[0][i].evaluate(y); for(SizeType j=1; j!=n; ++j) { rx[i]=rx[i]+x[j]*dg[j][i].evaluate(y); } }
        for(SizeType j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + t + z[j]; }
        for(SizeType j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        rt=1-x[0]; for(SizeType j=1; j!=n; ++j) { rt=rt-x[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
        res[n+m+n]=rt;
    }
};



struct ConstrainedFeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    SizeType m;
    SizeType n;
    ExactIntervalVectorType d;
    Array<ValidatedScalarMultivariateFunction> g;
    ExactIntervalVectorType c;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    ConstrainedFeasibilityKuhnTuckerFunctionBody(ExactBoxType D, ValidatedVectorMultivariateFunction _g, ExactBoxType C) {
        m=_g.argument_size();
        n=_g.result_size();
        d=D; c=C;
        g.resize(n); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return 5*m+4*n+1u; }
    SizeType argument_size() const { return 5*m+4*n+1u; }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream& os) const { return os << "KuhnTuckerFunctionBody"; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const X zero=arg[0].zero_element();
        const SizeType l=2*(m+n);
        assert(arg.size()==l+m+l+1);
        Vector<X> x(project(arg,range(0u,l)));
        Vector<X> y(project(arg,range(l,l+m)));
        Vector<X> z(project(arg,range(l+m,l+m+l)));
        X t(arg[l+m+l]);
        Vector<X> rx(m,zero), rz(l,zero), rs(l,zero); X rt(zero);
        Vector<X> gy(n);
        for(SizeType j=0; j!=n; ++j) { gy[j]=g[j].evaluate(y); }
        Matrix<X> dgy(n,m);
        for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { dgy[j][i]=dg[j][i].evaluate(y); } }

        for(SizeType i=0; i!=m; ++i) {
            for(SizeType j=0; j!=n; ++j) { rx[i]+=x[j]*(dgy[j][i]-c[j].upper_bound()); rx[i]+=x[n+j]*(c[j].lower_bound()-dgy[j][i]); }
            rx[i]+=x[2*n+i]*(y[i]-d[i].upper_bound())-x[2*n+m+i]*(d[i].lower_bound()-y[i]);
        }
        for(SizeType j=0; j!=n; ++j) { rz[j]=gy[j] + t + z[j]; rz[n+j]=t+z[n+j]-gy[j]; }
        for(SizeType i=0; i!=m; ++i) { rz[2*n+i]=y[i]+t+z[2*n+i]; rz[2*n+m+i]=y[i]+t+z[2*n+m+i]; }
        for(SizeType k=0; k!=l; ++k) { rs[k]=x[k]*z[k]; }
        rt+=1.0; for(SizeType j=0; j!=2*n; ++j) { rt=rt-x[j]; }
        project(res,range(0,l))=rz;
        project(res,range(l,l+m))=rx;
        project(res,range(l+m,l+m+l))=rs;
        res[l+m+l]=rt;
    }
};



ValidatedVectorType KrawczykOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


ValidatedKleenean KrawczykOptimiser::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    LOGGING_PRINTLN_AT(2,"KrawczykOptimiser::feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c)");
    LOGGING_PRINTLN_AT(2,"  d="<<d<<", g="<<g<<", c="<<c);

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());

    ExactIntervalType t; ExactIntervalVectorType x,y,z;
    setup_feasibility(d,g,c,x,y,z,t);


    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        LOGGING_PRINTLN_AT(4,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z);
        try {
            this->feasibility_step(d,g,c,x,y,z,t);
        }
        catch(const SingularMatrixException& e) {
            return indeterminate;
        }
        if(t.lower_bound()>t.upper_bound()) {
            LOGGING_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return indeterminate;
        }
        if(t.lower_bound()>0.0) {
            LOGGING_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            if(this->is_feasible_point(d,g,c,midpoint(y))) {
                return true;
            }
        }
        if(t.upper_bound()<0.0) {
            LOGGING_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return false;
        }
    }
    LOGGING_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
    if(this->is_infeasibility_certificate(d,g,c,midpoint(x))) {
        return false;
    }
    return indeterminate;
}



Void KrawczykOptimiser::setup_feasibility(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                          ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const SizeType m=g.argument_size();
    const SizeType n=g.result_size();
    const SizeType l=2*(m+n);
    x=ExactIntervalVectorType(l, ExactIntervalType(0,1)/l);
    y=d;
    z.resize(2*(m+n));
    compute_tz(d,g,c,y,t,z);
}


Void KrawczykOptimiser::compute_tz(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                   const ExactIntervalVectorType& y, ExactIntervalType& t, ExactIntervalVectorType& z) const
{
    ARIADNE_ASSERT(d.size()>0u);
    //static const double EPS=1.0/8;
    static const float min_float=std::numeric_limits<float>::min();

    const SizeType m=g.argument_size();
    const SizeType n=g.result_size();

    // Compute the image of y under the constraint function
    ExactIntervalVectorType gy=g(y);
    gy+=ExactIntervalVectorType(gy.size(),ExactIntervalType(-min_float,+min_float));
    ExactIntervalVectorType my=midpoint(y);
    ExactIntervalVectorType mgy=g(my);

    // Find the range of possible values of the optimal t
    // This range is too pessimistic
    t=ExactIntervalType(+inf,+inf);
    for(SizeType j=0; j!=n; ++j) {
        t=min(t,c[j]-gy[j]);
        t=min(t,gy[j]-c[j]);
    }
    for(SizeType i=0; i!=m; ++i) {
        t=min(t,d[i]-y[i]);
        t=min(t,y[i]-d[i]);
    }

    // Find the range of possible values of the optimal t
    FloatDP tmin=+inf;
    FloatDP tmax=+inf;
    for(SizeType j=0; j!=n; ++j) {
        tmax=min(tmax,sub(up,c[j].upper_bound(),gy[j].lower_bound()));
        tmax=min(tmax,sub(up,gy[j].upper_bound(),c[j].lower_bound()));
        tmin=min(tmin,sub(down,c[j].upper_bound(),mgy[j].upper_bound()));
        tmin=min(tmin,sub(down,mgy[j].lower_bound(),c[j].lower_bound()));
    }
    for(SizeType i=0; i!=m; ++i) {
        tmin=min(tmin,sub(up,d[i].upper_bound(),y[i].lower_bound()));
        tmax=min(tmax,sub(up,y[i].upper_bound(),d[i].lower_bound()));
        tmin=min(tmin,sub(down,d[i].upper_bound(),my[i].upper_bound()));
        tmax=min(tmax,sub(down,my[i].lower_bound(),d[i].lower_bound()));
    }
    tmin-=0.0625;
    t=ExactIntervalType(tmin,tmax);


    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(SizeType j=0; j!=n; ++j) {
        z[j]=max(c[j].upper_bound()-gy[j]-t,0.0);
        z[n+j]=max(gy[j]-c[j].lower_bound()-t,0.0);
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=max(d[i].upper_bound()-y[i]-t,0.0);
        z[2*n+m+i]=max(y[i]-d[i].lower_bound()-t,0.0);
    }

    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(SizeType j=0; j!=n; ++j) {
        z[j]=ExactIntervalType(0.0,c[j].upper_bound()-mgy[j].lower_bound()-tmin);
        z[n+j]=ExactIntervalType(0.0,mgy[j].upper_bound()-c[j].lower_bound()-tmin);
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=ExactIntervalType(0.0,d[i].upper_bound()-my[i].lower_bound()-tmin);
        z[2*n+m+i]=ExactIntervalType(0.0,my[i].upper_bound()-d[i].lower_bound()-tmin);
    }

    LOGGING_PRINTLN_AT(9,"  d="<<d<<", c="<<c<<", y="<<y<<", g(y)="<<gy<<", t="<<t<<", z="<<z);

}


Void KrawczykOptimiser::
minimisation_step(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g,
                  ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z) const
{
    const SizeType m=f.argument_size();
    const SizeType n=g.result_size();

    Differential<UpperIntervalType> ddf=f.evaluate(Differential<UpperIntervalType>::variables(2,y));
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    ExactIntervalMatrixType H(m,m);
    set_hessian(H,ddf);
    for(SizeType j=0; j!=n; ++j) { add_hessian(H,-x[j],ddg[j]); }

    ExactIntervalMatrixType A(m,n);
    set_jacobian_transpose(A,ddg);

    LOGGING_PRINTLN_AT(9,"f="<<f<<"\ng="<<g<<"\nx="<<x<<" y="<<y<<" z="<<z);
    LOGGING_PRINTLN_AT(9,"A="<<A<<"\nH="<<H);

    ARIADNE_NOT_IMPLEMENTED;

}



Void KrawczykOptimiser::feasibility_step(const ValidatedVectorMultivariateFunction& g,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    ARIADNE_NOT_IMPLEMENTED;
    const SizeType m=y.size();
    const SizeType n=x.size();

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    LOGGING_PRINTLN_AT(9,"A="<<A);

    // H is the Hessian matrix Hik = xj*dgj/dyidyk
    ExactIntervalMatrixType H(m,m);
    for(SizeType j=0; j!=n; ++j) {
        add_hessian(H,x[j],ddg[j]);
    }
    LOGGING_PRINTLN_AT(9," H="<<H);

    FloatDPMatrix mA=midpoint(A);
    LOGGING_PRINTLN_AT(9," mA="<<mA);
    FloatDPMatrix mH=midpoint(H);
    LOGGING_PRINTLN_AT(9," mH="<<mH);

    RawFloatDPVector mD(n);
    for(SizeType j=0; j!=n; ++j) { mD[j]=midpoint(x[j])/midpoint(z[j]); }
    LOGGING_PRINTLN_AT(9," mD="<<mD);

    FloatDPMatrix& mS=mH;
    adat(mS,mA,mD);
    LOGGING_PRINTLN_AT(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    LOGGING_PRINTLN_AT(9,"mSinv="<<mSinv);
}

// Feasibility step for dual (inequality constrained) problem without using slack variables
// FIXME: Do we need a slackness parameter mu? Probably not; hopefully the infinities are kept in check...
// This method has the advantage of not needing to update the primal variables
Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& y, ExactIntervalType& t) const
{
    const SizeType m=d.size();
    const SizeType n=c.size();

    // Compute function values
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // gy is the vector of values of g(y)
    ExactIntervalVectorType gy(n);
    for(SizeType j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }

    // z is the vector of slack variables z[k]=cu[k]-gy[k]-t or z[k]=gy[k]-cl[k]-t
    ExactIntervalVectorType z(2*(m+n));
    for(SizeType j=0; j!=n; ++j) { z[j]=d[j].upper_bound()-gy[j]-t; z[n+j]=gy[j]-d[j].lower_bound()-t; }
    for(SizeType i=0; i!=m; ++i) { z[i]=c[2*n+i].upper_bound()-y[i]-t; z[2*n+m+i]=y[i]-c[i].lower_bound()-t; }

    ExactIntervalVectorType zr(2*(m+n));
    for(SizeType k=0; k!=2*(m+n); ++k) { zr[k]=1.0/z[k]; }

    ExactIntervalVectorType D(2*(m+n));
    for(SizeType k=0; k!=2*(m+n); ++k) { D[k]=zr[k]*zr[k]; }

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { A[i][j]=ddg[j][i]; } }

    // A is the sum of scaled Hessian matrices hi1i2=zj*ddgj/dyi1yi2
    ExactIntervalMatrixType H(m,m);

    ExactIntervalMatrixType SE(m+1,m+1);
    // SE[0:m][0:m] is the matrix H+/-A(D1+D2)AT+(D3+D4) where D=z.z
    for(SizeType i1=0; i1!=m; ++i1) { for(SizeType i2=0; i2!=m; ++i2) { SE[i1][i2]=H[i1][i2];
        for(SizeType j=0; j!=n; ++j) { SE[i1][i2]+=A[i1][j]*(D[j]+D[n+j])*A[i2][j]; }
    } }
    for(SizeType i=0; i!=m; ++i) { SE[i][i]+=(D[2*n+i]+D[2*n+m+i]); }
    // SE[m][0:m]=SE[0:m][m] is the vector A(D1-D2)e+(D3-D4)e
    for(SizeType i=0; i!=m; ++i) { SE[i][m]=D[2*n+i]-D[2*n+m+i];
        for(SizeType j=0; j!=n; ++j) { SE[i][m]+=A[i][j]*(D[j]-D[n+j]); }
        SE[m][i]=SE[i][m];
    }
    // SE[m][m] is the scalar eT(D1+D2)e+eT(D3+D4)e
    for(SizeType k=0; k!=2*(m+n); ++k) { SE[m][m]+=D[k]; }

    // Vector of residuals
    ExactIntervalVectorType re(m+1);
    for(SizeType i=0; i!=m; ++i) { re[i]+=(zr[2*n+i]-zr[2*n+m+i]);
        for(SizeType j=0; j!=n; ++j) { re[i]+=A[i][j]*(zr[j]-zr[n+j]); }
    }
    for(SizeType j=0; j!=n; ++j) { re[m]+=(zr[j]+zr[n+n]); }
    for(SizeType i=0; i!=m; ++i) { re[m]+=(zr[2*n+i]+zr[2*n+m+i]); }

    // Compute inverse Jacobian matrix
    ExactIntervalMatrixType JE;
    try {
        JE=inverse(midpoint(SE));
    }
    catch(const SingularMatrixException& e) {
        ARIADNE_WARN("Matrix S="<<midpoint(SE)<<" is not invertible");
        LOGGING_PRINTLN_AT(1,"WARNING: Matrix S="<<midpoint(SE)<<" is not invertible");
        throw e;
    }

    // Krawczyk step
    ExactIntervalVectorType dyt=prod(JE,ExactIntervalVectorType(midpoint(re)))+prod(ExactIntervalMatrixType::identity(m+1)-prod(JE,SE),re-midpoint(re));

    // Extract y and t
    ExactIntervalVectorType yt=join(y,t);
    ExactIntervalVectorType nyt=yt+dyt;

    yt=intersection(yt,nyt);
    y=project(yt,range(0,m));
    t=yt[m];
}




Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const SizeType m=d.size();
    const SizeType n=c.size();
    const SizeType o=2*(m+n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==o);
    ARIADNE_ASSERT(y.size()==m);
    ARIADNE_ASSERT(z.size()==o);

    ExactIntervalVectorType yt=join(y,t);
    LOGGING_PRINTLN_AT(9,"m="<<m<<" n="<<n);
    LOGGING_PRINTLN_AT(9,"x="<<x<<" yt="<<yt<<" z="<<z);

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));
    LOGGING_PRINTLN_AT(9,"  ddg="<<ddg);

    // gy is the vector of values of g(y)
    UpperIntervalVectorType gy(n); for(SizeType j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }
    LOGGING_PRINTLN_AT(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    UpperIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    LOGGING_PRINTLN_AT(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    UpperIntervalMatrixType H(m,m);
    for(SizeType j=0; j!=n; ++j) {
        add_hessian(H,x[j]-x[n+j],ddg[j]);
    }
    LOGGING_PRINTLN_AT(9," H="<<H);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    UpperIntervalVectorType gye(o);
    for(SizeType j=0; j!=n; ++j) { gye[j]=gy[j]-c[j].upper_bound()+t; gye[n+j]=c[j].lower_bound()-gy[j]+t; }
    for(SizeType i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper_bound()+t; gye[2*n+m+i]=d[i].lower_bound()-y[i]+t; }
    LOGGING_PRINTLN_AT(9,"  GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    UpperIntervalMatrixType AE(m+1,o);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    for(SizeType i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    for(SizeType k=0; k!=o; ++k) { AE[m][k]=1; }
    UpperIntervalMatrixType AET=transpose(AE);

    FloatDPMatrix mA=midpoint(A);
    FloatDPMatrix mAE=midpoint(AE);
    FloatDPMatrix mAET=midpoint(AET);
    FloatDPMatrix mH=midpoint(H);
    RawFloatDPVector mx=midpoint(x);
    RawFloatDPVector myt=midpoint(yt);
    RawFloatDPVector mz=midpoint(z);
    RawFloatDPVector mDE=ediv(mx,mz);


    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //LOGGING_PRINTLN_AT(9,"S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //LOGGING_PRINTLN_AT(9,"S="<<S);
    FloatDPMatrix mS=feasibility_adat(mH,mA,mDE);
    LOGGING_PRINTLN_AT(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    LOGGING_PRINTLN_AT(9,"mSinv="<<mSinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    UpperIntervalVectorType rx=emul(mx,mz);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType ryt=-feasibility_mul(mA,mx); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType rz=midpoint(gye)+mz;
    LOGGING_PRINTLN_AT(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    // Construct the errors on the residuals ([M]-M)([x]-x)
    UpperIntervalVectorType ex=x-mx;
    UpperIntervalVectorType eyt=yt-myt;
    UpperIntervalVectorType ez=z-mz;
    UpperIntervalMatrixType eA=A-mA;
    UpperIntervalMatrixType eH=H-mH;

    UpperIntervalVectorType erx=2.0*emul(ex,ez);
    UpperIntervalVectorType eryt=UpperIntervalMatrixType(AE-mAE)*ex;
    UpperIntervalVectorType erz=UpperIntervalMatrixType(AET-mAET)*eyt;
    LOGGING_PRINTLN_AT(9,"erx="<<erx<<" eryt="<<eryt<<" erz="<<erz);

    rx+=2.0*emul(ex,ez);
    ryt+=UpperIntervalMatrixType(AE-mAE)*ex;
    rz+=UpperIntervalMatrixType(AET-mAET)*eyt;
    LOGGING_PRINTLN_AT(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;

    // Compute the error differences
    UpperIntervalVectorType erxdz=ediv(erx,mz);
    UpperIntervalVectorType edyt=(mSinv*mAE)*erxdz + mSinv*eyt - (mSinv*(mAE*DiagonalMatrix<FloatDP>(mDE))) * ez;
    UpperIntervalVectorType edz=-erz-feasibility_trmul(mA,edyt);
    UpperIntervalVectorType edx=-ediv(UpperIntervalVectorType(erx+emul(mx,edz)),mz);
    LOGGING_PRINTLN_AT(9,"edx="<<edx<<" edyt="<<edyt<<" edz="<<edz);

    // Compute the error differences
    UpperIntervalVectorType eerr=prod(mAE,ediv(esub(erx,emul(mx,erz)),mz))-eryt;
    LOGGING_PRINTLN_AT(9,"  eerr="<<eerr);
    UpperIntervalVectorType eedyt=prod(mSinv,eerr);
    UpperIntervalVectorType eedz=-erz-feasibility_trmul(mA,eedyt);
    UpperIntervalVectorType eedx=-ediv(UpperIntervalVectorType(erx+emul(mx,eedz)),mz);
    LOGGING_PRINTLN_AT(9,"eedx="<<eedx<<" eedyt="<<eedyt<<" eedz="<<eedz);


    // Compute the differences
    UpperIntervalVectorType rr=prod(mAE,ediv(esub(rx,emul(mx,rz)),mz))-ryt;
    UpperIntervalVectorType dyt=prod(mSinv,rr);
    UpperIntervalVectorType dz=-rz-feasibility_trmul(mA,dyt);
    UpperIntervalVectorType dx=-ediv(UpperIntervalVectorType(rx+emul(mx,dz)),mz);
    LOGGING_PRINTLN_AT(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz<<"\n");

    UpperIntervalVectorType nx,ny,nyt,nz; FloatDP nt;
    nx=mx+dx;
    nyt=myt+dyt;
    nz=mz+dz;

    x=intersection(x,nx);
    yt=intersection(yt,nyt);
    z=intersection(z,nz);

    y=project(yt,range(0,m));
    t=yt[m];
}

*/


} // namespace Ariadne
