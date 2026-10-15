/***************************************************************************
 *            solving/interval_newton_method.cpp
 *
 *  Copyright  2006-20  Pieter Collins
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

#include "function/functional.hpp"
#include "interval/interval.hpp"
#include "function/function_model.hpp"
#include "solving/interval_newton_method.hpp"
#include "logging/logging.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "function/function_model.hpp"

namespace Ariadne {

namespace {

template<class PR> auto
operator*(const Matrix<FloatBounds<PR>>& A,const ValidatedVectorMultivariateFunctionModel<PR>& v)
    -> ValidatedVectorMultivariateFunctionModel<PR>
{
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorMultivariateFunctionModel<PR> r(A.row_size(),factory(v).create_zero());
    for(SizeType i=0; i!=r.size(); ++i) {
        ValidatedScalarMultivariateFunctionModel<PR> t=r[i];
        for(SizeType j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

template<class X> auto
operator*(const Matrix<X>& A,const ValidatedVectorMultivariateFunctionPatch& v)
    -> ValidatedVectorMultivariateFunctionPatch
{
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorMultivariateFunctionPatch r(A.row_size(),factory(v).create_zero());
    for(SizeType i=0; i!=r.size(); ++i) {
        ValidatedScalarMultivariateFunctionPatch t=r[i];
        for(SizeType j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

template<class X> auto
linear_solve(const Matrix<X>& A,const ValidatedVectorMultivariateFunctionPatch& v)
    -> ValidatedVectorMultivariateFunctionPatch
{
    return inverse(A)*v;
}

template<class PR> FloatError<PR> sup_error(const ValidatedVectorMultivariateFunctionModel<PR>& x) {
    ARIADNE_PRECONDITION(x.size()>0);
    FloatError<PR> r(0u,x[0].error().precision());
    for(SizeType i=0; i!=x.size(); ++i) { r=max(r,x[i].error()); }
    return r;
}

PositiveValidatedUpperNumber sup_error(ValidatedVectorMultivariateFunctionPatch const& fp) {
    return fp.error();
}

} // namespace

auto IntervalNewtonSolver::step(const ValidatedVectorMultivariateFunction& f,
                                        const Vector<ValidatedNumericType>& x) const
    -> Vector<ValidatedNumericType>
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("Testing for root in "<<x);
    LOGGING_PRINTLN_AT(1,"e="<<sup_error(x)<<", x="<<x);
    Vector<FloatDP> m(cast_exact(x));
    LOGGING_PRINTLN_AT(1,"m="<<m);
    Vector<ValidatedNumericType> im(m);
    Vector<ValidatedNumericType> w=f.evaluate(im);
    LOGGING_PRINTLN_AT(1,"f(m)="<<w);
    Matrix<ValidatedNumericType> A=f.jacobian(x);
    LOGGING_PRINTLN_AT(1,"Df(r)="<<A);
    Matrix<ValidatedNumericType> Ainv=inverse(A);
    LOGGING_PRINTLN_AT(1,"inverse(Df(r))="<<Ainv);
    Vector<ValidatedNumericType> dx=Ainv*w;
    LOGGING_PRINTLN_AT(1,"dx="<<dx);
    Vector<ValidatedNumericType> nx= m - dx;
    LOGGING_PRINTLN_AT(1,"nx="<<nx);
    return nx;
}

ValidatedVectorMultivariateFunctionPatch
IntervalNewtonSolver::implicit_step(const ValidatedVectorMultivariateFunction& f,
                                    const ValidatedVectorMultivariateFunctionPatch& id,
                                    const ValidatedVectorMultivariateFunctionPatch& h) const
{
    LOGGING_SCOPE_CREATE;
    const SizeType m=id.result_size();
    const SizeType n=h.result_size();
    DP pr;

    LOGGING_PRINTLN("f="<<f);
    LOGGING_PRINTLN("h="<<h);
    ValidatedVectorMultivariateFunctionPatch mh=h;
    mh.clobber();
    LOGGING_PRINTLN("midpoint(h)="<<mh);

    ValidatedScalarMultivariateFunction zero_function(f.domain());
    Matrix<ValidatedScalarMultivariateFunction> D2f(n,n,zero_function);
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            D2f[i][j]=f[i].derivative(m+j);
        }
    }
    LOGGING_PRINTLN("D2f="<<D2f);

    ValidatedNumericType zero(0,pr);
    ValidatedScalarMultivariateFunctionPatch z=h[0]*zero;
    ValidatedVectorMultivariateFunctionPatch idh=join(id,h);

    Matrix<ValidatedScalarMultivariateFunctionPatch> J(n,n,z);
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            J[i][j]=compose(D2f[i][j],idh);
        }
    }
    LOGGING_PRINTLN("J="<<J);

    Matrix<UpperIntervalType> rngJ(n,n,pr);
    for(SizeType i=0; i!=n; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            UpperIntervalType D2fij=UpperIntervalType(unchecked_evaluate(D2f[i][j],cast_singleton(product(id.range(),h.range()))));
            rngJ[i][j]=intersection(J[i][j].range(),D2fij);
        }
    }
    LOGGING_PRINTLN("rngJ="<<rngJ);

    ValidatedVectorMultivariateFunctionPatch fidmh=compose(f,join(id,mh));
    LOGGING_PRINTLN_AT(1,"compose(f,join(id,midpoint(h)))="<<fidmh);

    ValidatedVectorMultivariateFunctionPatch dh(n,z);
    if(n==1) {
        if(possibly(contains(rngJ[0][0],FloatDP(0.0_x,dp)))) {
            ARIADNE_THROW(SingularJacobianException,"IntervalNewtonSolver","D2f(P,X)="<<rngJ[0][0]<<" which contains zero.");
        }
        if(possibly(contains(J[0][0].range(),FloatDP(0.0_x,dp)))) {
            dh[0]=fidmh[0]/cast_singleton(rngJ[0][0]);
        } else {
            dh[0]=fidmh[0]/J[0][0];
        }
    } else {
        dh=Ariadne::linear_solve(cast_singleton(rngJ),fidmh);
    }
    LOGGING_PRINTLN("dh="<<dh);

    ValidatedVectorMultivariateFunctionPatch nh=mh-dh;
    LOGGING_PRINTLN("nh="<<nh);
    return nh;
}


Void IntervalNewtonSolver::_write(OutputStream& os) const
{
    os << "IntervalNewtonSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}


} // namespace Ariadne
