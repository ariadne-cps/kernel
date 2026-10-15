/***************************************************************************
 *            solving/krawczyk_method.cpp
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
#include "solving/krawczyk_method.hpp"
#include "logging/logging.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/differential.hpp"
#include "algebra/algebra.hpp"
#include "function/taylor_model.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "function/function_model.hpp"

#include "algebra/expansion.inl.hpp"
#include "algebra/evaluate.hpp"
#include "algebra/evaluate.tpl.hpp"

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

} // namespace

auto KrawczykSolver::step(const ValidatedVectorMultivariateFunction& f,
                          const Vector<ValidatedNumericType>& x) const
    -> Vector<ValidatedNumericType>
{
    LOGGING_SCOPE_CREATE
    Matrix<ValidatedNumericType> I=Matrix<ValidatedNumericType>::identity(x.size(),x.zero_element());
    LOGGING_PRINTLN("Testing for root in "<<x);
    LOGGING_PRINTLN_AT(1,"e="<<sup_error(x)<<", x="<<x);
    Vector<FloatDP> m(cast_exact(x));
    LOGGING_PRINTLN_AT(1,"m="<<m);
    Vector<ValidatedNumericType> im(m);
    Vector<ValidatedNumericType> fm=f.evaluate(im);
    LOGGING_PRINTLN_AT(1,"f(m)="<<fm);
    Matrix<ValidatedNumericType> J=f.jacobian(x);
    LOGGING_PRINTLN_AT(1,"Df(r)="<<J);
    Matrix<ValidatedNumericType> M=inverse(midpoint(J));
    LOGGING_PRINTLN_AT(1,"inverse(Df(m))="<<M);
    Vector<ValidatedNumericType> dx=M*fm-(I-M*J)*(x-m);
    LOGGING_PRINTLN_AT(1,"dx="<<dx);
    Vector<ValidatedNumericType> nx= m - dx;
    LOGGING_PRINTLN_AT(1,"nx="<<nx);
    Vector<ValidatedNumericType> nr(nx);
    LOGGING_PRINTLN_AT(1,"nr="<<nr);
    return nr;
}


auto FactoredKrawczykSolver::step(const ValidatedVectorMultivariateFunction& f,
                                  const Vector<ValidatedNumericType>& x) const
    -> Vector<ValidatedNumericType>
{
    LOGGING_SCOPE_CREATE;
    Matrix<ValidatedNumericType> I=Matrix<ValidatedNumericType>::identity(x.size(),x.zero_element());
    LOGGING_PRINTLN("Testing for root in "<<x);
    LOGGING_PRINTLN_AT(1,"e="<<sup_error(x)<<", x="<<x);
    Vector<FloatDP> m(cast_exact(x));
    LOGGING_PRINTLN_AT(1,"m="<<m);
    Vector<ValidatedNumericType> im(m);
    Vector<ValidatedNumericType> fm=f.evaluate(im);
    LOGGING_PRINTLN_AT(1,"f(m)="<<fm);
    Matrix<ValidatedNumericType> J=f.jacobian(x);
    LOGGING_PRINTLN_AT(1,"Df(r)="<<J);
    Matrix<ValidatedNumericType> mJ(midpoint(J));
    Matrix<ValidatedNumericType> M=inverse(mJ);
    LOGGING_PRINTLN_AT(1,"inverse(Df(m))="<<M);
    Vector<ValidatedNumericType> dx=M*(fm+(J-mJ)*(x-m));
    LOGGING_PRINTLN_AT(1,"dx="<<dx);
    Vector<ValidatedNumericType> nx= m - dx;
    LOGGING_PRINTLN_AT(1,"nx="<<nx);
    Vector<ValidatedNumericType> nr(nx);
    LOGGING_PRINTLN_AT(1,"nr="<<nr);
    return nr;
}


ValidatedVectorMultivariateFunctionPatch
KrawczykSolver::implicit_step(const ValidatedVectorMultivariateFunction& f,
                              const ValidatedVectorMultivariateFunctionPatch& p,
                              const ValidatedVectorMultivariateFunctionPatch& x) const
{
    LOGGING_SCOPE_CREATE;
    const SizeType np=p.size();
    const SizeType nx=x.size();
    Matrix<ValidatedNumericType> I=Matrix<ValidatedNumericType>::identity(nx,dp);
    LOGGING_PRINTLN("Contracting x="<<x);
    LOGGING_PRINTLN("p="<<p);
    LOGGING_PRINTLN("f="<<f);
    ValidatedVectorMultivariateFunctionPatch mx(x);
    for(SizeType i=0; i!=mx.size(); ++i) { mx[i].clobber(); }
    LOGGING_PRINTLN("mx="<<mx);
    Vector<FloatDPError> ex(nx,dp);
    for(SizeType i=0; i!=nx; ++i) { ex[i]=x[i].error(); }
    Vector<ValidatedNumericType> eix=make_bounds(ex);
    LOGGING_PRINTLN("ex="<<ex);
    ValidatedVectorMultivariateFunctionPatch fm=compose(f,join(p,mx));
    LOGGING_PRINTLN_AT(1,"f(p,mx)="<<fm);
    Vector<ValidatedNumericType> rp(np,dp);
    for(SizeType i=0; i!=np; ++i) { rp[i]=cast_singleton(p[i].range()); }
    Vector<ValidatedNumericType> rx(nx,dp);
    for(SizeType i=0; i!=nx; ++i) { rx[i]=cast_singleton(x[i].range()); }
    Matrix<ValidatedNumericType> J=project(f.jacobian(join(rp,rx)),range(0,nx),range(np,np+nx));
    LOGGING_PRINTLN("D2f(r)=J="<<J);
    Matrix<ValidatedNumericType> M=inverse(midpoint(J));
    LOGGING_PRINTLN("inverse(D2f(m))=M="<<M);
    LOGGING_PRINTLN("M*f(p,mx)="<<M*fm);
    LOGGING_PRINTLN("(I-M*J)="<<(I-M*J));
    LOGGING_PRINTLN("(I-M*J) * (ex*ValidatedNumericType(-1,+1))="<<(I-M*J)<<"*"<<eix<<"="<<(I-M*J) * eix);
    ValidatedVectorMultivariateFunctionPatch dx= M*fm - (I-M*J) * eix;
    LOGGING_PRINTLN("dx="<<dx);
    ValidatedVectorMultivariateFunctionPatch nwx= mx - dx;
    LOGGING_PRINTLN("nwx="<<nwx);
    return nwx;
}


Void KrawczykSolver::_write(OutputStream& os) const
{
    os << "KrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

Void FactoredKrawczykSolver::_write(OutputStream& os) const
{
    os << "FactoredKrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}



} // namespace Ariadne
