/***************************************************************************
 *            solving/nonlinear_programming_internal.hpp
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
#ifndef ARIADNE_NONLINEAR_PROGRAMMING_INTERNAL_HPP
#define ARIADNE_NONLINEAR_PROGRAMMING_INTERNAL_HPP

#include "function/functional.hpp"
#include <limits>

#include "macro/general.hpp"
#include "logging/logging.hpp"
#include "utility/tuple.hpp"
#include "foundation/logical.hpp"
#include "numeric/numeric.hpp"
#include "algebra/linear_algebra.decl.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"
#include "algebra/differential.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "function/procedure.hpp"

#include "solving/nonlinear_programming.hpp"
#include "solving/nonlinear_interior_point_method.hpp"
#include "solving/nonlinear_penalty_function_method.hpp"
#include "algebra/multi_index-noaliasing.hpp"
#include "solving/constraint_solver.hpp"

#include "algebra/expansion.inl.hpp"

namespace Ariadne {

inline Sweeper<FloatDP> default_sweeper() { return Sweeper<FloatDP>(); }

typedef DiagonalMatrix<FloatDPBounds> FloatDPBoundsDiagonalMatrix;

typedef Vector<FloatDPApproximation> FloatDPApproximationVector;
typedef VectorRange<FloatDPApproximationVector> FloatDPApproximationVectorRange;
typedef Matrix<FloatDPApproximation> FloatDPApproximationMatrix;
typedef DiagonalMatrix<FloatDPApproximation> FloatDPApproximationDiagonalMatrix;
typedef Differential<FloatDPApproximation> FloatDPApproximationDifferential;

typedef Vector<FloatDPBounds> FloatDPBoundsVector;
typedef Matrix<FloatDPBounds> FloatDPBoundsMatrix;
typedef Differential<FloatDPBounds> FloatDPBoundsDifferential;
typedef Vector<FloatDP> ExactFloatDPVectorType;

typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Matrix<UpperIntervalType> UpperIntervalMatrixType;

typedef FloatDPApproximation ApproximateNumericType;

inline Vector<FloatDP> exact_vector(Vector<FloatDPApproximation> const& v) {
    Vector<FloatDP> result(v.size(),dp);
    for(SizeType i=0; i!=v.size(); ++i) {
        result[i]=Ariadne::cast_exact(v[i]);
    }
    return result;
}

Matrix<ApproximateNumericType> join(Matrix<ApproximateNumericType> const&, Matrix<ApproximateNumericType> const&, Matrix<ApproximateNumericType> const&);

inline Vector<Differential<RawFloatDP>>const& cast_raw(Vector<Differential<FloatDPApproximation>>const& v) {
    return reinterpret_cast<Vector<Differential<RawFloatDP>>const&>(v);
}

inline Vector<FloatDPApproximation>& cast_approximate(Vector<RawFloatDP>& v) {
    return reinterpret_cast<Vector<FloatDPApproximation>&>(v);
}

inline Vector<FloatDPApproximation>& cast_approximate(Vector<FloatDPApproximation>& v) {
    return v;
}
inline Vector<Differential<FloatDPApproximation>>const& cast_approximate(Vector<Differential<RawFloatDP>>const& v) {
    return reinterpret_cast<Vector<Differential<FloatDPApproximation>>const&>(v);
}

inline UpperIntervalType dot(Vector<UpperIntervalType> const& bx1, Vector<ExactIntervalType> const& bx2) {
    return dot(bx1,Vector<UpperIntervalType>(bx2));
}

template<class X1, class IVL2> inline decltype(auto) dot(Vector<X1> const& x1, Box<IVL2> const& bx2) {
    return dot(x1,cast_vector(bx2)); }

template<class X> inline
DiagonalMatrix<X> diagonal_matrix(const Vector<X>& v) {
    return DiagonalMatrix<X>(v.array());
}

template<class X> inline
Bool epos(const Vector<X>& x) {
    for(SizeType i=0; i!=x.size(); ++i) { if(x[i]<=0) { return false; } } return true;
}

template<class X> inline
Bool eneg(const Vector<X>& x) {
    for(SizeType i=0; i!=x.size(); ++i) { if(x[i]>=0) { return false; } } return true;
}

template<class X, class XX> inline
Bool egtr(const Vector<X>& x, const XX& s) {
    for(SizeType i=0; i!=x.size(); ++i) { if(decide(x[i]<=s)) { return false; } } return true;
}

template<class X, class XX> inline
Bool elss(const Vector<X>& x, const XX& s) {
    for(SizeType i=0; i!=x.size(); ++i) { if(decide(x[i]>=s)) { return false; } } return true;
}

template<class X> inline
Vector<X> eadd(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]+y[i]; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const X& s) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]-s; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]-y[i]; } return r;
}

template<class X> inline
Vector<X> emul(const Vector<X>& x, const Vector<X>& z) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]*z[i]; } return r;
}

template<class X, class XX> inline
Vector<X> ediv(const Vector<X>& x, const Vector<XX>& z) {
    Vector<X> r(x.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=x[i]/z[i]; } return r;
}

template<class X> inline
Vector<X> ediv(const X& s, const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=s/z[i]; } return r;
}

template<class X> inline
Vector<X> erec(const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=rec(z[i]); } return r;
}

template<class X> inline
Vector<X> esqr(const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(SizeType i=0; i!=r.size(); ++i) { r[i]=sqr(z[i]); } return r;
}

inline
ExactIntervalType eivl(const RawFloatDPVector& x) {
    ARIADNE_ASSERT(x.size()>0); ExactIntervalType r=ExactIntervalType(FloatDP(x[0]));
    for(SizeType i=1; i!=x.size(); ++i) { r=hull(r,FloatDP(x[i])); } return r;
}

inline Matrix<ApproximateNumericType> join(Matrix<ApproximateNumericType> const& A1, Matrix<ApproximateNumericType> const& A2, Matrix<ApproximateNumericType> const& A3) {
    SizeType m=A1.row_size(); SizeType n1=A1.column_size(); SizeType n2=A2.column_size(); SizeType n3=A3.column_size();
    Matrix<ApproximateNumericType> A123(m,n1+n2+n3,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
    project(A123,range(0,m),range(0,n1))=A1;
    project(A123,range(0,m),range(n1,n1+n2))=A2;
    project(A123,range(0,m),range(n1+n2,n1+n2+n3))=A3;
    return A123;
}

template<class X> Matrix<X> cojoin(Matrix<X> const& A1, Matrix<X> const& A2, Matrix<X> const& A3) {
    SizeType n=A1.column_size(); SizeType m1=A1.row_size(); SizeType m2=A2.row_size(); SizeType m3=A3.row_size();
    Matrix<X> A123(m1+m2+m3,n,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
    project(A123,range(0,m1),range(0,n))=A1;
    project(A123,range(m1,m1+m2),range(0,n))=A2;
    project(A123,range(m1+m2,m1+m2+m3),range(0,n))=A3;
    return A123;
}


// Compute S+=ADA^T, where D is diagonal and S is symmetric.
template<class X>
Void adat(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            X ADij=A[i1][j]*D[j];
            for(SizeType i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(SizeType i1=1; i1!=m; ++i1) {
        for(SizeType i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
Void atda(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    ARIADNE_PRECONDITION(S.row_size()==S.column_size());
    ARIADNE_PRECONDITION(S.column_size()==A.column_size());
    ARIADNE_PRECONDITION(D.size()==A.row_size());

    const SizeType m=A.column_size();
    const SizeType n=A.row_size();
    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            X ATDij=A[j][i1]*D[j];
            for(SizeType i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ATDij*A[j][i2];
            }
        }
    }
    for(SizeType i1=1; i1!=m; ++i1) {
        for(SizeType i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
Void atda(Matrix<X>& S, const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    atda(S,A,D.diagonal());
}

// Compute S=ADA^T, where D is diagonal.
template<class X>
Matrix<X> adat(const Matrix<X>& A, const Vector<X>& D)
{
    const SizeType m=A.row_size();
    Matrix<X> S=Matrix<X>::zero(m,m);
    adat(S,A,D);
    return S;
}

// Compute S+=AA^T
template<class X>
Matrix<X> amulat(const Matrix<X>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    Matrix<X> S(m,m);
    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            for(SizeType i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=A[i1][j]*A[i2][j];
            }
        }
    }
    for(SizeType i1=1; i1!=m; ++i1) {
        for(SizeType i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
    return S;
}

template<class X> inline Bool all_greater(const Vector<X>& x, const X& e) {
    for(SizeType i=0; i!=x.size(); ++i) { if(x[i]<=e) { return false; } } return true;
}






template<class X> Vector< Differential<X> > second_derivative(const ValidatedVectorMultivariateFunction& f, const Vector<X>& x) {
    Vector< Differential<X> > d=Differential<X>::variables(f.result_size(),f.argument_size(),2,x);
    return f.evaluate(d);
}

template<class Vec, class Diff> Void set_gradient(Vec& g, const Diff& D) {
    SizeType i=0;
    typename Diff::ConstIterator iter=D.begin();
    if(iter!=D.end() && iter->index().degree()==0) { ++iter; }
    while(iter!=D.end() && iter->index().degree()<=2) {
        while(iter->index()[i]==0) { ++i; }
        g[i]=iter->coefficient();
        ++iter;
    }
}

template<class Mx, class Diff> Void set_jacobian_transpose(Mx& A, const Vector<Diff>& D) {
    for(SizeType j=0; j!=A.column_size(); ++j) {
        for(SizeType i=0; i!=A.row_size(); ++i) {
            A[i][j]=D[j][i];
        }
    }
}

template<class Mx, class Diff> Void set_hessian(Mx& H, const Diff& D) {
    typedef typename Diff::ValueType X;
    SizeType i=0; SizeType j=1;
    typename Diff::ConstIterator iter=D.begin();
    while(iter!=D.end() && iter->index().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->index().degree()<=2) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        while(a[i]==0) { ++i; j=i+1; }
        if(a[i]==2) { H[i][i]=c; }
        else { while(a[j]==0) { ++j; } H[i][j]=c; H[j][i]=c; }
        ++iter;
    }
}

template<class Mx, class S, class Diff> Void add_hessian(Mx& H, const S& s, const Diff& D) {
    typedef typename Diff::ValueType X;
    typename Diff::ConstIterator iter=D.begin();
    while(iter!=D.end() && iter->index().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->index().degree()==2) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        SizeType i=0;
        while(a[i]==0) { ++i; }
        if(a[i]==2) { H[i][i]+=s*c; }
        else { SizeType j=i+1; while(a[j]==0) { ++j; } H[i][j]+=s*c; H[j][i]+=s*c; }
        ++iter;
    }
}

// Compute the product (A -A I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_mul(const Matrix<XX>& A, const Vector<X>& v)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(v.size()==2*(m+n));
    Vector<X> r(m+1u,v.zero_element());
    for(SizeType i=0; i!=m; ++i) {
        r[i]=v[2*n+i]-v[2*n+m+i];
        for(SizeType j=0; j!=n; ++j) {
            r[i]+=A[i][j]*(v[j]-v[n+j]);
        }
    }
    for(SizeType k=0; k!=2*(m+n); ++k) {
        r[m]+=v[k];
    }
    return r;
}

// Compute the product (AT 1 \\ -AT 1 \\ I 1 \\ -I 1) I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_trmul(const Matrix<XX>& A, const Vector<X>& w)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(w.size()==m+1);
    Vector<X> r(2*(m+n),w.zero_element());
    for(SizeType j=0; j!=n; ++j) {
        r[j]=0;
        for(SizeType i=0; i!=m; ++i) {
            r[j]+=A[i][j]*w[i];
        }
        r[n+j]=-r[j];
        r[j]+=w[m];
        r[n+j]+=w[m];
    }
    for(SizeType i=0; i!=m; ++i) {
        r[2*n+i]=w[i]+w[m];
        r[2*n+m+i]=-w[i]+w[m];
    }
    return r;
}


// Compute the product \f$\hat{A}^T \hat{D} \hat{A} + \hat{H}\f$ where \f$\hat{A}=\left(\begin{matrix}A&-A&I&-I\\1&1&1&1\end{matrix}\right)\f$ and \f$\hat{D}=D\f$ is diagonal.
template<class X> Matrix<X> feasibility_adat(const Matrix<X>& H, const Matrix<X>& A, const Vector<X>& D)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    ARIADNE_ASSERT(H.row_size()==m);
    ARIADNE_ASSERT(H.column_size()==m);
    ARIADNE_ASSERT(D.size()==2*(m+n));
    Matrix<X> S(m+1,m+1,H.zero_element());

    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=m; ++j) { S[i][j] = H[i][j]; } }
    for(SizeType i=0; i!=m; ++i) { S[i][m]=0; S[m][i]=0; } S[m][m]=0;

    for(SizeType i1=0; i1!=m; ++i1) {
        for(SizeType j=0; j!=n; ++j) {
            X ADij=A[i1][j]*(D[j]+D[n+j]);
            for(SizeType i2=0; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(SizeType i=0; i!=m; ++i) {
        S[i][i]+=(D[2*n+i]+D[2*n+m+i]);
    }
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            S[i][m]+=A[i][j]*(D[j]-D[n+j]);
        }
        S[i][m]+=(D[2*n+i]-D[2*n+m+i]);
        S[m][i]=S[i][m];
    }
    for(SizeType k=0; k!=2*(m+n); ++k) {
        S[m][m]+=D[k];
    }

    return S;
}






template<class R>
class ConstrainedFeasibilityMatrix {
    ConstrainedFeasibilityMatrix(const Vector<R>& x, const Vector<R>& z, const Matrix<R>& a, const Matrix<R>& h)
        : X(x), Z(z), D(ediv(x,z)), A(a), H(h) { }

    template<class RR> Tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    mul(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Vector<RR> nx=Z*x+X*z;
        Vector<RR> nyt=H*yt-A*x;
        Vector<RR> nz=H*yt-A*x;
        return make_tuple(nx,nyt,nz);
    }

    template<class RR> Tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    solve(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Sinv=inverse(feasibility_adat(H,A,D));
        Vector<RR> rx=Z.solve(x);
        Vector<RR> ryt=yt+feasibility_mul(A,rx-D*z);
        Vector<RR> rz=z;
        ryt=Sinv*ryt;
        rz=rz-feasibility_trmul(A,ryt);
        rx=rz-D*rz;
        return make_tuple(rx,ryt,rz);
    }

    DiagonalMatrix<R> X;
    DiagonalMatrix<R> Z;
    DiagonalMatrix<R> D;
    const Matrix<R>& A;
    const Matrix<R>& H;
    Matrix<R> Sinv;
};


inline ExactBoxType cast_exact_widen(ExactBoxType const& bx, RawFloatDP e) {
    ExactBoxType r(bx);
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=ExactIntervalType(sub(down,bx[i].lower_bound(),e),add(up,bx[i].upper_bound(),e));
    }
    return bx;
}



typedef NonlinearProgrammingSolverBase::ValidatedVectorType ValidatedVectorType;

} // namespace Ariadne

#endif
