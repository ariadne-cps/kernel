/***************************************************************************
 *            solving/equation_solver.cpp
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

#include "solving/equation_solver.hpp"

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

PositiveValidatedUpperNumber sup_error(ValidatedVectorMultivariateFunctionPatch const& fp) {
    return fp.error();
}

} // namespace




static const Bool ALLOW_PARTIAL_FUNCTION = true;

FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory();


EquationSolverBase::EquationSolverBase(ApproximateDouble max_error, CounterType max_steps)
  : _max_error(cast_exact(max_error)), _max_steps(max_steps), _function_factory_ptr(make_taylor_function_patch_factory())
{
}

Void
EquationSolverBase::set_function_factory(const FunctionPatchFactoryInterface<ValidatedTag>& factory)
{
    this->_function_factory_ptr=std::shared_ptr< FunctionPatchFactoryInterface<ValidatedTag> >(factory.clone());
}

const FunctionPatchFactoryInterface<ValidatedTag>&
EquationSolverBase::function_factory() const
{
    return *this->_function_factory_ptr;
}


auto EquationSolverBase::solve(const ValidatedVectorMultivariateFunction& f,
                       const Vector<ValidatedNumericType>& ix) const
    -> Vector<ValidatedNumericType>
{
    ExactBoxType bx=cast_exact_box(ix);
    return this->solve(f,bx);
}

auto EquationSolverBase::solve(const ValidatedVectorMultivariateFunction& f,
                       const ExactBoxType& bx) const
    -> Vector<ValidatedNumericType>
{
    Set< Vector<ValidatedNumericType> > r = this->solve_all(f,bx);
    if(r.size()==0u) { ARIADNE_THROW(NoSolutionException,"EquationSolverBase::solve","no solution in solve("<<f<<","<<bx<<")"); }
    if(r.size()!=1u) { ARIADNE_THROW(EquationSolverException,"EquationSolverBase::solve","non-unique solution in solve("<<f<<","<<bx<<")"); }
    return *r.begin();
}

Void
solve_all(Set< Vector<EquationSolverInterface::ValidatedNumericType> >& r,
          const EquationSolverInterface& s,
          const ValidatedVectorMultivariateFunction& f,
          const ExactBoxType& ix);

auto EquationSolverBase::solve_all(const ValidatedVectorMultivariateFunction& f,
                           const ExactBoxType& bx) const
    -> Set< Vector<ValidatedNumericType> >
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("f="<<f<<", ix="<<bx);

    // Create result set
    Set< Vector<ValidatedNumericType> > r;

    Vector<FloatDPBounds> x=cast_singleton(bx);

    // Test for no solution
    const Vector<ValidatedNumericType> z(bx.size(),dp);
    if(inconsistent(f.evaluate(x),z)) {
        return r;
    }

    Bool invertible_jacobian=true;
    //Vector<ValidatedNumericType> nx=2*ix-cast_exact(bx);
    Vector<ValidatedNumericType> nx=x;
    try {
        Matrix<ValidatedNumericType> Jinv=inverse(f.jacobian(nx));
    }
    catch(const SingularMatrixException&) {
        invertible_jacobian=false;
    }

    Bool need_to_split=true;

    if(invertible_jacobian) {
        //std::cerr<<"Nonsingular matrix -- applying contractor\n";
        try {
            Vector<ValidatedNumericType> y=this->zero(f,bx);
            Bool is_new=true;
            for(Set<Vector<ValidatedNumericType> >::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
                if(consistent(y,*iter)) {
                    is_new=false;
                    break;
                }
            }
            if(is_new) {
                r.insert(y);
            }
            need_to_split=false;
        }
        catch(const NoSolutionException&) {
            //ARIADNE_WARN("NoSolutionException exception: "<<e.what());
            need_to_split=false;
        }
        catch(const EquationSolverException& e) {
            ARIADNE_WARN("EquationSolverException exception: "<<e.what());
            need_to_split=true;
            // No solution found, try splitting
        }
    } else {
        need_to_split=true;
    }

    if(need_to_split) {
        // If sup_error is too small, assume solution is not verified
        if(definitely(sup_error(cast_singleton(bx))<this->maximum_error())) {
            if(!invertible_jacobian) {
                ARIADNE_WARN("Cannot verify solution in "<<bx<<" with f="<<f(cast_singleton(bx))<<"); "
                             <<"Jacobian "<<f.jacobian(nx)<<" is not invertible; "
                             <<"approximate inverse="<<inverse(midpoint(f.jacobian(nx))));
            } else {
                ARIADNE_WARN("Cannot verify or falsify solution in "<<bx<<"; f("<<bx<<")="<<f(cast_singleton(bx))<<".");
            }
            return r;
        }

        //std::cerr<<"  Splitting "<<bx<<"\n";
        Pair< Box<ExactIntervalType>, Box<ExactIntervalType> > splt=split(bx);
        r.adjoin(this->solve_all(f,splt.first));
        r.adjoin(this->solve_all(f,splt.second));
    }

    return r;
}





auto EquationSolverBase::zero(const ValidatedVectorMultivariateFunction& f,
                      const ExactBoxType& bx) const
    -> Vector<ValidatedNumericType>
{
    LOGGING_SCOPE_CREATE;

    const ExactDouble e=this->maximum_error();
    CounterType n=this->maximum_number_of_steps();
    Vector<ValidatedNumericType> r=cast_singleton(bx);
    Vector<ValidatedNumericType> nr(r.size(),dp);
    Bool has_solution=false;
    while(n>0) {
        nr=this->step(f,r);
        LOGGING_PRINTLN_AT(1,"nr="<<nr);

        if(!has_solution && refines(nr,r)) {
            has_solution=true;
        }

        if(has_solution && definitely(sup_error(nr) < e)) {
            return nr;
        }

        if(!consistent(nr,r)) {
            ARIADNE_THROW(NoSolutionException,"EquationSolverBase::zero","No result found in "<<bx<<"; "<<nr<<" is inconsistent with "<<r);
        }
        r=refinement(nr,r);
        n=n-1;
    }
    if(!consistent(f.evaluate(r),Vector<ValidatedNumericType>(f.result_size(),dp))) {
        ARIADNE_THROW(NoSolutionException,"EquationSolverBase::zero","No result found in "<<bx<<"; f("<<r<<") is inconsistent with zero");
    } else {
        FloatDPError widen=FloatDPError(FloatDP::eps(dp))*sup_error(r);
        r+=Vector<ValidatedNumericType>(r.size(),ValidatedNumericType(-widen,+widen));
        nr=this->step(f,r);
        if(refines(nr,r)) {
            return nr;
        }
        ARIADNE_THROW(EquationSolverException,"EquationSolverBase::zero","No result verified in "<<bx<<"; maximum number of steps reached with approximation "<<r<<" which cannot be robustly checked");
    }
}



auto EquationSolverBase::fixed_point(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const
    -> Vector<ValidatedNumericType>
{
    ValidatedVectorMultivariateFunction id=ValidatedVectorMultivariateFunction::identity(f.argument_size());
    return this->solve(f-id,bx);
}


auto EquationSolverBase::implicit(const ValidatedVectorMultivariateFunction& f,
                          const ExactBoxType& ip,
                          const ExactBoxType& ix) const
    -> ValidatedVectorMultivariateFunctionPatch
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1,"f="<<f);
    ARIADNE_ASSERT(f.result_size()==ix.size());
    ARIADNE_ASSERT(f.argument_size()==ip.size()+ix.size());

    const SizeType n=ix.size();
    const ExactDouble err=this->maximum_error();

    ValidatedVectorMultivariateFunctionPatch id(this->function_factory().create_identity(ip));
    ValidatedVectorMultivariateFunctionPatch h(this->function_factory().create_constants(ip,cast_singleton(ix)));
    ValidatedVectorMultivariateFunctionPatch nh(this->function_factory().create_zeros(n,ip));
    ValidatedVectorMultivariateFunctionPatch fnh(this->function_factory().create_zeros(n,ip));

    CounterType steps_remaining=this->maximum_number_of_steps();
    SizeType number_unrefined=n;
    Array<Bool> is_refinement(n,false);

    while(steps_remaining>0) {
        LOGGING_PRINTLN_AT(1,"step="<<this->maximum_number_of_steps()-steps_remaining);
        nh=this->implicit_step(f,id,h);
        fnh=compose(f,join(id,nh));
        LOGGING_PRINTLN_AT(2,"nh="<<nh);
        LOGGING_PRINTLN_AT(2,"fnh="<<fnh);

        if(ALLOW_PARTIAL_FUNCTION) {
            for(SizeType i=0; i!=n; ++i) {
                if(!is_refinement[i]) {
                    LOGGING_PRINTLN_AT(3,"refines(nh["<<i<<"],x["<<i<<"])="<<refines(nh[i],h[i]));
                    if(refines(nh[i],h[i])) {
                        is_refinement[i]=true;
                        --number_unrefined;
                    } else if(definitely(disjoint(ValidatedScalarMultivariateFunctionPatch(nh[i]).range(),ix[i]))) {
                        ARIADNE_THROW(NoSolutionException,"EquationSolverBase::implicit","No result found in "<<ix<<"; function "<<nh<<" is disjoint from "<<h<<" for at least one point.");
                    // TODO: Avoid cast to Dyadic below
                    } else if (cast_exact(nh[i].error())>static_cast<Dyadic>(cast_exact(mag(h[i].range())))) {
                       ARIADNE_THROW(UnknownSolutionException,"EquationSolverBase::implicit","No convergence looking for implicit function with domain "<<h.domain()<<" and codomain "<<ix<<"; "<<nh<<" error bound is larger than range of "<<h);
                    }
                }
            }

            LOGGING_PRINTLN_AT(2,"is_refinement="<<is_refinement);
            LOGGING_PRINTLN_AT(2,"nh.range()="<<nh.range());
            LOGGING_PRINTLN_AT(2,"sup_error(nh)="<<sup_error(nh));
            LOGGING_PRINTLN_AT(2,"sup_error(fnh)="<<sup_error(fnh));

            h=nh;

        } else { // !ALLOW_PARTIAL_FUNCTION
            for(SizeType i=0; i!=n; ++i) {
                if(!is_refinement[i]) {
                    if(refines(nh[i],h[i])) {
                        is_refinement[i]=true;
                        --number_unrefined;
                    } else if(inconsistent(nh[i],h[i])) {
                        ARIADNE_THROW(NoSolutionException,"EquationSolverBase::implicit","No result found in "<<ix<<"; "<<nh<<" is disjoint from "<<h);
                    // TODO: Avoid cast to Dyadic below
                    } else if (cast_exact(nh[i].error())>static_cast<Dyadic>(cast_exact(mag(h[i].range())))) {
                        ARIADNE_THROW(UnknownSolutionException,"EquationSolverBase::implicit","No convergence looking for implicit function with domain "<<h.domain()<<" and codomain "<<ix<<"; "<<nh<<" error bound is larger than range of "<<h);
                    }
                }
            }

            h=refinement(nh,h);
        }

        steps_remaining=steps_remaining-1;
        if( (number_unrefined==0) && ( (steps_remaining==0) || (definitely(sup_error(nh)<err) && definitely(sup_error(fnh)<err)) ) ) {
            return ValidatedVectorMultivariateFunctionPatch(h);
        }
    }

    ARIADNE_THROW(NoSolutionException,"EquationSolverBase::implicit","Could not prove existence of a solution in "<<ix<<".");
}

ValidatedScalarMultivariateFunctionPatch
EquationSolverBase::implicit(const ValidatedScalarMultivariateFunction& f,
                     const ExactBoxType& ip,
                     const ExactIntervalType& ix) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN_AT(1,"f="<<f);
    ValidatedVectorMultivariateFunctionPatch res=this->implicit(ValidatedVectorMultivariateFunction(List<ValidatedScalarMultivariateFunction>(1u,f)),ip,ExactBoxType(1u,ix));
    return res[0];
}

auto EquationSolverBase::continuation(const ValidatedVectorMultivariateFunction&,
                              const Vector<ApproximateNumericType>&,
                              const ExactBoxType&,
                              const ExactBoxType&) const
    -> ValidatedVectorMultivariateFunctionModelType
{
    ARIADNE_NOT_IMPLEMENTED;
}





} // namespace Ariadne
