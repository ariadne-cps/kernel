/***************************************************************************
 *            function_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "pybind11.hpp"

#include "kernel-utilities.hpp"
#include "numeric_submodule.hpp"

#include <iostream>
#include <iomanip>

#include "utility/array.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/differential.hpp"
#include "function/formula.hpp"
#include "algebra/polynomial.hpp"
#include "function/affine.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/function_expression.hpp"
#include "solving/constraint_solver.hpp"

#include "function/function_mixin.tpl.hpp"

namespace Ariadne {

template<> struct PythonTemplateName<Procedure> { static std::string get() { return "Procedure"; } };
template<> struct PythonTemplateName<Function> { static std::string get() { return "Function"; } };

template<class P, class ARG> pybind11::object function_from_python(typename DomainTraits<ARG>::EntireDomainType dom, pybind11::object pyf) {
    Function<P,ARG(ARG)> id = Function<P,ARG(ARG)>::identity(dom);
    pybind11::object pyobjf = pyf(id);
    try {
        VectorFunction<P,ARG> vf=pybind11::cast<VectorFunction<P,ARG>>(pyobjf);
        return pybind11::cast(vf);
    } catch (...) { }
    try {
        pybind11::list pylstf = pybind11::cast<pybind11::list>(pyobjf);
        SizeType rs = len(pylstf);
        //VectorFunction<P,ARG>(Vector<ScalarFunction<P,ARG>>(rs, [&id,&lstf](SizeType i){return pybind11::cast<ScalarFunction<P,ARG>>(lstf[i](id));}));
        List<ScalarFunction<P,ARG>> lsf;
        for (SizeType i=0; i!=rs; ++i) {
            lsf.append(pybind11::cast<ScalarFunction<P,ARG>>(pylstf[i]));
        }
        return pybind11::cast(VectorFunction<P,ARG>(lsf));
    } catch (...) { }
    return pybind11::cast(pybind11::cast<ScalarFunction<P,ARG>>(pyobjf));
}

pybind11::object univariate_function_from_python(pybind11::object pyf) {
    return function_from_python<EffectiveTag,RealScalar>(RealDomain(),pyf);
}
pybind11::object multivariate_function_from_python(SizeType as, pybind11::object pyf) {
    return function_from_python<EffectiveTag,RealVector>(EuclideanDomain(as),pyf);
}



template<class X, class D>
inline Matrix<X> get_jacobian(const Vector<D>& d) {
    const Nat rs=d.size(); const Nat as=d[0].argument_size();
    Matrix<X> J(rs,as);
    for(Nat i=0; i!=rs; ++i) {
        for(Nat j=0; j!=as; ++j) {
            J[i][j]=d[i][j];
        }
    }
    return J;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<MultiIndex>& arepr) {
    MultiIndex const& a=arepr.reference(); os << "("; for(SizeType i=0; i!=a.size(); ++i) { if(i!=0) { os << ','; } os << a[i]; } os << ")"; return os;
}

template<class P, class SIG> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Function<P,SIG> >& frepr) {
    Function<P,SIG> const& f = frepr.reference();
    os << "Function(";
    if constexpr (Same<decltype(f.argument_size()),SizeOne>) { }
    else { os << frepr.reference().argument_size() << ","; }
    return os << " lambda x: " << frepr.reference() << " )";
}

}// namespace Ariadne

using namespace Ariadne;

typedef Vector<FloatDPApproximation> FV;
typedef Vector<ExactIntervalType> IV;
typedef Matrix<FloatDPApproximation> FMx;
typedef Matrix<ExactIntervalType> IMx;
typedef Vector< Differential<FloatDPApproximation> > FSDV;
typedef Vector< Differential<ExactIntervalType> > ISDV;
typedef Vector<ValidatedTaylorModelDP> TMV;
typedef ValidatedVectorMultivariateTaylorFunctionModelDP TFM;
typedef ValidatedTaylorModelDP TM;


template<class AT, class F> void define_call(pybind11::class_<F>& f) {
    typedef typename F::template Argument<AT> A; typedef ResultOf<F(A)> R; f.def("__call__", (R(F::*)(A const&)const) &F::operator()); }

template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<ScalarUnivariateFunction<P>>& function_class) {
    typedef ScalarUnivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(slope(declval<F>(),declval<A>())) R;
    module.def("slope", (R(*)(F const&,A const&)) &slope);
    function_class.def("slope", (R(F::*)(A const&)const) &F::slope);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<ScalarMultivariateFunction<P>>& function_class) {
    typedef ScalarMultivariateFunction<P> F; typedef typename F::template Argument<AT> A;
    typedef decltype(gradient(declval<F>(),declval<A>())) R;
    module.def("gradient", (R(*)(F const&,A const&)) &gradient);
    function_class.def("gradient", (R(F::*)(A const&)const) &F::gradient);
    typedef decltype(hessian(declval<F>(),declval<A>())) H;
    module.def("hessian", (H(*)(F const&,A const&)) &hessian);
    //function_class.def("hessian", (H(F::*)(A const&)const) &F::hessian);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<VectorUnivariateFunction<P>>& function_class) {
    typedef VectorUnivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(tangent(declval<F>(),declval<A>())) R;
    module.def("tangent", (R(*)(F const&,A const&)) &tangent);
    function_class.def("tangent", (R(F::*)(A const&)const) &F::tangent);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<VectorMultivariateFunction<P>>& function_class) {
    typedef VectorMultivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(jacobian(declval<F>(),declval<A>())) R;
    module.def("jacobian", (R(*)(F const&,A const&)) &jacobian);
    function_class.def("jacobian", (R(F::*)(A const&)const) &F::jacobian);
}


template<class F, class T> using ArgumentType = typename F::template Argument<T>;

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, ApproximateTag) {
    define_call<FloatDPApproximation>(function_class);
    define_call<FloatMPApproximation>(function_class);
    define_call<Differential<FloatDPApproximation>>(function_class);
    define_call<Differential<FloatMPApproximation>>(function_class);

    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatDPApproximation>>);
    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatMPApproximation>>);
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, ValidatedTag) {
    define_call<FloatDPBounds>(function_class);
    define_call<FloatMPBounds>(function_class);
    define_call<Differential<FloatDPBounds>>(function_class);
    define_call<Differential<FloatMPBounds>>(function_class);

    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatDPBounds>>);
    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatMPBounds>>);

    export_function_evaluation(module,function_class,ApproximateTag());
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, EffectiveTag) {
    export_function_evaluation(module,function_class,ValidatedTag());
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class)
{
    export_function_evaluation(module, function_class, Paradigm<F>());
}


template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class, ApproximateTag) {
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatDPApproximation>,DegreeType>);
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatMPApproximation>,DegreeType>);

    module.def("differential", _differential_<F,ArgumentType<F,FloatDPApproximation>,DegreeType>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatMPApproximation>,DegreeType>);

    define_named_derivatives<FloatDPApproximation>(module,function_class);
    define_named_derivatives<FloatMPApproximation>(module,function_class);
}

template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class, ValidatedTag) {
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatDPBounds>,DegreeType>);
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatMPBounds>,DegreeType>);

    module.def("differential", _differential_<F,ArgumentType<F,FloatDPBounds>,DegreeType>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatMPBounds>,DegreeType>);

    define_named_derivatives<FloatDPBounds>(module,function_class);
    define_named_derivatives<FloatMPBounds>(module,function_class);

    export_function_derivatives(module,function_class,ApproximateTag());
}

template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class, EffectiveTag) {
    export_function_derivatives(module,function_class,ValidatedTag());
}

template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class)
{
    export_function_derivatives(module, function_class, Paradigm<F>());
}


Void export_domains(pybind11::module& module)
{
    pybind11::class_<RealDomain> real_domain_class(module,"RealDomain");
    real_domain_class.def(pybind11::init<>());
    //    real_domain_class.def(pybind11::init<SizeOne>());
    real_domain_class.def("__repr__", &__cstr__<RealDomain>);

    pybind11::class_<EuclideanDomain> euclidean_domain_class(module,"EuclideanDomain");
    euclidean_domain_class.def(pybind11::init<SizeType>());
    euclidean_domain_class.def("__repr__", &__cstr__<EuclideanDomain>);
}

template<class SIG> String signature_name() {
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    return String(Same<RES,RealScalar> ? "Scalar" : "Vector") + (Same<ARG,RealScalar> ? "Univariate" : "Multivariate");
}

template<class P, class SIG> Void export_function(pybind11::module& module) {
    using F=Function<P,SIG>;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using DomainType = typename SignatureTraits<SIG>::DomainType;
    using CodomainType = typename SignatureTraits<SIG>::CodomainType;
    using ResultSizeType = ElementSizeType<CodomainType>;
    using ArgumentSizeType = ElementSizeType<DomainType>;

    pybind11::class_<F> function_class(module,(class_name<P>()+signature_name<SIG>()+"Function").c_str());

    function_class.def(pybind11::init<F>());
    if constexpr (Same<P,ValidatedTag>) {
        function_class.def(pybind11::init<Function<EffectiveTag,SIG>>());
        pybind11::implicitly_convertible<Function<EffectiveTag,SIG>,Function<ValidatedTag,SIG>>();
    }

    function_class.def(pybind11::init<ResultSizeType,ArgumentSizeType>());
    function_class.def(pybind11::init<ResultSizeType,DomainType>());
    if constexpr (Same<RES,RealScalar>) {
        function_class.def(pybind11::init<ArgumentSizeType>());
        function_class.def(pybind11::init<DomainType>());
    }

    function_class.def(pybind11::init<typename SignatureTraits<SIG>::ArgumentSpaceType, typename SignatureTraits<SIG>::template Result<RealExpression>>());

    if constexpr (Same<RES,RealVector>) {
        // NOTE: This must go *after* the conversion constructor
        function_class.def(pybind11::init([](std::vector<Function<P,RealScalar(ARG)>> const& lst){return F(lst);}));
    }

    function_class.def("result_size", &F::result_size);
    function_class.def("argument_size", &F::argument_size);

    if constexpr (Same<RES,RealVector>) {
        function_class.def("__getitem__", &F::get);
        function_class.def("__setitem__", &F::set);
    }

    if constexpr (Same<ARG,RealScalar>) {
        function_class.def("derivative", (F(F::*)(IndexZero)const) &F::derivative);
        function_class.def("derivative", (F(*)(const F&)) &derivative);
        module.def("derivative", (F(F::*)(IndexZero)const) &F::derivative);
        module.def("derivative", (F(*)(const F&)) &derivative);
    } else if constexpr (Same<ARG,RealVector>) {
        function_class.def("derivative", (F(F::*)(SizeType)const) &F::derivative);
        module.def("derivative", (F(F::*)(SizeType)const) &F::derivative);
    }


    if constexpr (Same<RES,RealScalar>) { define_elementary_algebra<F,Number<P>>(module,function_class); }

    // TODO: Put these in C++ API
    // define_vector_algebra_arithmetic<VectorMultivariateFunction<P>,ScalarMultivariateFunction<P>,Number<P>>(module,vector_function_class);
    // FIXME: Define vector function operations for Validated and Approximate
    if constexpr (Same<RES,RealVector>) {
        if constexpr (Same<ARG,RealVector> and Same<P,EffectiveTag>) {
            define_vector_algebra_arithmetic<VectorFunction<P,ARG>,ScalarFunction<P,ARG>>(module,function_class);
        }
    }

    export_function_evaluation(module,function_class);
    export_function_derivatives(module,function_class);

    if constexpr (Same<SIG,RealScalar(RealScalar)>) {
        if constexpr (not Same<P,ApproximateTag>) {
            module.def("derivative", [](const F& f, const FloatDPBounds& x){return static_cast<FloatDPBounds>(differential(f,x,1u).gradient()[0]);} );
            module.def("derivative", [](const F& f, const FloatMPBounds& x){return static_cast<FloatMPBounds>(differential(f,x,1u).gradient()[0]);} );
            module.def("second_derivative", [](const F& f, const FloatDPBounds& x){return static_cast<FloatDPBounds>(differential(f,x,2u).hessian()[0][0]);} );
            module.def("second_derivative", [](const F& f, const FloatMPBounds& x){return static_cast<FloatMPBounds>(differential(f,x,2u).hessian()[0][0]);} );
        }
        module.def("derivative", [](const F& f, const FloatDPApproximation& x){return static_cast<FloatDPApproximation>(differential(f,x,1u).gradient()[0]);} );
        module.def("derivative", [](const F& f, const FloatMPApproximation& x){return static_cast<FloatMPApproximation>(differential(f,x,1u).gradient()[0]);} );
        module.def("second_derivative", [](const F& f, const FloatDPApproximation& x){return static_cast<FloatDPApproximation>(differential(f,x,2u).hessian()[0][0]);} );
        module.def("second_derivative", [](const F& f, const FloatMPApproximation& x){return static_cast<FloatMPApproximation>(differential(f,x,2u).hessian()[0][0]);} );
    }

    if constexpr (Same<SIG,RealScalar(RealScalar)>) {
        function_class.def_static("constant", (F(*)(Number<P>)) &F::constant);
        function_class.def_static("constant", (F(*)(RealDomain,Number<P>)) &F::constant);
        function_class.def_static("coordinate", (F(*)()) &F::coordinate);
        function_class.def_static("coordinate", (F(*)(RealDomain)) &F::coordinate);
        function_class.def_static("identity", (F(*)()) &F::identity);
        function_class.def_static("identity", (F(*)(RealDomain)) &F::identity);
    }
    if constexpr (Same<SIG,RealVector(RealScalar)>) {
        function_class.def_static("zeros", (F(*)(SizeType)) &F::zeros);
        function_class.def_static("zeros", (F(*)(SizeType,RealDomain)) &F::zeros);
    }
    if constexpr (Same<SIG,RealScalar(RealVector)>) {
        function_class.def_static("constant", (F(*)(SizeType,Number<P>)) &F::constant);
        function_class.def_static("constant", (F(*)(EuclideanDomain,Number<P>)) &F::constant);
        function_class.def_static("coordinate", (F(*)(SizeType,SizeType)) &F::coordinate);
        function_class.def_static("coordinate", (F(*)(EuclideanDomain,SizeType)) &F::coordinate);
    }
    if constexpr (Same<SIG,RealVector(RealVector)>) {
        function_class.def_static("zeros", (F(*)(SizeType,SizeType)) &F::zeros);
        function_class.def_static("zeros", (F(*)(SizeType,EuclideanDomain)) &F::zeros);
        function_class.def_static("identity", (F(*)(SizeType)) &F::identity);
        function_class.def_static("identity", (F(*)(EuclideanDomain)) &F::identity);
    }

    function_class.def("__str__", &__cstr__<F>);
    function_class.def("__repr__", &__repr__<F>);
    function_class.def("__crepr__", &__crepr__<F>);

    if constexpr (Same<RES,RealVector>) if constexpr (Same<ARG,RealVector>) {
        module.def("join", (VectorFunction<P,ARG>(*)(const ScalarFunction<P,ARG>&, const ScalarFunction<P,ARG>&)) &join);
        module.def("join", (VectorFunction<P,ARG>(*)(const VectorFunction<P,ARG>&, const ScalarFunction<P,ARG>&)) &join);
        module.def("join", (VectorFunction<P,ARG>(*)(const ScalarFunction<P,ARG>&, const VectorFunction<P,ARG>&)) &join);
        module.def("join", (VectorFunction<P,ARG>(*)(const VectorFunction<P,ARG>&, const VectorFunction<P,ARG>&)) &join);
    }

    module.def("compose", (Function<P,RES(ARG)>(*)(const Function<P,RES(RealScalar)>&,const Function<P,RealScalar(ARG)>&)) &compose);
    module.def("compose", (Function<P,RES(ARG)>(*)(const Function<P,RES(RealVector)>&,const Function<P,RealVector(ARG)>&)) &compose);

    if constexpr (Same<P,EffectiveTag> and Same<SIG,RealScalar(RealVector)>) {
        function_class.def("__eq__", &__eq__<EffectiveScalarMultivariateFunction,EffectiveNumber,Return<EffectiveConstraint>>);
        function_class.def("__le__", &__le__<EffectiveScalarMultivariateFunction,EffectiveNumber,Return<EffectiveConstraint>>);
        function_class.def("__ge__", &__ge__<EffectiveScalarMultivariateFunction,EffectiveNumber,Return<EffectiveConstraint>>);
    }

}


template<class SIG> Void export_functions(pybind11::module& module)
{
    export_function<EffectiveTag,SIG>(module);
    export_function<ValidatedTag,SIG>(module);
    export_function<ApproximateTag,SIG>(module);

    pybind11::implicitly_convertible<Function<EffectiveTag,SIG>,Function<ValidatedTag,SIG>>();
    pybind11::implicitly_convertible<Function<EffectiveTag,SIG>,Function<ApproximateTag,SIG>>();
    pybind11::implicitly_convertible<Function<ValidatedTag,SIG>,Function<ApproximateTag,SIG>>();

    if constexpr (Same<SIG,RealScalar(RealVector)>) {
        module.def("derivatives", (VectorMultivariateFunction<EffectiveTag>(*)(const ScalarMultivariateFunction<EffectiveTag>&)) &derivatives);
        module.def("lie_derivative", (ScalarMultivariateFunction<EffectiveTag>(*)(const ScalarMultivariateFunction<EffectiveTag>&,const VectorMultivariateFunction<EffectiveTag>&)) &lie_derivative);
    }
}


template<class Y, class X> Void define_procedure_evaluation(pybind11::module& module) {
    module.def("evaluate", (X(*)(Procedure<Y> const&, Vector<X> const&)) &evaluate);
    module.def("gradient", (Covector<X>(*)(Procedure<Y> const&, Vector<X> const&)) &gradient);
    module.def("hessian", (X(*)(Procedure<Y> const&, Vector<X> const&, Vector<X> const&)) &hessian);
}

template<class Y> pybind11::class_<Procedure<Y>> export_procedure(pybind11::module& module) {
    typedef Paradigm<Y> P;
    pybind11::class_<Procedure<Y>> procedure_class(module,(class_name<P>()+"Procedure").c_str());
    procedure_class.def("__str__", &__cstr__<Procedure<Y>>);
    module.def("make_procedure", (Procedure<Y>(*)(ScalarMultivariateFunction<P> const&)) &make_procedure);
    if constexpr (Same<Y,ApproximateNumber>) {
        define_procedure_evaluation<Y,FloatDPApproximation>(module);
        //define_procedure_evaluation<Y,FloatMPApproximation>(module);
    }
    if constexpr (Same<Y,ValidatedNumber>) {
        define_procedure_evaluation<Y,FloatDPBounds>(module);
        //define_procedure_evaluation<Y,FloatMPBounds>(module);
    }
    return procedure_class;
}

Void export_procedures(pybind11::module& module) {
    export_procedure<ApproximateNumber>(module);
    export_procedure<ValidatedNumber>(module);
}



template<class P, class SIG> Void export_function_patch(pybind11::module& module) {
    using FP=FunctionPatch<P,SIG>;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using DomainType = typename FunctionPatch<P,SIG>::DomainType;
    using ArgumentIndexType = ElementIndexType<DomainType>;

    pybind11::class_<FP> function_patch_class(module,(class_name<P>()+signature_name<SIG>()+"FunctionPatch").c_str());
    function_patch_class.def(pybind11::init<FunctionPatch<P,SIG>>());

    function_patch_class.def("result_size", &FP::result_size);
    function_patch_class.def("argument_size", &FP::argument_size);
    function_patch_class.def("domain", &FP::domain);
    function_patch_class.def("codomain", &FP::codomain);
    function_patch_class.def("range", &FP::range);

    function_patch_class.def("__str__", &__cstr__<FP>);

    if constexpr (Same<RES,RealScalar>) {
        define_elementary_algebra<FP,Number<P>>(module,function_patch_class);
    }
    if constexpr (Same<RES,RealVector>) {
        define_vector_algebra_arithmetic<FunctionPatch<P,RealVector(ARG)>,FunctionPatch<P,RealScalar(ARG)>,Number<P>>(module,function_patch_class);
    }

    // TODO: Export all function operators
    // export_function_evaluation(module,function_patch_class);
    // export_function_derivatives(module,function_patch_class);


    if constexpr (IsStronger<P,ApproximateTag>::value) {
        define_call<FloatDPApproximation>(function_patch_class);
        define_call<FloatMPApproximation>(function_patch_class);
    }
    if constexpr (IsStronger<P,ValidatedTag>::value) {
        define_call<FloatDPBounds>(function_patch_class);
        define_call<FloatMPBounds>(function_patch_class);
    }

    if  constexpr (Same<RES,RealVector>) {
        module.def("antiderivative", [](FunctionPatch<P,SIG>const& fp, ArgumentIndexType j){return _antiderivative_(fp,j);});
        module.def("antiderivative", [](FunctionPatch<P,SIG>const& fp, ArgumentIndexType j, Number<P> const& c){return _antiderivative_(fp,j,c);});
    }

    module.def("compose", [](const Function<P,RealScalar(RES)>& f, const FunctionPatch<P,SIG>& gp){return compose(f,gp);});
    module.def("compose", [](const Function<P,RealVector(RES)>& f, const FunctionPatch<P,SIG>& gp){return compose(f,gp);});
    module.def("compose", [](const FunctionPatch<P,RealScalar(RES)>& fp, const FunctionPatch<P,SIG>& gp){return compose(fp,gp);});
    module.def("compose", [](const FunctionPatch<P,RealVector(RES)>& fp, const FunctionPatch<P,SIG>& gp){return compose(fp,gp);});
    module.def("unchecked_compose", [](const FunctionPatch<P,RealScalar(RES)>& fp, const FunctionPatch<P,SIG>& gp){return unchecked_compose(fp,gp);});
    module.def("unchecked_compose", [](const FunctionPatch<P,RealVector(RES)>& fp, const FunctionPatch<P,SIG>& gp){return unchecked_compose(fp,gp);});

    //module.def("restriction", (FunctionPatch<P,SIG>(*)(Function<P,SIG> const&, DomainType const&)) &restriction);
}

template<class SIG> Void export_function_patches(pybind11::module& module) {
    export_function_patch<ValidatedTag,SIG>(module);
    export_function_patch<ApproximateTag,SIG>(module);
}

Void export_function_patches(pybind11::module& module) {
    // TODO: Correct export of univariate function patches.
    export_function_patches<RealScalar(RealScalar)>(module);
    // export_function_patches<RealVector(RealScalar)>(module);
    export_function_patches<RealScalar(RealVector)>(module);
    export_function_patches<RealVector(RealVector)>(module);
}


Void export_function_models(pybind11::module& module);

Void function_submodule(pybind11::module& module) {



    export_domains(module);

    export_functions<RealScalar(RealScalar)>(module);
    export_functions<RealVector(RealScalar)>(module);
    export_functions<RealScalar(RealVector)>(module);
    export_functions<RealVector(RealVector)>(module);

    export_function_patches(module);
    export_procedures(module);

    template_<Function> function_template(module,python_template_name<Function>().c_str());
    function_template.def_new([](RealVariable var, Scalar<RealExpression> e){return EffectiveScalarUnivariateFunction(var,e);});
    function_template.def_new([](RealVariable var, Vector<RealExpression> e){return EffectiveVectorUnivariateFunction(var,e);});
    function_template.def_new([](RealSpace spc, Scalar<RealExpression> e){return EffectiveScalarMultivariateFunction(spc,e);});
    function_template.def_new([](RealSpace spc, Vector<RealExpression> e){return EffectiveVectorMultivariateFunction(spc,e);});

    function_template.def_new([](pybind11::object pyf){return univariate_function_from_python(pyf);});
    function_template.def_new([](SizeType as, pybind11::object pyf){return multivariate_function_from_python(as,pyf);});

    export_function_models(module);
}

/***************************************************************************
 *            calculus_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include <type_traits>

#include "pybind11.hpp"
#include "kernel-utilities.hpp"
#include "numeric_submodule.hpp"

#include "algebra/expansion.tpl.hpp"
#include "algebra/algebra.hpp"
#include "function/function_interface.hpp"
#include "algebra/polynomial.hpp"
#include "function/function.hpp"
#include "function/procedure.hpp"

#include "function/function_model.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"

using namespace Ariadne;

namespace Ariadne {

template<class PR> std::string numeric_class_tag();

template<> struct PythonTemplateName<ApproximateTaylorModel> { static std::string get() { return "ApproximateTaylorModel"; } };
template<> struct PythonTemplateName<ValidatedTaylorModel> { static std::string get() { return "ValidatedTaylorModel"; } };
template<> struct PythonTemplateName<ValidatedScalarMultivariateFunctionModel> { static std::string get() { return "ValidatedScalarMultivariateFunctionModel"; } };
template<> struct PythonTemplateName<ValidatedVectorMultivariateFunctionModel> { static std::string get() { return "ValidatedVectorMultivariateFunctionModel"; } };
template<> struct PythonTemplateName<ValidatedScalarMultivariateTaylorFunctionModel> { static std::string get() { return "ValidatedScalarMultivariateTaylorFunctionModel"; } };
template<> struct PythonTemplateName<ValidatedVectorMultivariateTaylorFunctionModel> { static std::string get() { return "ValidatedVectorMultivariateTaylorFunctionModel"; } };

template<template<class...>class T, class PR> std::string python_template_tag_name() { return python_template_name<T>()+numeric_class_tag<PR>(); }
template<template<class...>class T, class FLT> std::string python_template_precision_tag_name() { return python_template_name<T>()+numeric_class_tag<typename FLT::PrecisionType>(); }

template<class PR> struct PythonClassName<ValidatedScalarMultivariateFunctionModel<PR>> { static std::string get() { return python_template_tag_name<ValidatedScalarMultivariateFunctionModel,PR>(); } };
template<class PR> struct PythonClassName<ValidatedVectorMultivariateFunctionModel<PR>> { static std::string get() { return python_template_tag_name<ValidatedVectorMultivariateFunctionModel,PR>(); } };


template<class FLT> struct PythonClassName<ApproximateTaylorModel<FLT>> {
    static std::string get() { return python_template_precision_tag_name<ApproximateTaylorModel,FLT>(); } };
template<class FLT> struct PythonClassName<ValidatedTaylorModel<FLT>> {
    static std::string get() { return python_template_precision_tag_name<ValidatedTaylorModel,FLT>(); } };

template<class FLT> struct PythonClassName<ValidatedScalarMultivariateTaylorFunctionModel<FLT>> {
    static std::string get() { return python_template_precision_tag_name<ValidatedScalarMultivariateTaylorFunctionModel,FLT>(); } };
template<class FLT> struct PythonClassName<ValidatedVectorMultivariateTaylorFunctionModel<FLT>> {
    static std::string get() { return python_template_precision_tag_name<ValidatedVectorMultivariateTaylorFunctionModel,FLT>(); } };


//ValidatedNumericType evaluate(const ValidatedScalarMultivariateFunctionModelDP& f, const Vector<ValidatedNumericType>& x) { return f(x); }
//Vector<ValidatedNumericType> evaluate(const ValidatedVectorMultivariateFunctionModelDP& f, const Vector<ValidatedNumericType>& x) { return f(x); }

//ValidatedScalarMultivariateFunctionModelDP partial_evaluate(const ValidatedScalarMultivariateFunctionModelDP&, SizeType, const ValidatedNumericType&);
//ValidatedVectorMultivariateFunctionModelDP partial_evaluate(const ValidatedVectorMultivariateFunctionModelDP&, SizeType, const ValidatedNumericType&);

//ValidatedScalarMultivariateFunctionModelDP antiderivative(const ValidatedScalarMultivariateFunctionModelDP&, SizeType, ValidatedNumericType);
//ValidatedVectorMultivariateFunctionModelDP antiderivative(const ValidatedVectorMultivariateFunctionModelDP&, SizeType, ValidatedNumericType);

//ValidatedScalarMultivariateFunctionModelDP compose(const ValidatedScalarMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedScalarMultivariateFunctionModelDP compose(const ValidatedScalarMultivariateFunction&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP compose(const ValidatedVectorMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP compose(const ValidatedVectorMultivariateFunction&, const ValidatedVectorMultivariateFunctionModelDP&);

//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedScalarMultivariateFunctionModelDP&, const ValidatedScalarMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedScalarMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedVectorMultivariateFunctionModelDP&, const ValidatedScalarMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedVectorMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);

template<class P, class SIG, class PR, class PRE> OutputStream& operator<<(OutputStream& os, const Representation< FunctionModel<P,SIG,PR,PRE> >& frepr) {
    static_cast<const FunctionModelInterface<P,SIG,PR,PRE>&>(frepr.reference())._repr(os); return os;
}

ValidatedVectorMultivariateTaylorFunctionModelDP __getslice__(const ValidatedVectorMultivariateTaylorFunctionModelDP& tf, Int start, Int stop) {
    Int rs = static_cast<Int>(tf.result_size());
    if(start<0) { start+=rs; }
    if(stop<0) { stop+=rs; }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&stop<=rs,
            "result_size="<<rs<<", start="<<start<<", stop="<<stop);
    return ValidatedVectorMultivariateTaylorFunctionModelDP(tf.domain(),Vector<ValidatedTaylorModelDP>(project(tf.models(),range(static_cast<SizeType>(start),static_cast<SizeType>(stop)))));
}






template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Vector<X> >& repr) {
    const Vector<X>& vec=repr.reference();
    os << "[";
    for(SizeType i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << python_representation(vec[i]);
    }
    os << "]";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation< ExactBoxType >& bx) {
    return os << PythonRepresentation< Vector<ExactIntervalType> >(cast_vector(bx.reference())); }

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedScalarMultivariateTaylorFunctionModelDP>& repr) {
    const ValidatedScalarMultivariateTaylorFunctionModelDP& stf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedScalarMultivariateTaylorFunctionModelDP"
       << "(" << python_representation(stf.domain())
       << "," << python_representation(stf.expansion())
       << "," << python_representation(stf.error())
       << "," << python_representation(stf.properties())
       << ")";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedVectorMultivariateTaylorFunctionModelDP>& repr) {
    const ValidatedVectorMultivariateTaylorFunctionModelDP& vtf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedVectorMultivariateTaylorFunctionModelDP"
       << "(" << python_representation(vtf.domain())
       << "," << python_representation(vtf.expansions())
       << "," << python_representation(vtf.errors())
       << "," << python_representation(vtf.properties())
       << ")";
    return os;
}

template<class P, class F> List<MultiIndex> keys(const TaylorModel<P,F>& tm) {
    List<MultiIndex> r(tm.argument_size());
    for(auto iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
        r.append(iter->index());
    }
    return r;
}

ValidatedScalarMultivariateFunction unrestrict(const ValidatedScalarMultivariateFunctionModelDP& fm) {
    return ValidatedScalarMultivariateFunction(fm.raw_pointer()->_clone());
}

ValidatedVectorMultivariateFunction unrestrict(const ValidatedVectorMultivariateFunctionModelDP& fm) {
    return ValidatedVectorMultivariateFunction(fm.raw_pointer()->_clone());
}

} // namespace Ariadne

static constexpr auto py_self = pybind11::detail::self;

Expansion<MultiIndex,FloatDP>const& get_expansion(ValidatedTaylorModelDP const& tm) { return tm.expansion(); }
Expansion<MultiIndex,FloatDPApproximation>const& get_expansion(ApproximateTaylorModelDP const& tm) { return tm.expansion(); }


template<class FLT> Void export_validated_taylor_model(pybind11::module& module)
{
    typedef ValidatedTaylorModel<FLT> ModelType;
    typedef NumericType<ModelType> NumericType;
    typedef GenericType<NumericType> GenericNumericType;

    Tag<GenericNumericType> generic_number_tag;

    pybind11::class_<ValidatedTaylorModel<FLT>> taylor_model_class(module,python_class_name<ValidatedTaylorModel<FLT>>().c_str());
    taylor_model_class.def(pybind11::init<ValidatedTaylorModel<FLT>>());
    taylor_model_class.def(pybind11::init< SizeType,Sweeper<FLT> >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModel<FLT>&))&keys);
    taylor_model_class.def("value", (const FloatDP(ValidatedTaylorModel<FLT>::*)()const) &ValidatedTaylorModel<FLT>::value);
    taylor_model_class.def("error", (const FloatDPError&(ValidatedTaylorModel<FLT>::*)()const) &ValidatedTaylorModel<FLT>::error);
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FLT>&(*)(ValidatedTaylorModel<FLT> const&)) &get_expansion);
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModel<FLT>::*)(const FloatDPError&)) &ValidatedTaylorModel<FLT>::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModel<FLT>::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModel<FLT>::domain);
    taylor_model_class.def("range", &ValidatedTaylorModel<FLT>::range);
    taylor_model_class.def("norm", &ValidatedTaylorModel<FLT>::norm);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModel<FLT>::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModel<FLT>::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModel<FLT>&(ValidatedTaylorModel<FLT>::*)()) &ValidatedTaylorModel<FLT>::sweep, pybind11::return_value_policy::reference);
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModel<FLT>,MultiIndex,FLT>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModel<FLT>,MultiIndex,FLT>);

    define_elementary_algebra(module,taylor_model_class);
    define_inplace_algebra(module,taylor_model_class);
    define_mixed_arithmetic(module, taylor_model_class, generic_number_tag);

    define_lattice(module,taylor_model_class);


    taylor_model_class.def("__str__", &__cstr__<ValidatedTaylorModel<FLT>>);

    taylor_model_class.def_static("constant",(ValidatedTaylorModel<FLT>(*)(SizeType, const FloatDPBounds&,Sweeper<FLT>))&ValidatedTaylorModel<FLT>::constant);
    taylor_model_class.def_static("coordinate",(ValidatedTaylorModel<FLT>(*)(SizeType, SizeType,Sweeper<FLT>))&ValidatedTaylorModel<FLT>::coordinate);

    taylor_model_class.def("range", (UpperIntervalType(ValidatedTaylorModel<FLT>::*)()const) &ValidatedTaylorModel<FLT>::range);

    module.def("evaluate",&_evaluate_<ValidatedTaylorModel<FLT>,Vector<FloatDPBounds>>);
    //module.def("split",&_split_<ValidatedTaylorModel<FLT>,SizeType,SplitPart>);

    //to_python< Vector<ValidatedTaylorModel<FLT>> >();
    //from_python< Vector<ValidatedTaylorModel<FLT>> >();
}

template<class FLT> Void export_approximate_taylor_model(pybind11::module& module)
{
    typedef ApproximateTaylorModel<FLT> ModelType;

    pybind11::class_<ApproximateTaylorModel<FLT>> taylor_model_class(module,python_class_name<ApproximateTaylorModel<FLT>>().c_str());
    taylor_model_class.def(pybind11::init<ModelType>());
    taylor_model_class.def(pybind11::init< SizeType,Sweeper<FLT> >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ModelType&))&keys);
    taylor_model_class.def("value", (const Approximation<FLT>(ModelType::*)()const) &ModelType::value);
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,Approximation<FLT>>&(*)(ModelType const&)) &get_expansion);
    taylor_model_class.def("argument_size", &ModelType::argument_size);
    taylor_model_class.def("domain", &ModelType::domain);
    taylor_model_class.def("range", &ModelType::range);
    taylor_model_class.def("norm", &ModelType::norm);
    taylor_model_class.def("set_sweeper", &ModelType::set_sweeper);
    taylor_model_class.def("sweeper", &ModelType::sweeper);
    taylor_model_class.def("sweep", (ModelType&(ModelType::*)()) &ModelType::sweep, pybind11::return_value_policy::reference);
    taylor_model_class.def("__getitem__", &__getitem__<ModelType,MultiIndex,Approximation<FLT>>);
    taylor_model_class.def("__setitem__",&__setitem__<ModelType,MultiIndex,Approximation<FLT>>);

    define_elementary_algebra(module,taylor_model_class);
    define_inplace_algebra(module,taylor_model_class);
    define_lattice(module,taylor_model_class);

    taylor_model_class.def("__str__",&__cstr__<ModelType>);

    taylor_model_class.def_static("constant",(ModelType(*)(SizeType, const Approximation<FLT>&,Sweeper<FLT>))&ModelType::constant);
    taylor_model_class.def_static("coordinate",(ModelType(*)(SizeType, SizeType,Sweeper<FLT>))&ModelType::coordinate);

    //FIXME: Not in C++ API
    //module.def("evaluate",&_evaluate_<ModelType,Vector<ApproximateNumericType>>);
    //module.def("split",&_split_<ModelType,SizeType,SplitPart>);

    //from_python< Vector<ModelType> >();
    //to_python< Vector<ModelType> >();
}



template<class PR> Void export_scalar_function_model(pybind11::module& module)
{
    using FLT = RawFloatType<PR>;

    typedef typename ValidatedScalarMultivariateFunctionModel<PR>::NumericType NumericType;
    typename NumericType::PrecisionType pr;

    pybind11::class_<ValidatedScalarMultivariateFunctionModel<PR>> scalar_function_model_class(module,python_class_name<ValidatedScalarMultivariateFunctionModel<PR>>().c_str());
    scalar_function_model_class.def(pybind11::init<ValidatedScalarMultivariateFunctionModel<PR>>());
    scalar_function_model_class.def(pybind11::init<ValidatedScalarMultivariateTaylorFunctionModel<FLT>>());
    scalar_function_model_class.def("argument_size", &ValidatedScalarMultivariateFunctionModel<PR>::argument_size);
    scalar_function_model_class.def("domain", &ValidatedScalarMultivariateFunctionModel<PR>::domain);
    scalar_function_model_class.def("codomain", &ValidatedScalarMultivariateFunctionModel<PR>::codomain);
    scalar_function_model_class.def("range", &ValidatedScalarMultivariateFunctionModel<PR>::range);
    scalar_function_model_class.def("clobber", &ValidatedScalarMultivariateFunctionModel<PR>::clobber);
    scalar_function_model_class.def("error", &ValidatedScalarMultivariateFunctionModel<PR>::error);
    scalar_function_model_class.def("__call__", (FloatBounds<PR>(ValidatedScalarMultivariateFunctionModel<PR>::*)(const Vector<FloatBounds<PR>>&)const) &ValidatedScalarMultivariateFunctionModel<PR>::operator());
    scalar_function_model_class.def(py_self+py_self);
    scalar_function_model_class.def(py_self-py_self);
    scalar_function_model_class.def(py_self*py_self);
    scalar_function_model_class.def(py_self/py_self);
    scalar_function_model_class.def(py_self+NumericType(pr));
    scalar_function_model_class.def(py_self-NumericType(pr));
    scalar_function_model_class.def(py_self*NumericType(pr));
    scalar_function_model_class.def(py_self/NumericType(pr));
    scalar_function_model_class.def(NumericType(pr)+py_self);
    scalar_function_model_class.def(NumericType(pr)-py_self);
    scalar_function_model_class.def(NumericType(pr)*py_self);
    scalar_function_model_class.def(NumericType(pr)/py_self);
    scalar_function_model_class.def("__str__", &__cstr__<ValidatedScalarMultivariateFunctionModel<PR>>);
    scalar_function_model_class.def("__repr__", &__crepr__<ValidatedScalarMultivariateFunctionModel<PR>>);

    module.def("norm", &_norm_<ValidatedScalarMultivariateFunctionModel<PR>>);

    module.def("evaluate", &_evaluate_<ValidatedScalarMultivariateFunctionModel<PR>,Vector<NumericType>>);
    module.def("partial_evaluate", &_partial_evaluate_<ValidatedScalarMultivariateFunctionModel<PR>,SizeType,NumericType>);

    module.def("compose", _compose_<ValidatedScalarMultivariateFunctionModel<PR>,ValidatedVectorMultivariateFunctionModel<PR>>);
    module.def("compose", _compose_<ValidatedScalarMultivariateFunction,ValidatedVectorMultivariateFunctionModel<PR>>);

    module.def("unrestrict", (ValidatedScalarMultivariateFunction(*)(const ValidatedScalarMultivariateFunctionModel<PR>&)) &unrestrict);

    module.def("antiderivative",  &_antiderivative_<ValidatedScalarMultivariateFunctionModel<PR>,SizeType,NumericType>);
}

template<class PR> Void export_vector_function_model(pybind11::module& module)
{
    using FLT = RawFloatType<PR>;

    typedef typename ValidatedScalarMultivariateFunctionModel<PR>::NumericType NumericType;
    //using VectorMultivariateFunctionModelType = ValidatedVectorMultivariateFunctionModel<PR>;
    //using ScalarMultivariateFunctionModelType = ValidatedScalarMultivariateFunctionModel<PR>;

    pybind11::class_<ValidatedVectorMultivariateFunctionModel<PR>> vector_function_model_class(module,python_class_name<ValidatedVectorMultivariateFunctionModel<PR>>().c_str());
    vector_function_model_class.def(pybind11::init<ValidatedVectorMultivariateFunctionModel<PR>>());
    vector_function_model_class.def(pybind11::init<ValidatedVectorMultivariateTaylorFunctionModel<FLT>>());
    vector_function_model_class.def("result_size", &ValidatedVectorMultivariateFunctionModel<PR>::result_size);
    vector_function_model_class.def("argument_size", &ValidatedVectorMultivariateFunctionModel<PR>::argument_size);
    vector_function_model_class.def("domain", &ValidatedVectorMultivariateFunctionModel<PR>::domain);
    vector_function_model_class.def("codomain", &ValidatedVectorMultivariateFunctionModel<PR>::codomain);
    vector_function_model_class.def("range", &ValidatedVectorMultivariateFunctionModel<PR>::range);
    //vector_function_model_class.def("__getslice__", (ValidatedVectorMultivariateTaylorFunctionModel<PR>(*)(const ValidatedVectorMultivariateTaylorFunctionModel<PR>&,Int,Int))&__getslice__);
    vector_function_model_class.def("__getitem__", &__getitem__<ValidatedVectorMultivariateFunctionModel<PR>,SizeType,ValidatedScalarMultivariateFunctionModel<PR>>);
    vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorMultivariateFunctionModel<PR>,SizeType,ValidatedScalarMultivariateFunctionModel<PR>>);
    //vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorMultivariateFunctionModel<PR>,SizeType,ValidatedScalarMultivariateFunction>);
    vector_function_model_class.def("__call__", (Vector<FloatBounds<PR>>(ValidatedVectorMultivariateFunctionModel<PR>::*)(const Vector<FloatBounds<PR>>&)const) &ValidatedVectorMultivariateFunctionModel<PR>::operator());

    // NOTE: Not all operations are exported in C++ API.
    //define_vector_algebra_arithmetic<VectorMultivariateFunctionModelType,ScalarMultivariateFunctionModelType,NumericType>(module,vector_function_model_class);
    //define_vector_arithmetic<VectorMultivariateFunctionModelType,ScalarMultivariateFunctionModelType>(module,vector_function_model_class);

    vector_function_model_class.def("__str__", &__cstr__<ValidatedVectorMultivariateFunctionModel<PR>>);
    vector_function_model_class.def("__repr__", &__crepr__<ValidatedVectorMultivariateFunctionModel<PR>>);
    //export_vector_function_model.def(pybind11::module& module, "__repr__",&__repr__<ValidatedVectorMultivariateFunctionModel<PR>>);


//    module.def("evaluate", (Vector<ValidatedNumericType>(*)(const ValidatedVectorMultivariateFunctionModel<PR>&,const Vector<ValidatedNumericType>&)) &evaluate);
    module.def("partial_evaluate", &_partial_evaluate_<ValidatedVectorMultivariateFunctionModel<PR>,SizeType,NumericType>);

    module.def("compose", &_compose_<ValidatedVectorMultivariateFunctionModel<PR>,ValidatedVectorMultivariateFunctionModel<PR>>);
    module.def("compose", &_compose_<ValidatedVectorMultivariateFunction,ValidatedVectorMultivariateFunctionModel<PR>>);

    module.def("unrestrict", (ValidatedVectorMultivariateFunction(*)(const ValidatedVectorMultivariateFunctionModel<PR>&)) &unrestrict);

    module.def("norm", &_norm_<ValidatedVectorMultivariateFunctionModel<PR>>);

    module.def("join", &_join_<ValidatedScalarMultivariateFunctionModel<PR>,ValidatedScalarMultivariateFunctionModel<PR>>);
    module.def("join", &_join_<ValidatedScalarMultivariateFunctionModel<PR>,ValidatedVectorMultivariateFunctionModel<PR>>);
    module.def("join", &_join_<ValidatedVectorMultivariateFunctionModel<PR>,ValidatedScalarMultivariateFunctionModel<PR>>);
    module.def("join", &_join_<ValidatedVectorMultivariateFunctionModel<PR>,ValidatedVectorMultivariateFunctionModel<PR>>);

    module.def("combine", &_combine_<ValidatedVectorMultivariateFunctionModel<PR>,ValidatedVectorMultivariateFunctionModel<PR>>);
        module.def("combine", &_combine_<ValidatedVectorMultivariateFunctionModel<PR>,ValidatedScalarMultivariateFunctionModel<PR>>);

    module.def("antiderivative", &_antiderivative_<ValidatedVectorMultivariateFunctionModel<PR>,SizeType,NumericType>);
    module.def("antiderivative", &_antiderivative_<ValidatedVectorMultivariateFunctionModel<PR>,SizeType,ValidatedNumber>);

//    to_python< List<ValidatedVectorMultivariateFunctionModel<PR>> >();
}


template<class FLT> Void export_scalar_taylor_function(pybind11::module& module)
{
    using FunctionModelType = ValidatedScalarMultivariateTaylorFunctionModel<FLT>;
//    using FunctionType = ScalarMultivariateFunction<Paradigm<FunctionModelType>>;
    using GenericFunctionType = ScalarMultivariateFunction<ValidatedTag>;
    using NumericType = NumericType<FunctionModelType>;
    using GenericNumericType = GenericType<NumericType>;

    typedef ValidatedScalarUnivariateFunction SFN;
    typedef ValidatedVectorUnivariateFunction VFN;
    typedef ValidatedScalarMultivariateTaylorFunctionModel<FLT> F;
    typedef ValidatedVectorMultivariateTaylorFunctionModel<FLT> VF;
    typedef typename F::DomainType D;
    typedef typename F::NumericType X;
    typedef Vector<X> VX;
    typedef SizeType I;

    Tag<GenericNumericType> generic_number_tag;
    Tag<GenericFunctionType> generic_function_tag;

    pybind11::class_<ValidatedScalarMultivariateTaylorFunctionModel<FLT>> scalar_taylor_function_class(module,python_class_name<ValidatedScalarMultivariateTaylorFunctionModel<FLT>>().c_str());
    scalar_taylor_function_class.def(pybind11::init<ValidatedScalarMultivariateTaylorFunctionModel<FLT>>());
    scalar_taylor_function_class.def(pybind11::init<ExactBoxType,ValidatedTaylorModel<FLT>>());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType,Sweeper<FLT> >());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType, const EffectiveScalarMultivariateFunction&,Sweeper<FLT> >());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType, Expansion<MultiIndex,FLT>, Error<FLT>, Sweeper<FLT> >());
    scalar_taylor_function_class.def("error", (const Error<FLT>(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::error);
    scalar_taylor_function_class.def("set_error", (Void(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)(const Error<FLT>&)) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::set_error);
    scalar_taylor_function_class.def("argument_size", &F::argument_size);
    scalar_taylor_function_class.def("domain", &F::domain);
    scalar_taylor_function_class.def("codomain", &F::codomain);
    scalar_taylor_function_class.def("range", &F::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModel<FLT>&(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const)&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::model);
    scalar_taylor_function_class.def("polynomial", (MultivariatePolynomial<NumericType>(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const)&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const)&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::number_of_nonzeros);
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ValidatedScalarMultivariateTaylorFunctionModel<FLT>,MultiIndex,FLT>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ValidatedScalarMultivariateTaylorFunctionModel<FLT>,MultiIndex,FLT>);

    define_elementary_algebra(module,scalar_taylor_function_class);
    define_inplace_algebra(module,scalar_taylor_function_class);
    define_mixed_arithmetic(module,scalar_taylor_function_class,generic_number_tag);
    define_mixed_arithmetic(module,scalar_taylor_function_class,generic_function_tag);
    define_transcendental(module,scalar_taylor_function_class);
    define_lattice(module,scalar_taylor_function_class);

    scalar_taylor_function_class.def("__str__", &__cstr__<F>);
    scalar_taylor_function_class.def("__repr__", &__crepr__<F>);

    scalar_taylor_function_class.def("value", (const FLT(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::value);
    scalar_taylor_function_class.def("clobber", (Void(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::clobber);
    scalar_taylor_function_class.def("set_properties",&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::set_properties);
    scalar_taylor_function_class.def("properties",&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::properties);
    scalar_taylor_function_class.def("__call__", (Approximation<FLT>(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)(const Vector<Approximation<FLT>>&)const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::operator());
    scalar_taylor_function_class.def("__call__", (Bounds<FLT>(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)(const Vector<Bounds<FLT>>&)const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::operator());
//      scalar_taylor_function_class.def("gradient", (Covector<Bounds<FLT>>(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)(const Vector<Bounds<FLT>>&)const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::gradient);
    scalar_taylor_function_class.def("function", (ValidatedScalarMultivariateFunction(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::function);
    scalar_taylor_function_class.def("polynomial", (MultivariatePolynomial<Bounds<FLT>>(ValidatedScalarMultivariateTaylorFunctionModel<FLT>::*)()const) &ValidatedScalarMultivariateTaylorFunctionModel<FLT>::polynomial);
    scalar_taylor_function_class.def("restriction", &_restriction_<F,D>);
    module.def("restrict", &_restriction_<F,D>);
//    scalar_taylor_function_class.def("extension",&_extension_<F,D>);

    scalar_taylor_function_class.def_static("zero",(ValidatedScalarMultivariateTaylorFunctionModel<FLT>(*)(const ExactBoxType&,Sweeper<FLT>))&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::zero);
    scalar_taylor_function_class.def_static("constant",(ValidatedScalarMultivariateTaylorFunctionModel<FLT>(*)(const ExactBoxType&,const NumericType&,Sweeper<FLT>))&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::constant);
    scalar_taylor_function_class.def_static("constant",[](const ExactBoxType& bx, const ValidatedNumber& c,Sweeper<FLT> swp){return ValidatedScalarMultivariateTaylorFunctionModel<FLT>::constant(bx,NumericType(c,swp.precision()),swp);});
    scalar_taylor_function_class.def_static("coordinate",(ValidatedScalarMultivariateTaylorFunctionModel<FLT>(*)(const ExactBoxType&,SizeType,Sweeper<FLT>))&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::coordinate);
    scalar_taylor_function_class.def_static("unit_ball",(ValidatedScalarMultivariateTaylorFunctionModel<FLT>(*)(const ExactBoxType&,Sweeper<FLT>))&ValidatedScalarMultivariateTaylorFunctionModel<FLT>::unit_ball);


    module.def("restriction", &_restriction_<F,D>);
    module.def("join",(VF(*)(F,F)) &_join_<F,F>);
//    module.def("extension",&_extension_<F,D>);
    module.def("embed", &_embed_<D,F,D>);
//    module.def("split", &_split_<F,I>);
    module.def("evaluate", &_evaluate_<F,VX>);
    module.def("partial_evaluate", &_partial_evaluate_<F,I,X>);
    module.def("midpoint", &_midpoint_<F>);
    module.def("derivative", &_derivative_<F,I>);
    module.def("antiderivative", &_antiderivative_<F,I>);
    module.def("antiderivative", &_antiderivative_<F,I,X>);
    module.def("norm", &_norm_<F>);

    module.def("inconsistent", &_inconsistent_<F,F>);
    module.def("refines", &_refines_<F,F>);
    module.def("refinement", &_refinement_<F,F>);

    module.def("compose", &_compose_<SFN,F>);
    module.def("compose", &_compose_<VFN,F>);

//    to_python< Vector<ValidatedScalarMultivariateTaylorFunctionModel<FLT>> >();
}


template<class FLT> Void export_vector_taylor_function(pybind11::module& module)
{
    using PR = typename FLT::PrecisionType;
    using ScalarMultivariateFunctionModelType = ValidatedScalarMultivariateFunctionModel<PR>;
    using NumericType = NumericType<ScalarMultivariateFunctionModelType>;

    typedef SizeType I;
    typedef ValidatedScalarMultivariateFunction SFN;
    typedef ValidatedVectorMultivariateFunction VFN;
    typedef ValidatedScalarMultivariateTaylorFunctionModel<FLT> SF;
    typedef ValidatedVectorMultivariateTaylorFunctionModel<FLT> VF;
    typedef typename VF::DomainType D;
    typedef typename D::IntervalType Di;
    typedef typename VF::NumericType X;
    typedef Vector<X> VX;

    Tag<ValidatedScalarMultivariateFunctionModel<PR>> scalar_taylor_function_tag;
    Tag<Vector<NumericType>> number_vector_tag;

    pybind11::class_<ValidatedVectorMultivariateTaylorFunctionModel<FLT>> vector_taylor_function_class(module,python_class_name<ValidatedVectorMultivariateTaylorFunctionModel<FLT>>().c_str());
    vector_taylor_function_class.def( pybind11::init([](List<ValidatedScalarMultivariateTaylorFunctionModel<FLT>> lst){return ValidatedVectorMultivariateTaylorFunctionModel<FLT>(Vector<ValidatedScalarMultivariateTaylorFunctionModel<FLT>>(lst));}));
    vector_taylor_function_class.def( pybind11::init< SizeType, ExactBoxType, Sweeper<FLT> >());
    vector_taylor_function_class.def( pybind11::init< ExactBoxType,const EffectiveVectorMultivariateFunction&,Sweeper<FLT> >());
    vector_taylor_function_class.def(pybind11::init< ExactBoxType, Vector< Expansion<MultiIndex,FLT> >, Vector<Error<FLT>>, Sweeper<FLT> >());
    vector_taylor_function_class.def( pybind11::init< Vector<ValidatedScalarMultivariateTaylorFunctionModel<FLT>> >());
    vector_taylor_function_class.def("_len_", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::result_size);
    vector_taylor_function_class.def("result_size", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::result_size);
    vector_taylor_function_class.def("argument_size", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::argument_size);
    vector_taylor_function_class.def("domain", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::domain);
    vector_taylor_function_class.def("codomain", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::codomain);
    // FIXME: Omitted since const and non-const versions
    // vector_taylor_function_class.def("models", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::models);
    vector_taylor_function_class.def("centre", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::centre);
    vector_taylor_function_class.def("range", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::range);
    vector_taylor_function_class.def("errors", &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::errors);
    vector_taylor_function_class.def("clobber", (Void(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)()) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::clobber);
    vector_taylor_function_class.def("set_properties",&ValidatedVectorMultivariateTaylorFunctionModel<FLT>::set_properties);
    vector_taylor_function_class.def("properties",&ValidatedVectorMultivariateTaylorFunctionModel<FLT>::properties);
    //FIXME: Omitted since project(...) fails
    //vector_taylor_function_class.def("__getslice__", &__getslice__<ValidatedVectorMultivariateTaylorFunctionModel<FLT>,SizeType,SizeType,ValidatedVectorMultivariateTaylorFunctionModel<FLT>>);
    vector_taylor_function_class.def("__getitem__", &__getitem__<ValidatedVectorMultivariateTaylorFunctionModel<FLT>,SizeType,ValidatedScalarMultivariateTaylorFunctionModel<FLT>>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<ValidatedVectorMultivariateTaylorFunctionModel<FLT>,SizeType,ValidatedScalarMultivariateTaylorFunctionModel<FLT>>);

    define_vector_arithmetic(module,vector_taylor_function_class, scalar_taylor_function_tag);
    //FIXME: Not completely defined in C++ API
    //define_inplace_vector_arithmetic(module,vector_taylor_function_class, scalar_taylor_function_tag);
    define_inplace_mixed_vector_arithmetic(module,vector_taylor_function_class, number_vector_tag);

    vector_taylor_function_class.def("__str__", &__cstr__<ValidatedVectorMultivariateTaylorFunctionModel<FLT>>);
    vector_taylor_function_class.def("__repr__", &__crepr__<ValidatedVectorMultivariateTaylorFunctionModel<FLT>>);
    vector_taylor_function_class.def("clobber", (Void(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)()) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::clobber);
    vector_taylor_function_class.def("__call__", (Vector<Approximation<FLT>>(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)(const Vector<Approximation<FLT>>&)const) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::operator());
    vector_taylor_function_class.def("__call__", (Vector<Bounds<FLT>>(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)(const Vector<Bounds<FLT>>&)const) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::operator());
     vector_taylor_function_class.def("jacobian", (Matrix<Bounds<FLT>>(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)(const Vector<Bounds<FLT>>&)const) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< MultivariatePolynomial<Bounds<FLT>> >(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)()const) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::polynomials);
    vector_taylor_function_class.def("function", (ValidatedVectorMultivariateFunction(ValidatedVectorMultivariateTaylorFunctionModel<FLT>::*)()const) &ValidatedVectorMultivariateTaylorFunctionModel<FLT>::function);

    vector_taylor_function_class.def_static("constant",(ValidatedVectorMultivariateTaylorFunctionModel<FLT>(*)(const ExactBoxType&, const Vector<NumericType>&,Sweeper<FLT>))&ValidatedVectorMultivariateTaylorFunctionModel<FLT>::constant);
    vector_taylor_function_class.def_static("constant",[](const ExactBoxType& bx, const Vector<ValidatedNumber>& c,Sweeper<FLT> swp){return ValidatedVectorMultivariateTaylorFunctionModel<FLT>::constant(bx,Vector<NumericType>(c,swp.precision()),swp);});
    vector_taylor_function_class.def_static("identity",(ValidatedVectorMultivariateTaylorFunctionModel<FLT>(*)(const ExactBoxType&,Sweeper<FLT>))&ValidatedVectorMultivariateTaylorFunctionModel<FLT>::identity);

    module.def("norm", &_norm_<VF>);

    module.def("inconsistent", &_inconsistent_<VF,VF>);
    module.def("refinement", &_refinement_<VF,VF>);
    module.def("refines", &_refines_<VF,VF>);

    module.def("join", &_join_<SF,SF>);
    module.def("join", &_join_<SF,VF>);
    module.def("join", &_join_<VF,SF>);
    module.def("join", &_join_<VF,VF>);
    module.def("combine", &_combine_<SF,SF>);
    module.def("combine", &_combine_<SF,VF>);
    module.def("combine", &_combine_<VF,SF>);
    module.def("combine", &_combine_<VF,VF>);
    module.def("embed", &_embed_<VF,D>);
        module.def("embed", &_embed_<D,VF,D>);

    module.def("restriction", &_restriction_<VF,I,Di>);
    module.def("restrict", &_restriction_<VF,I,Di>);
//    module.def("split", &_split_<VF,I>);

    module.def("evaluate", &_evaluate_<VF,VX>);
    module.def("partial_evaluate", &_partial_evaluate_<VF,I,X>);
    module.def("compose", &_compose_<VF,VF>);
    module.def("compose", &_compose_<SF,VF>);
    module.def("compose", &_compose_<SFN,VF>);
    module.def("compose", &_compose_<VFN,VF>);
    module.def("unchecked_compose", &_compose_<SF,VF>);
    module.def("unchecked_compose", &_compose_<VF,VF>);
    module.def("derivative", &_derivative_<VF,I>);
    module.def("antiderivative", &_antiderivative_<VF,I>);
    module.def("antiderivative", &_antiderivative_<VF,I,X>);

    module.def("compose", &_compose_<EffectiveScalarMultivariateFunction,VF>);
    module.def("compose", &_compose_<EffectiveVectorMultivariateFunction,VF>);

}


Void export_function_models(pybind11::module& module)
{

    export_approximate_taylor_model<FloatDP>(module);
    export_validated_taylor_model<FloatDP>(module);

    export_scalar_function_model<DP>(module);
    export_vector_function_model<DP>(module);
    export_scalar_taylor_function<FloatDP>(module);
    export_vector_taylor_function<FloatDP>(module);

    template_<ValidatedScalarMultivariateFunctionModel> scalar_function_model_template(module);
    scalar_function_model_template.instantiate<DP>();
    template_<ValidatedVectorMultivariateFunctionModel> vector_function_model_template(module);
    vector_function_model_template.instantiate<DP>();

    template_<ValidatedScalarMultivariateTaylorFunctionModel> scalar_taylor_function_model_template(module);
    scalar_taylor_function_model_template.instantiate<FloatDP>();
    scalar_taylor_function_model_template.def_new([](BoxDomainType dom,ValidatedScalarMultivariateFunction f,Sweeper<FloatDP> swp){return ValidatedScalarMultivariateTaylorFunctionModel<FloatDP>(dom,f,swp);});
    template_<ValidatedVectorMultivariateTaylorFunctionModel> vector_taylor_function_model_template(module);
    vector_taylor_function_model_template.instantiate<FloatDP>();
    vector_taylor_function_model_template.def_new([](BoxDomainType dom,ValidatedVectorMultivariateFunction f,Sweeper<FloatDP> swp){return ValidatedVectorMultivariateTaylorFunctionModel<FloatDP>(dom,f,swp);});



}
