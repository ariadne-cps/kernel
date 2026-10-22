/***************************************************************************
 *            formula.tpl.hpp
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

#ifndef ARIADNE_FORMULA_TPL_HPP
#define ARIADNE_FORMULA_TPL_HPP

#include "function/operation_node.hpp"

namespace Ariadne {


template<class Y> using ConstantFormulaNode = OperationNode<Cnst,Y>;
template<class Y> class IndexFormulaNode : public OperationNode<Var,Index> { public: using OperationNode<Var,Index>::OperationNode; };
template<class Y> using UnaryFormulaNode = OperationNode<UnaryElementaryOperator,Formula<Y>>;
template<class Y> using BinaryFormulaNode = OperationNode<BinaryElementaryOperator,Formula<Y>,Formula<Y>>;
template<class Y> using GradedFormulaNode = OperationNode<GradedElementaryOperator,Formula<Y>,Int>;
template<class Y> using ScalarFormulaNode = OperationNode<BinaryElementaryOperator,Y,Formula<Y>>;

template<class Y> using FormulaNodeVariantType = Variant<ConstantFormulaNode<Y>,IndexFormulaNode<Y>, UnaryFormulaNode<Y>,BinaryFormulaNode<Y>,GradedFormulaNode<Y>,ScalarFormulaNode<Y>>;

template<class Y>
class FormulaNode
    : public FormulaNodeVariantType<Y>
{
  public:
    template<class FN> FormulaNode(FN&& fn) : FormulaNodeVariantType<Y>(std::forward<FN>(fn)) { }
    template<class VIS> decltype(auto) accept(VIS&& vis) const {
        return std::visit(std::forward<VIS>(vis),static_cast<FormulaNodeVariantType<Y> const&>(*this)); }
    FormulaNodeVariantType<Y> const& base() const { return *this; }
};



namespace FormulaNodeDetail {

template<class F, class J> inline auto derivative_impl(Abs, F const&, J) -> F {
    ARIADNE_THROW(std::runtime_error,"derivative(abs(f))","Cannot take derivative of non-smooth function"); }
template<class OP, class F, class J> inline decltype(auto) derivative_impl(OP op, F const& f, J j) {
    return op.derivative(f,derivative(f,j)); }
template<class F1, class F2, class J> inline decltype(auto) derivative_impl(Add, F1 const& f1, F2 const& f2, J j) {
    return derivative(f1,j)+derivative(f2,j); }
template<class F1, class F2, class J> inline decltype(auto) derivative_impl(Sub, F1 const& f1, F2 const& f2, J j) {
    return derivative(f1,j)-derivative(f2,j); }
template<class F1, class F2, class J> inline decltype(auto) derivative_impl(Mul, F1 const& f1, F2 const& f2, J j) {
    return derivative(f1,j)*f2+f1*derivative(f2,j); }
template<class F1, class F2, class J> inline decltype(auto) derivative_impl(Div, F1 const& f1, F2 const& f2, J j) {
    return (derivative(f1,j)-derivative(f2,j)*(f1/f2))/f2; }
template<class F1, class F2, class J> inline auto derivative_impl(Max, F1 const& f1, F2 const& f2, J) -> decltype(max(f1,f2)) {
    ARIADNE_THROW(std::runtime_error,"derivative(max(f1,f2))","Cannot take derivative of non-smooth function."); }
template<class F1, class F2, class J> inline auto derivative_impl(Min, F1 const& f1, F2 const& f2, J) -> decltype(min(f1,f2)) {
    ARIADNE_THROW(std::runtime_error,"derivative(min(f1,f2))","Cannot take derivative of non-smooth function."); }
template<class F, class N, class J> inline decltype(auto) derivative_impl(Pow op, F const& f, N const& n, J j) {
    return op.derivative(f,derivative(f,j),n); }
template<class F, class J> decltype(auto) derivative_impl(UnaryElementaryOperator ops, F const& f, J j) {
    return ops.accept([&](auto op){return derivative_impl(op,f,j);}); }
template<class F1, class F2, class J> decltype(auto) derivative_impl(BinaryElementaryOperator ops, F1 const& f1, F2 const& f2, J j) {
    return ops.accept([&](auto op){return derivative_impl(op,f1,f2,j);}); }
template<class F, class N, class J> decltype(auto) derivative_impl(GradedElementaryOperator ops, F const& f, N n, J j) {
    return ops.accept([&](auto op){return derivative_impl(op,f,n,j);}); }

template<class Y, class J> Formula<Y> derivative_node(ConstantFormulaNode<Y> const&, J) {
    return Formula<Y>::constant(Y(0)); }
template<class Y, class J> Formula<Y> derivative_node(IndexFormulaNode<Y> const& node, J j) {
    return Formula<Y>::constant(Y(Nat(node._ind)==Nat(j)?1:0)); }
template<class Y, class J> Formula<Y> derivative_node(UnaryFormulaNode<Y> const& node, J j) {
    return derivative_impl(node._op,node._arg,j); }
template<class Y, class J> Formula<Y> derivative_node(BinaryFormulaNode<Y> const& node, J j) {
    return derivative_impl(node._op,node._arg1,node._arg2,j); }
template<class Y, class J> Formula<Y> derivative_node(GradedFormulaNode<Y> const& node, J j) {
    return derivative_impl(node._op,node._arg,node._num,j); }
template<class Y, class J> Formula<Y> derivative_node(ScalarFormulaNode<Y> const& node, J j) {
    return derivative_impl(node._op,Formula<Y>::constant(node._cnst),node._arg,j); }

template<class Y> Formula<Y> substitute_node(ConstantFormulaNode<Y> const& node, Nat, Formula<Y> const&) {
    return Formula<Y>::constant(node._val); }
template<class Y> Formula<Y> substitute_node(IndexFormulaNode<Y> const& node, Nat i, Formula<Y> const& s) {
    return Nat(node._ind)==i ? s : Formula<Y>::coordinate(Nat(node._ind)); }
template<class Y> Formula<Y> substitute_node(UnaryFormulaNode<Y> const& node, Nat i, Formula<Y> const& s) {
    return Formula<Y>::unary(node._op,substitute(node._arg,i,s)); }
template<class Y> Formula<Y> substitute_node(BinaryFormulaNode<Y> const& node, Nat i, Formula<Y> const& s) {
    return Formula<Y>::binary(node._op,substitute(node._arg1,i,s),substitute(node._arg2,i,s)); }
template<class Y> Formula<Y> substitute_node(GradedFormulaNode<Y> const& node, Nat i, Formula<Y> const& s) {
    return Formula<Y>::graded(node._op,substitute(node._arg,i,s),node._num); }
template<class Y> Formula<Y> substitute_node(ScalarFormulaNode<Y> const& node, Nat i, Formula<Y> const& s) {
    return Formula<Y>::scalar(node._op,node._cnst,substitute(node._arg,i,s)); }

template<class E> inline E simpl(Add op, E const& e1, E const& e2) {
    if (is_constant(e1,0)) return e2;
    if (is_constant(e2,0)) return e1;
    return op(e1,e2); }
template<class E> inline E simpl(Sub op, E const& e1, E const& e2) {
    if (is_constant(e1,0)) return simplify(neg(e2));
    if (is_constant(e2,0)) return e1;
    if (identical(e1,e2)) return E::constant(0);
    return op(e1,e2); }
template<class E> inline E simpl(Mul op, E const& e1, E const& e2) {
    if (is_constant(e1,0) || is_constant(e2,1)) return e1;
    if (is_constant(e2,0) || is_constant(e1,1)) return e2;
    return op(e1,e2); }
template<class E> inline E simpl(Div op, E const& e1, E const& e2) {
    if (is_constant(e1,0) || is_constant(e2,1)) return e1;
    if (is_constant(e1,1)) return simplify(rec(e2));
    if (identical(e1,e2)) return E::constant(1);
    return op(e1,e2); }
template<class E> inline E simpl(Max op, E const& e1, E const& e2) { return op(e1,e2); }
template<class E> inline E simpl(Min op, E const& e1, E const& e2) { return op(e1,e2); }
template<class E> inline E simpl(Pow, E const& e, Int n) {
    switch(n) {
        case -1: return simplify(rec(e));
        case 0: return E::constant(1);
        case 1: return e;
        case 2: return simplify(sqr(e));
        default: return pow(e,n);
    }
}

template<class Y> Formula<Y> simplify_node(ConstantFormulaNode<Y> const& node) {
    return Formula<Y>::constant(node._val); }
template<class Y> Formula<Y> simplify_node(IndexFormulaNode<Y> const& node) {
    return Formula<Y>::coordinate(Nat(node._ind)); }
template<class Y> Formula<Y> simplify_node(UnaryFormulaNode<Y> const& node) {
    auto arg=simplify(node._arg);
    return node._op.accept([&](auto op) -> Formula<Y> {
        if constexpr (Same<decltype(op),Pos>) {
            return arg;
        } else if (arg.code()==Cnst::code()) {
            return Formula<Y>(op(arg.val()));
        } else {
            auto const* inverse_node=std::get_if<UnaryFormulaNode<Y>>(&arg.node_ref().base());
            if (inverse_node && are_inverses(node._op,inverse_node->_op)) {
                return inverse_node->_arg;
            }
            return Formula<Y>(op(arg));
        }
    }); }
template<class Y> Formula<Y> simplify_node(BinaryFormulaNode<Y> const& node) {
    auto arg1=simplify(node._arg1);
    auto arg2=simplify(node._arg2);
    return node._op.accept([&](auto op){return simpl(op,arg1,arg2);}); }
template<class Y> Formula<Y> simplify_node(GradedFormulaNode<Y> const& node) {
    auto arg=simplify(node._arg);
    return node._op.accept([&](auto op){return simpl(op,arg,node._num);}); }
template<class Y> Formula<Y> simplify_node(ScalarFormulaNode<Y> const& node) {
    auto arg=simplify(node._arg);
    Formula<Y> cnst=Formula<Y>::constant(node._cnst);
    return node._op.accept([&](auto op){return simpl(op,cnst,arg);}); }

template<class Y> Bool identical_node(ConstantFormulaNode<Y> const& a, ConstantFormulaNode<Y> const& b) {
    return same(a._val,b._val); }
template<class Y> Bool identical_node(IndexFormulaNode<Y> const& a, IndexFormulaNode<Y> const& b) {
    return a._ind==b._ind; }
template<class Y> Bool identical_node(UnaryFormulaNode<Y> const& a, UnaryFormulaNode<Y> const& b) {
    return a._op.code()==b._op.code() && identical(a._arg,b._arg); }
template<class Y> Bool identical_node(BinaryFormulaNode<Y> const& a, BinaryFormulaNode<Y> const& b) {
    return a._op.code()==b._op.code() && identical(a._arg1,b._arg1) && identical(a._arg2,b._arg2); }
template<class Y> Bool identical_node(GradedFormulaNode<Y> const& a, GradedFormulaNode<Y> const& b) {
    return a._op.code()==b._op.code() && a._num==b._num && identical(a._arg,b._arg); }
template<class Y> Bool identical_node(ScalarFormulaNode<Y> const& a, ScalarFormulaNode<Y> const& b) {
    return a._op.code()==b._op.code() && same(a._cnst,b._cnst) && identical(a._arg,b._arg); }
template<class A, class B> Bool identical_node(A const&, B const&) { return false; }

template<class Y> Bool is_constant_node(ConstantFormulaNode<Y> const&) { return true; }
template<class Y> Bool is_constant_node(IndexFormulaNode<Y> const&) { return false; }
template<class Y> Bool is_constant_node(UnaryFormulaNode<Y> const&) { return false; }
template<class Y> Bool is_constant_node(BinaryFormulaNode<Y> const&) { return false; }
template<class Y> Bool is_constant_node(GradedFormulaNode<Y> const&) { return false; }
template<class Y> Bool is_constant_node(ScalarFormulaNode<Y> const&) { return false; }

template<class Y> Bool is_constant_in_node(ConstantFormulaNode<Y> const&, Set<Nat> const&) { return true; }
template<class Y> Bool is_constant_in_node(IndexFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return !vars.contains(Nat(node._ind)); }
template<class Y> Bool is_constant_in_node(UnaryFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return is_constant_in(node._arg,vars); }
template<class Y> Bool is_constant_in_node(BinaryFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return is_constant_in(node._arg1,vars) && is_constant_in(node._arg2,vars); }
template<class Y> Bool is_constant_in_node(GradedFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return node._num==0 || is_constant_in(node._arg,vars); }
template<class Y> Bool is_constant_in_node(ScalarFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return is_constant_in(node._arg,vars); }

template<class E, class VARS> inline Bool affine_impl(Atan, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Acos, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Asin, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Cos, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Sin, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Tan, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Tanh, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Abs, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Sqr, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Log, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Sqrt, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Exp, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Rec, E const& e, VARS const& vars) { return is_constant_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Hlf, E const& e, VARS const& vars) { return is_affine_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Pos, E const& e, VARS const& vars) { return is_affine_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Neg, E const& e, VARS const& vars) { return is_affine_in(e,vars); }
template<class E, class VARS> inline Bool affine_impl(Nul, E const&, VARS const&) { return true; }
template<class E1, class E2, class VARS> inline Bool affine_impl(Add, E1 const& a, E2 const& b, VARS const& vars) {
    return is_affine_in(a,vars) && is_affine_in(b,vars); }
template<class E1, class E2, class VARS> inline Bool affine_impl(Sub, E1 const& a, E2 const& b, VARS const& vars) {
    return is_affine_in(a,vars) && is_affine_in(b,vars); }
template<class E1, class E2, class VARS> inline Bool affine_impl(Mul, E1 const& a, E2 const& b, VARS const& vars) {
    return (is_affine_in(a,vars)&&is_constant_in(b,vars)) || (is_constant_in(a,vars)&&is_affine_in(b,vars)); }
template<class E1, class E2, class VARS> inline Bool affine_impl(Div, E1 const& a, E2 const& b, VARS const& vars) {
    return is_affine_in(a,vars) && is_constant_in(b,vars); }
template<class E1, class E2, class VARS> inline Bool affine_impl(Max, E1 const& a, E2 const& b, VARS const& vars) {
    return is_constant_in(a,vars) && is_constant_in(b,vars); }
template<class E1, class E2, class VARS> inline Bool affine_impl(Min, E1 const& a, E2 const& b, VARS const& vars) {
    return is_constant_in(a,vars) && is_constant_in(b,vars); }
template<class E, class N, class VARS> inline Bool affine_impl(Pow, E const& e, N n, VARS const& vars) {
    return n==0 || (n==1 && is_affine_in(e,vars)) || is_constant_in(e,vars); }

template<class Y> Bool is_affine_in_node(ConstantFormulaNode<Y> const&, Set<Nat> const&) { return true; }
template<class Y> Bool is_affine_in_node(IndexFormulaNode<Y> const&, Set<Nat> const&) { return true; }
template<class Y> Bool is_affine_in_node(UnaryFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return node._op.accept([&](auto op){return affine_impl(op,node._arg,vars);}); }
template<class Y> Bool is_affine_in_node(BinaryFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return node._op.accept([&](auto op){return affine_impl(op,node._arg1,node._arg2,vars);}); }
template<class Y> Bool is_affine_in_node(GradedFormulaNode<Y> const& node, Set<Nat> const& vars) {
    return node._op.accept([&](auto op){return affine_impl(op,node._arg,node._num,vars);}); }
template<class Y> Bool is_affine_in_node(ScalarFormulaNode<Y> const& node, Set<Nat> const& vars) {
    Formula<Y> cnst=Formula<Y>::constant(node._cnst);
    return node._op.accept([&](auto op){return affine_impl(op,cnst,node._arg,vars);}); }

} // namespace FormulaNodeDetail


template<class Y> inline Formula<Y> make_formula(const Y& c) {
    return Formula<Y>::constant(c); }
template<class Y> inline Formula<Y> make_formula(Cnst, const Y& c) {
    return Formula<Y>::constant(c); }
template<class Y> inline Formula<Y> make_formula(Ind, SizeType j) {
    return Formula<Y>::index(j); }
template<class Y> inline Formula<Y> make_formula(const UnaryElementaryOperator& op, const Formula<Y>& arg) {
    return Formula<Y>::unary(op,arg); }
template<class Y> inline Formula<Y> make_formula(const BinaryElementaryOperator& op, const Formula<Y>& arg1, const Formula<Y>& arg2) {
    return Formula<Y>::binary(op,arg1,arg2); }
template<class Y> inline Formula<Y> make_formula(const BinaryElementaryOperator& op, const Formula<Y>& arg1, const Y& arg2) {
    return Formula<Y>::binary(op,arg1,make_formula(arg2)); }
template<class Y> inline Formula<Y> make_formula(const GradedElementaryOperator& op, const Formula<Y>& arg, Int num) {
    return Formula<Y>::graded(op,arg,num); }
template<class Y, class OP> inline Formula<Y> make_formula(const OP& op, Y const& cnst, const Formula<Y>& arg) {
    OperatorCode op_code=static_cast<OperatorCode>((char)op.code()+((char)OperatorCode::SADD-(char)OperatorCode::ADD));
    return Formula<Y>::scalar(Operator(op_code,OperatorKind::SCALAR),cnst,arg); }

template<class Y> Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg) {
    return make_formula(UnaryElementaryOperator(op.code()),arg); }
template<class Y> Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg1, const Formula<Y>& arg2) {
    return make_formula(BinaryElementaryOperator(op.code()),arg1,arg2); }
template<class Y> Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg1, Int n2) {
    return make_formula(GradedElementaryOperator(op.code()),arg1,n2); }
template<class Y> Formula<Y> make_formula(const Operator& op, const Y& c1, const Formula<Y>& arg2) {
    return make_formula(BinaryElementaryOperator(op.code()),c1,arg2); }


namespace FormulaDetail {
template<class X, class Y> decltype(auto) evaluate(OperationNode<Cnst,Y> c, Vector<X> const& v) {
    return make_constant(c._val,v); }
template<class X, class I> decltype(auto) evaluate(OperationNode<Var,I> i, Vector<X> const& v) {
    return v[i._ind]; }
template<class X, class OP, class A> decltype(auto) evaluate(OperationNode<OP,A> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg,v)); }
template<class X, class OP, class A1, class A2> decltype(auto) evaluate(OperationNode<OP,A1,A2> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg1,v),evaluate(s._arg2,v)); }
template<class X, class OP, class A1, class A2, class A3> decltype(auto) evaluate(OperationNode<OP,A1,A2,A3> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg1,v),evaluate(s._arg2,v),evaluate(s._arg3,v)); }
template<class X, class OP, class A> decltype(auto) evaluate(OperationNode<OP,A,Formula<A>> s, Vector<X> const& v) {
    return s._op(evaluate(OperationNode<Cnst,A>(s._cnst),v),evaluate(s._arg,v)); }
template<class X, class OP, class A> decltype(auto) evaluate(OperationNode<OP,A,Int> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg,v),s._num); }

template<class X, class Y> decltype(auto) _cached_evaluate(OperationNode<Cnst,Y> const& c, Vector<X> const& v, Map<const Void*,X>&) {
    return make_constant(c._val,v); }
template<class X, class I> decltype(auto) _cached_evaluate(OperationNode<Var,I> const& i, Vector<X> const& v, Map<const Void*,X>&) {
    return v[i._ind]; }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(OperationNode<OP,A> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg,v,cache)); }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(OperationNode<OP,A,Formula<A>> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(evaluate(ConstantFormulaNode<A>(s._cnst),v),cached_evaluate(s._arg,v,cache)); }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(OperationNode<OP,Formula<A>,A> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(s._cnst,cached_evaluate(s._arg,v,cache)); }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(OperationNode<OP,A,Int> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg,v,cache),s._num); }
template<class X, class OP, class A1, class A2> decltype(auto) _cached_evaluate(OperationNode<OP,A1,A2> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg1,v,cache),cached_evaluate(s._arg2,v,cache)); }
template<class X, class OP, class A1, class A2, class A3> decltype(auto) _cached_evaluate(OperationNode<OP,A1,A2,A3> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg1,v,cache),cached_evaluate(s._arg2,v,cache),cached_evaluate(s._arg3,v,cache)); }

template<class X, class OP, class... AS> decltype(auto) _cached_evaluate(OperationNode<OP,AS...> const&, Vector<X> const& v, Map<const Void*,X>&) {
    return v[0]; }

} // namespace FormulaDetail


template<class X, class Y> X direct_evaluate(const Formula<Y>& f, const Vector<X>& x) {
    return f.node_ref().accept([&x](auto s){return FormulaDetail::evaluate(s,x);});
}


//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
template<class X, class Y> const X& cached_evaluate(const Formula<Y>& f, const Vector<X>& x, Map<const Void*,X>& cache) {
    const FormulaNode<Y>* fptr=f.node_ptr();
    if(cache.has_key(fptr)) { return cache.get(fptr); }
    else { X r=f.node_ref().accept([&x,&cache](auto s){return FormulaDetail::_cached_evaluate(s,x,cache);}); return insert(cache,fptr,r); }
}

template<class X, class Y> inline X cached_evaluate(const Formula<Y>& f, const Vector<X>& v) {
    Map<const Void*,X> cache;
    return cached_evaluate(f,v,cache);
}

template<class X, class Y> Vector<X> cached_evaluate(const Vector<Formula<Y>>& f, const Vector<X>& v) {
    Vector<X> r(f.size(),zero_element(v));
    Map<const Void*,X> cache;
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=cached_evaluate(f[i],v,cache);
    }
    return r;
}

template<class X, class Y> X evaluate(const Formula<Y>& f, const Vector<X>& x) {
    return cached_evaluate(f,x);
}
template<class X, class Y> Vector<X> evaluate(const Vector<Formula<Y>>& f, const Vector<X>& x) {
    return cached_evaluate(f,x);
}




//! \ingroup FunctionModule
//! \brief Convert a power-series expansion into a formula using a version of Horner's rule.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
//!
template<class X> Formula<X> formula(const Expansion<MultiIndex,X>& e)
{
    Vector<Formula<X>> identity(e.argument_size());
    for(SizeType i=0; i!=identity.size(); ++i) { identity[i]=Formula<X>::coordinate(i); }
    return horner_evaluate(e,identity);
}

namespace {
template<class X, class Y> inline const Formula<Y>& _substitute_coordinate(const Nat& ie, const Nat& is, const Formula<Y>& a, const Formula<X>& s) {
    ARIADNE_ASSERT_MSG(ie!=is,"Cannot substitute formula "<<s<<" for coordinate "<<ie);
    return a; }
template<class Y> inline const Formula<Y>& _substitute_coordinate(const Nat& ie, const Nat& is, const Formula<Y>& e, const Formula<Y>& s) {
    return ie==is ? s : e; }
} // namespace

template<class Y> Formula<Y> substitute(const Formula<Y>& a, const Nat& i, const Formula<Y>& is) {
    return a.node_ref().accept([&i,&is](auto fn){return FormulaNodeDetail::substitute_node(fn,i,is);});
}


template<class Y> inline Bool is_constant(const Formula<Y>& f, const SelfType<Y>& c) {
    auto const* cfn = std::get_if<ConstantFormulaNode<Y>>(&f.node_ref()); return cfn && same(cfn->val(),c);
}

template<class Y> inline Formula<Y> simplify(const Formula<Y>& a) {
    return a.node_ref().accept([](auto fn){return FormulaNodeDetail::simplify_node(fn);});
}



inline Bool same(EffectiveNumber const& v1, EffectiveNumber const& v2) {
    // FIXME: Use symbolic approach
    DoublePrecision pr;
    FloatDPBounds x1(v1,pr);
    FloatDPBounds x2(v2,pr);
    return x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw();
}

template<class Y> Bool identical(const Formula<Y>& a1, const Formula<Y>& a2)
{
    return a1.node_ref().accept([&a2](auto const& n1){
        return a2.node_ref().accept([&n1](auto const& n2){return FormulaNodeDetail::identical_node(n1,n2);});
    });
}

template<class Y> Vector<Formula<Y>> simplify(const Vector<Formula<Y>>& f) {
    return Vector<Formula<Y>>(f.size(),[&f](SizeType i){return simplify(f[i]);},f.zero_element());
}


template<class Y> Bool identical(const Vector<Formula<Y>>& f1, const Vector<Formula<Y>>& f2) {
    if (f1.size() != f2.size()) { return false; }
    for (auto i : range(f1.size())) {
        if (not identical(f1[i],f2[i])) { return false; }
    }
    return true;
}

template<class Y> Bool is_constant(Formula<Y> const& f) { return f.node_ref().accept([](auto fn){return FormulaNodeDetail::is_constant_node(fn);}); }

template<class Y> Bool is_constant_in(const Formula<Y>& f, const Set<Nat>& is) {
    return f.node_ref().accept([&is](auto fn){return FormulaNodeDetail::is_constant_in_node(fn,is);});
}


template<class Y> Bool is_affine_in(const Formula<Y>& f, const Set<Nat>& is) {
    return f.node_ref().accept([&is](auto fn){return FormulaNodeDetail::is_affine_in_node(fn,is);});
}

template<class Y> Bool is_affine_in(const Vector<Formula<Y>>& fs, const Set<Nat>& is) {
    for (auto idx : range(fs.size())) {
        if (not is_affine_in(fs[idx],is)) { return false; }
    }
    return true;
}

template<class Y> Bool is_additive_in(const Vector<Formula<Y>>& fs, const Set<Nat>& is) {
    // We treat the vector of formulas as additive in is if each variable in is appears at most once in all expressions,
    // with a constant multiplier
    // (FIXME: this simplifies the case of a diagonalisable matrix of constant multipliers)

    for (auto i : is) {
        Bool already_found = false;
        for (auto idx : range(fs.size())) {
            const Formula<Y>& f = fs[idx];
            if (not is_constant_in(f,{i})) {
                if (already_found) {
                    return false;
                } else {
                    already_found = true;
                    auto dfi = simplify(derivative(f, i));
                    if (not is_constant(dfi)) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}


} // namespace Ariadne

#endif // ARIADNE_FORMULA_TPL_HPP
