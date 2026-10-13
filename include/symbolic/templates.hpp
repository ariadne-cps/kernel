/***************************************************************************
 *            symbolic/templates.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file symbolic/templates.hpp
 *  \brief
 */


#ifndef ARIADNE_EXPRESSION_TEMPLATES_HPP
#define ARIADNE_EXPRESSION_TEMPLATES_HPP

#include "function/operation_node.hpp"
#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"

namespace Ariadne {


template<class T> class Expression;

struct Get {
    static constexpr OperatorCode code() { return OperatorCode::GET; }
    static OperatorKind kind() { return OperatorKind::COORDINATE; }
    template<class V, class I> decltype(auto) operator() (V const& v, I i) { return v.operator[](i); }
    friend OutputStream& operator<<(OutputStream& os, Get const&) { return os << "get"; }
};

struct Vec {
    static constexpr OperatorCode code() { return OperatorCode::VEC; }
    static OperatorKind kind() { return OperatorKind::NULLARY; }
    template<class X> decltype(auto) operator() (Array<X> const& a) { return Vector<X>(a); }
    friend OutputStream& operator<<(OutputStream& os, Vec const&) { return os << "Vector"; }
};

/************ Symbolic expression extensions ***************************/

template<class OP, class... AS> using Symbolic = OperationNode<OP,AS...>;

template<class O, class... Args> OperationNode<O,Args...> make_expression_template(O o, Args... as) {
    return OperationNode<O,Args...>(o,as...); }

template<class... SS> class SymbolicVariant : public Variant<SS...> {
    template<class VIS> decltype(auto) accept(VIS&& vis) const { return std::visit(std::forward<VIS>(vis),static_cast<Variant<SS...>const&>(*this)); }
};

struct While { };
struct Iterate { };
struct IfThnEls { };

template<class O> concept AGetter = Same<O,Get> or Same<O,OperatorVariant<Get>>;

template<class O, class V, class I> requires AGetter<O> struct OperationNode<O,V,I> {
    O _op; V _vec; I _ind;
    O op() const { return _op; } V arg() const { return _vec; } V arg1() const { return _vec; } I arg2() const { return _ind; }
    V vec() const { return _vec; } I ind() const { return _ind; }
    OperationNode(O o, V v, I i) : _op(o), _vec(v), _ind(i) { }
//    template<class T> explicit operator T() const { return _op(static_cast<T>(_vec),_ind); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_vec(vals...),_ind)) {
        return _op(_vec(vals...),_ind); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode<O,V,I> const& expr) {
        return os << expr._op.code() << "(" << expr._vec << "," << expr._ind << ")"; }
};

template<class C, class A> struct OperationNode<IfThnEls,C,A> {
    IfThnEls _op; C _cnd; A _atru; A _afls;
    OperationNode(IfThnEls op, C cnd, A atru, A afls) : _op(op), _cnd(cnd), _atru(atru), _afls(afls) { }
    template<class... VS> decltype(auto) operator () (VS... vals) const {
        if(_cnd(vals...)) { return  _atr(vals...); } else { return _afls(vals...); } }
    friend OutputStream& operator<<(OutputStream& os, OperationNode<IfThnEls,C,A> const& expr) {
        return os << "if(" << expr._cnd << ") { " << expr._atru << " } else { " << expr._afls << " }";
    }
};
template<class C, class F> struct OperationNode<While,C,F> {
    While _op; C _cnd; F _fn;
    OperationNode(While op, C cnd, F fn) : _op(op), _cnd(cnd), _fn(fn) { }
    template<class S> S& operator () (S& s) const {
        while(_cnd(s)) { s=_fn(s); } return s; }
    friend OutputStream& operator<<(OutputStream& os, OperationNode<While,C,F> const& expr) {
        return os << "while(" << expr._cnd << ") { " << expr._fn << " }";
    }
};

template<class F> struct OperationNode<Iterate,F> {
    F _fn;
    OperationNode(Iterate, F fn) : _fn(fn) { }
    template<class S> Sequence<S> operator () (S const& s) const {
        // TODO: Use cache so that sequence does not always have to be recomputed
        return Sequence<S>([&](Natural const& n){ S t=s; for(Natural i=0u; i!=n; ++i) { t=_fn(t); } return t; }); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode<Iterate,F> const& expr) {
        return os << "iterate { " << expr._fn << " }";
    }
};


template<class O, class... Args> struct TemporaryExpression;

template<class O, class A> struct TemporaryExpression<O,A> {
    O _op; A const& _arg;
    TemporaryExpression(O o, A const& a) : _op(o), _arg(a) { }
};

template<class O, class A1, class A2> struct TemporaryExpression<O,A1,A2> {
    O _op; A1 const& _arg1; A2 const& _arg2;
    TemporaryExpression(O o, A1 const& a1, A2 const& a2) : _op(o), _arg1(a1), _arg2(a2) { }
};

template<class O, class A1, class A2, class A3> struct TemporaryExpression<O,A1,A2,A3> {
    O _op; A1 const& _arg1; A2 const& _arg2; A3 const& _arg3;
    TemporaryExpression(O o, A1 const& a1, A2 const& a2, A3 const& a3) : _op(o), _arg1(a1), _arg2(a2), _arg3(a3) { }
};


} // namespace Ariadne

#endif
