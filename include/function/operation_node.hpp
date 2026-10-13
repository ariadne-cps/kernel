/***************************************************************************
 *            function/operation_node.hpp
 *
 *  Copyright  2013-26  Ariadne contributors
 *
 ****************************************************************************/

#ifndef ARIADNE_FUNCTION_OPERATION_NODE_HPP
#define ARIADNE_FUNCTION_OPERATION_NODE_HPP

#include "numeric/operators.hpp"

namespace Ariadne {

template<class OP, class... AS> struct OperationNode;

template<class C> struct OperationNode<Cnst,C> {
    Cnst _op; C _val;
    explicit OperationNode(C c) : _op(), _val(c) { }
    OperationNode(Cnst o, C c) : _op(o), _val(c) { }
    C val() const { return _val; }
    template<class T> explicit operator T() const { return static_cast<T>(_val); }
    template<class... AS> auto operator() (AS...) const -> C { return _val; }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) { return os << node._val; }
};

template<class I> struct OperationNode<Var,I> {
    Var _op; I _ind;
    explicit OperationNode(I i) : _op(), _ind(i) { }
    OperationNode(Var o, I i) : _op(o), _ind(i) { }
    I ind() const { return _ind; }
    template<class AS> auto operator() (AS vals) const -> decltype(vals[_ind]) { return vals[_ind]; }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) { return os << node._ind; }
};

template<class O, class A> struct OperationNode<O,A> {
    O _op; A _arg;
    O op() const { return _op; }
    A arg() const { return _arg; }
    OperationNode(O o, A a) : _op(o), _arg(a) { }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...))) {
        return _op(_arg(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) {
        return os << node._op.code() << "(" << node._arg << ")"; }
};

template<class O, class A1, class A2> struct OperationNode<O,A1,A2> {
    O _op; A1 _arg1; A2 _arg2;
    O op() const { return _op; }
    A1 arg1() const { return _arg1; }
    A2 arg2() const { return _arg2; }
    OperationNode(O o, A1 a1, A2 a2) : _op(o), _arg1(a1), _arg2(a2) { }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg1(vals...),_arg2(vals...))) {
        return _op(_arg1(vals...),_arg2(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) {
        return os << node._op.code() << "(" << node._arg1 << "," << node._arg2 << ")"; }
};

template<class O> concept AGraded = Same<O,Pow> or Same<O,GradedElementaryOperator>;

template<class O, class A, class N> requires AGraded<O> struct OperationNode<O,A,N> {
    O _op; A _arg; N _num;
    O op() const { return _op; }
    A arg() const { return _arg; }
    N num() const { return _num; }
    A arg1() const { return _arg; }
    N arg2() const { return _num; }
    OperationNode(O o, A a, N n) : _op(o), _arg(a), _num(n) { }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...),_num)) {
        return _op(_arg(vals...),_num); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) {
        return os << node._op.code() << "(" << node._arg << "," << node._num << ")"; }
};

template<class O, class A, template<class>class E> struct OperationNode<O,A,E<A>> {
    O _op; A _cnst; E<A> _arg;
    OperationNode(O o, A c, E<A> a) : _op(o), _cnst(c), _arg(a) { }
    O op() const { return _op; }
    A cnst() const { return _cnst; }
    E<A> arg() const { return _arg; }
    template<class T> operator T() const { return _op(static_cast<T>(_cnst),static_cast<T>(_arg)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_cnst,_arg(vals...))) {
        return _op(_cnst,_arg(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) {
        return os << node._op.code() << "(" << node._cnst << "," << node._arg << ")"; }
};

template<class O, class A1, class A2, class A3> struct OperationNode<O,A1,A2,A3> {
    O _op; A1 _arg1; A2 _arg2; A3 _arg3;
    OperationNode(O o, A1 a1, A2 a2, A3 a3) : _op(o), _arg1(a1), _arg2(a2), _arg3(a3) { }
    template<class T> explicit operator T() const {
        return _op(static_cast<T>(_arg1),static_cast<T>(_arg2),static_cast<T>(_arg3)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg1(vals...),_arg2(vals...),_arg3(vals...))) {
        return _op(_arg1(vals...),_arg2(vals...),_arg3(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, OperationNode const& node) {
        return os << node._op.code() << "(" << node._arg1 << "," << node._arg2 << "," << node._arg3 << ")"; }
};

} // namespace Ariadne

#endif
