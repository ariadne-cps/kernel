#ifndef ARIADNE_LABELLED_SET_HPP
#define ARIADNE_LABELLED_SET_HPP

#include "symbolic/expression_set.hpp"
#include "labelled/graphics_interface.hpp"

namespace Ariadne {

//! \brief A set formed by labelling the variables of a Euclidean set of type \a S.
template<class S> class LabelledSet {
    RealSpace _spc;
    S _set;
  public:
    typedef S EuclideanSetType;

    LabelledSet(const RealSpace& spc, const EuclideanSetType& set) : _spc(spc), _set(set) { ARIADNE_ASSERT_MSG(spc.dimension()==set.dimension(),"spc="<<spc<<", set.dimension()="<<set.dimension()); }
    template<class X> requires Same<S,Point<X>>
    explicit LabelledSet(VariablesPoint<X> const& point)
        : _spc(List<RealVariable>(point.variables())), _set(point.euclidean_set(_spc)) { }
    template<class IVL> requires Same<S,Box<IVL>>
    explicit LabelledSet(VariablesBox<IVL> const& box)
        : _spc(List<RealVariable>(box.variables())), _set(box.euclidean_set(_spc)) { }
    Set<RealVariable> variables() const { return _spc.variables(); }
    RealSpace const& space() const { return this->_spc; }
    EuclideanSetType const& continuous_set() const { return this->_set; }
    EuclideanSetType const& euclidean_set() const { return this->_set; }
    EuclideanSetType& euclidean_set() { return this->_set; }
    EuclideanSetType euclidean_set(const RealSpace& vars) const {
        Array<SizeType> prj(vars.dimension()); for(SizeType i=0; i!=vars.dimension(); ++i) { prj[i]=this->_spc.index(vars.variable(i)); } return project(this->_set,prj); }
    DimensionType dimension() const { return this->euclidean_set().dimension(); }
    template<class T> Void adjoin(LabelledSet<T> const& other) const {
        ARIADNE_ASSERT(this->space()==other.space()); this->euclidean_set().adjoin(other.euclidean_set()); }
    decltype(auto) bounding_box() const { return LabelledSet(this->space(),this->euclidean_set().bounding_box()); }
    decltype(auto) bounding_boxes() const { return LabelledSet(this->space(),this->euclidean_set().bounding_boxes()); }
    Void draw(CanvasInterface& cnvs, Variables2d const& vars) const;

    template<class E> friend decltype(auto) widen(const LabelledSet<S>& set, const E& e) {
        return LabelledSet(set.space(),widen(set.euclidean_set(),e)); }
    friend OutputStream& operator<<(OutputStream& os, const LabelledSet<S>& eset) {
        return os << "LabelledSet( space=" << eset.space() << ", set=" << eset.continuous_set() << " )"; }
};
template<class S> LabelledSet(const RealSpace&, const S&) -> LabelledSet<S>;

template<class IVL> class LabelledBox : public LabelledSet<Box<IVL>> {
  public:
    LabelledBox(const RealSpace& spc, const Box<IVL>& set) : LabelledSet<Box<IVL>>(spc,set) { }
    IVL operator[](const RealVariable& v) const { return this->euclidean_set()[this->space().index(v)]; }
};

template<class S> EqualsType<S> operator==(LabelledSet<S> const& eset1, LabelledSet<S> const& eset2) {
    if(eset1.variables()!=eset2.variables()) { return false; }
    ARIADNE_ASSERT(eset1.space()==eset2.space());
    return eset1.euclidean_set()==eset2.euclidean_set();
}

template<class UB> using LabelledInterval = VariableInterval<UB>;
template<class BS> using LabelledListSet = LabelledSet<ListSet<BS>>;


} // namespace Ariadne

#endif
