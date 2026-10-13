#ifndef ARIADNE_LABELLED_FIGURE_HPP
#define ARIADNE_LABELLED_FIGURE_HPP

#include "io/figure.hpp"
#include "labelled/graphics_interface.hpp"
#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression_set.hpp"

namespace Ariadne {

using ApproximateDoubleVariableInterval = VariableInterval<ApproximateDouble>;
using DoubleVariableInterval = VariableInterval<Double>;

struct Variables2d {
    Identifier _x,_y;
    Variables2d(const RealVariable& x, const RealVariable& y);
    RealVariable x() const;
    RealVariable y() const;
};

struct Axes2d {
    Axes2d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval& y);
    Axes2d(ApproximateDouble xl, const RealVariable& x, ApproximateDouble xu,
           ApproximateDouble yl, const RealVariable& y, ApproximateDouble yu);
    Variables2d variables;
    Map<RealVariable,ApproximateDoubleInterval> bounds;
};

Projection2d projection(const RealSpace& space, const Variables2d& variables);

class LabelledFigure {
  public:
    ~LabelledFigure();

    LabelledFigure(const Axes2d& axes);
    LabelledFigure(const Variables2d& vars, VariablesBox<ApproximateIntervalType> const& bnds);

    Void set_axes(const Axes2d& axes);
    Void set_bounds(const RealVariable& x, const ApproximateDoubleInterval& ivl);
    Void set_bounds(const Map<RealVariable,ApproximateDoubleInterval>& b);
    Void set_bounding_box(VariablesBox<ApproximateIntervalType> const& bnds);

    GraphicsProperties& properties();
    GraphicsProperties const& properties() const;

    LabelledFigure& draw(const LabelledDrawable2dInterface& shape);
    LabelledFigure& clear();
    LabelledFigure& set_animated(Bool b);

    Void display() const;
    Void write(const Char* filename, Nat nx, Nat ny) const;
    Void write(const Char* filename) const;

    Void _paint_all(CanvasInterface& canvas) const;

  public:
    struct Data;
  private:
    Data* _data;
};

Void draw(LabelledFigure& fig, const LabelledDrawable2dInterface& shape);
LabelledFigure& operator<<(LabelledFigure& fig, const LabelledDrawable2dInterface& shape);

inline LabelledFigure& operator<<(LabelledFigure& g, const LineStyle& ls) { g.properties().set_line_style(ls); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const LineWidth& lw) { g.properties().set_line_width(lw); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const LineColour& lc) { g.properties().set_line_colour(lc); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const FillStyle& fs) { g.properties().set_fill_style(fs); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const FillOpacity& fo) { g.properties().set_fill_opacity(fo); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const FillColour& fc) { g.properties().set_fill_colour(fc); return g; }
inline LabelledFigure& operator<<(LabelledFigure& g, const Animated& animation) { g.properties().set_animated(animation); return g; }

template<class S> class LabelledDrawableWrapper : public LabelledDrawable2dInterface {
    LabelledSet<S> _lset;
  public:
    LabelledDrawableWrapper(LabelledSet<S> lset) : _lset(lset) { }
    LabelledDrawableWrapper<S>* labelled_clone() const override { return new LabelledDrawableWrapper<S>(*this); }
    Void draw(CanvasInterface& cnvs, Variables2d const& vars) const override {
        Projection2d prj(this->_lset.euclidean_set().dimension(),
                         this->_lset.space()[vars.x()],this->_lset.space()[vars.y()]);
        this->_lset.euclidean_set().draw(cnvs,prj);
    }
};

template<class S> inline LabelledFigure& operator<<(LabelledFigure& g, const LabelledSet<S>& lset) {
    return g << LabelledDrawableWrapper<S>(lset);
}

} // namespace Ariadne

#endif
