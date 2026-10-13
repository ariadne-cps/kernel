#include "labelled/figure.hpp"

#include <memory>

namespace Ariadne {

namespace {

ApproximateIntervalType coordinate_interval(ApproximateDoubleInterval const& ivl) {
    return ApproximateIntervalType(
        FloatDPApproximation(cast_exact(ivl.lower_bound()),dp),
        FloatDPApproximation(cast_exact(ivl.upper_bound()),dp));
}

Bool has_axis_bounds(Variables2d const& vars, Map<RealVariable,ApproximateDoubleInterval> const& bounds) {
    return bounds.find(vars.x())!=bounds.end() && bounds.find(vars.y())!=bounds.end();
}

ApproximateBoxType coordinate_box(Variables2d const& vars, Map<RealVariable,ApproximateDoubleInterval> const& bounds) {
    ARIADNE_ASSERT(has_axis_bounds(vars,bounds));
    return ApproximateBoxType({
        coordinate_interval(bounds[vars.x()]),
        coordinate_interval(bounds[vars.y()])
    });
}

class LabelledDrawableAdapter final : public Drawable2dInterface {
  public:
    LabelledDrawableAdapter(LabelledDrawable2dInterface const& shape, Variables2d const& variables)
        : _shape(shape.labelled_clone()), _variables(variables) { }

    LabelledDrawableAdapter* clone() const override { return new LabelledDrawableAdapter(*this); }
    DimensionType dimension() const override { return 2u; }

    Void draw(CanvasInterface& canvas, Projection2d const&) const override {
        _shape->draw(canvas,_variables);
    }

  private:
    std::shared_ptr<const LabelledDrawable2dInterface> _shape;
    Variables2d _variables;
};

} // namespace

Variables2d::Variables2d(const RealVariable& x, const RealVariable& y) : _x(x.name()), _y(y.name()) { }
RealVariable Variables2d::x() const { return RealVariable(_x); }
RealVariable Variables2d::y() const { return RealVariable(_y); }

Axes2d::Axes2d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval& y)
    : variables(x.variable(),y.variable()), bounds() {
    bounds.insert(x.variable(),x.interval());
    bounds.insert(y.variable(),y.interval());
}

Axes2d::Axes2d(ApproximateDouble xl, const RealVariable& x, ApproximateDouble xu,
               ApproximateDouble yl, const RealVariable& y, ApproximateDouble yu)
    : variables(x,y), bounds() {
    bounds.insert(x,ApproximateDoubleInterval(xl,xu));
    bounds.insert(y,ApproximateDoubleInterval(yl,yu));
}

Bool valid_axis_variables(const RealSpace& space, const Variables2d& variables) {
    return ((variables.x().name()==TimeVariable().name()) || space.contains(variables.x()))
        && space.contains(variables.y());
}

Projection2d projection(const RealSpace& space, const Variables2d& variables) {
    ARIADNE_ASSERT(valid_axis_variables(space,variables));
    DimensionType x_index=static_cast<DimensionType>(
        (variables.x()==TimeVariable() && !space.contains(variables.x()))
            ? space.dimension() : space.index(variables.x()));
    DimensionType y_index=static_cast<DimensionType>(space.index(variables.y()));
    return Projection2d(space.dimension(),x_index,y_index);
}

struct LabelledFigure::Data {
    Data(Variables2d const& vars, Map<RealVariable,ApproximateDoubleInterval> const& bnds)
        : variables(vars), bounds(bnds),
          figure(coordinate_box(vars,bnds),Projection2d(2u,0u,1u)) {
        figure.set_axis_labels(vars.x().name(),vars.y().name());
    }

    Void synchronise() {
        figure.set_axis_labels(variables.x().name(),variables.y().name());
        if(has_axis_bounds(variables,bounds)) {
            figure.set_bounding_box(coordinate_box(variables,bounds));
        }
    }

    Variables2d variables;
    Map<RealVariable,ApproximateDoubleInterval> bounds;
    Figure figure;
};

LabelledFigure::~LabelledFigure() { delete _data; }

LabelledFigure::LabelledFigure(const Axes2d& axes)
    : _data(new Data(axes.variables,axes.bounds)) { }

LabelledFigure::LabelledFigure(const Variables2d& vars, const VariablesBox<ApproximateIntervalType>& bbx)
    : _data(nullptr) {
    ARIADNE_ASSERT_MSG(bbx.bounds().size()>=2,"The box must be at least two-dimensional.")
    Map<RealVariable,ApproximateDoubleInterval> bounds=VariablesBox<ApproximateDoubleInterval>(bbx).bounds();
    _data=new Data(vars,bounds);
}

Void LabelledFigure::set_axes(const Axes2d& axes) {
    _data->variables=axes.variables;
    _data->bounds=axes.bounds;
    _data->synchronise();
}

Void LabelledFigure::set_bounds(const RealVariable& x, const ApproximateDoubleInterval& ivl) {
    _data->bounds.insert(x,ivl);
    _data->synchronise();
}

Void LabelledFigure::set_bounds(const Map<RealVariable,ApproximateDoubleInterval>& b) {
    _data->bounds=b;
    _data->synchronise();
}

Void LabelledFigure::set_bounding_box(const VariablesBox<ApproximateIntervalType>& bx) {
    _data->bounds=VariablesBox<ApproximateDoubleInterval>(bx).bounds();
    _data->synchronise();
}

GraphicsProperties& LabelledFigure::properties() { return _data->figure.properties(); }
GraphicsProperties const& LabelledFigure::properties() const { return _data->figure.properties(); }

LabelledFigure& LabelledFigure::draw(const LabelledDrawable2dInterface& shape) {
    LabelledDrawableAdapter adapter(shape,_data->variables);
    _data->figure.draw(adapter);
    return *this;
}

LabelledFigure& LabelledFigure::clear() {
    _data->figure.clear();
    return *this;
}

LabelledFigure& LabelledFigure::set_animated(Bool b) {
    _data->figure.set_animated(b);
    return *this;
}

Void LabelledFigure::_paint_all(CanvasInterface& canvas) const {
    _data->figure._paint_all(canvas);
}

Void LabelledFigure::write(const Char* filename) const {
    _data->figure.write(filename);
}

Void LabelledFigure::write(const Char* filename, Nat nx, Nat ny) const {
    _data->figure.write(filename,nx,ny);
}

Void draw(LabelledFigure& fig, const LabelledDrawable2dInterface& shape) {
    fig.draw(shape);
}

LabelledFigure& operator<<(LabelledFigure& fig, const LabelledDrawable2dInterface& shape) {
    return fig.draw(shape);
}

} // namespace Ariadne
