#include <iostream>
#include <string>

#include "ariadne-kernel.hpp"

using namespace Ariadne;
using std::cout;
using std::endl;

namespace {

void section(const std::string& title) {
    cout << "\n== " << title << " ==\n";
}

template<class T>
void show(const std::string& label, const T& value) {
    cout << label << ": " << value << "\n";
}

} // namespace

int main() {
    section("Boxes as the basic geometric enclosure");

    ExactBoxType workspace({
        {-2.0_x, 2.0_x},
        {-1.5_x, 1.5_x}
    });
    ExactBoxType observation({
        {-0.75_x, 1.25_x},
        {-0.5_x, 1.0_x}
    });
    FloatDPPoint sample({ExactDouble(0.5), ExactDouble(0.25)}, dp);

    show("workspace", workspace);
    show("observation", observation);
    show("dimension", workspace.dimension());
    show("centre", workspace.centre());
    show("widths", workspace.widths());
    show("radii", workspace.radii());
    show("semiperimeter", workspace.semiperimeter());
    show("area", workspace.measure());
    show("sample", sample);
    show("workspace contains sample", workspace.contains(sample));
    show("observation inside workspace", observation.inside(workspace));
    show("workspace covers observation", workspace.covers(observation));

    // Boxes support the usual enclosure algebra.  intersection() keeps only
    // common points, hull() gives the smallest coordinate-aligned box
    // containing both operands, and split() is the primitive used by
    // subdivision algorithms.
    ExactBoxType overlap = intersection(workspace, observation);
    ExactBoxType envelope = hull(workspace, ExactBoxType({
        {1.0_x, 3.0_x},
        {-0.25_x, 0.5_x}
    }));
    auto halves = observation.split();

    show("intersection", overlap);
    show("hull", envelope);
    show("lower half after split", halves.first);
    show("upper half after split", halves.second);

    section("Sets defined by functions and constraints");

    // A geometric set need not be represented by a box.  Here parameters
    // u=(u0,u1) range over the unit box and are mapped to a curved subset
    // of R^2.  The set retains the parameter dependence instead of reducing
    // immediately to a coarse bounding box.
    auto u = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction parametrisation({
        u[0] + u[1],
        u[0] * u[1]
    });
    ExactBoxType parameter_domain = ExactBoxType::unit_box(2);
    ValidatedConstrainedImageSet image_set(parameter_domain, parametrisation);

    show("parameter domain", parameter_domain);
    show("parametrisation", parametrisation);
    show("image set", image_set);
    show("image dimension", image_set.dimension());
    show("number of parameters", image_set.number_of_parameters());
    show("bounding box", image_set.bounding_box());

    // Constraints can be imposed either in parameter space or after mapping
    // to the geometric space.  This is useful when a reachable or feasible
    // set has both a parametrisation and additional admissibility conditions.
    image_set.new_parameter_constraint(u[0] + u[1] <= 1.25_x);

    auto x = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveConstraint space_constraint = x[0] - x[1] <= 1.5_x;
    image_set.new_space_constraint(space_constraint);

    show("constrained image set", image_set);
    show("constraints", image_set.number_of_constraints());
    show("constrained bounding box", image_set.bounding_box());

    section("Affine over-approximation");

    // ValidatedAffineConstrainedImageSet is a cheaper geometric surrogate.
    // It is often useful when a downstream algorithm only needs a rigorous
    // affine enclosure of a more complicated constrained image.
    ValidatedAffineConstrainedImageSet affine_set(image_set.affine_over_approximation());
    show("affine over-approximation", affine_set);
    show("affine bounding box", affine_set.bounding_box());

    section("Grid cells and paving");

    // A Grid fixes the coordinate system used by discrete set
    // approximations.  The paving stores only active cells in a binary tree;
    // fineness controls how far the subdivision is refined.
    Grid grid(2);
    show("grid", grid);
    show("grid origin", grid.origin());
    show("grid lengths", grid.lengths());

    GridCell primary = GridCell::smallest_enclosing_primary_cell(
        image_set.bounding_box(), grid);
    show("smallest enclosing primary cell", primary);
    show("primary cell box", primary.box());

    GridTreePaving coarse(grid);
    image_set.adjoin_outer_approximation_to(coarse, 3);
    coarse.recombine();

    GridTreePaving fine(grid);
    image_set.adjoin_outer_approximation_to(fine, 5);
    fine.recombine();

    show("coarse paving cells", coarse.size());
    show("coarse paving measure", coarse.measure());
    show("coarse paving bounds", coarse.bounding_box());
    show("fine paving cells", fine.size());
    show("fine paving measure", fine.measure());
    show("fine paving bounds", fine.bounding_box());

    // Set operations are available directly on pavings.  Here the same set
    // is approximated at two resolutions; their intersection is therefore a
    // convenient way to demonstrate discrete set algebra without changing
    // the continuous set being represented.
    GridTreePaving common = intersection(coarse, fine);
    show("intersection of pavings", common);
    show("intersection cells", common.size());

    section("When to use each representation");
    cout
        << "Box: fast coordinate-aligned enclosure and predicates.\n"
        << "ValidatedConstrainedImageSet: keeps nonlinear parameter dependence and constraints.\n"
        << "ValidatedAffineConstrainedImageSet: rigorous affine surrogate for cheaper geometry.\n"
        << "GridTreePaving: finite cell representation for subdivision and set algebra.\n";

    return 0;
}
