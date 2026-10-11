#!/usr/bin/python3

from pyariadne import *


def section(title):
    print(f"\n== {title} ==")


def show(label, value):
    print(f"{label}: {value}")


def tutorial_geometry():
    dp = DoublePrecision()

    section("Points and boxes: cheap geometric enclosures")

    workspace = FloatDPExactBox([
        {exact(-2.0): exact(2.0)},
        {exact(-1.5): exact(1.5)},
    ])
    observation = FloatDPExactBox([
        {exact(-0.75): exact(1.25)},
        {exact(-0.5): exact(1.0)},
    ])
    sample = FloatDPPoint([FloatDP(0.5, dp), FloatDP(0.25, dp)])

    show("workspace", workspace)
    show("sample", sample)
    show("sample dimension", sample.dimension())
    show("sample bounding box", sample.bounding_box())
    show("workspace dimension", workspace.dimension())
    show("workspace centre", workspace.centre())
    show("workspace widths", workspace.widths())
    show("workspace radii", workspace.radii())
    show("workspace semiperimeter", workspace.semiperimeter())
    show("workspace area", workspace.measure())

    show("workspace contains sample", workspace.contains(sample))
    show("observation inside workspace", observation.inside(workspace))
    show("workspace covers observation", workspace.covers(observation))

    overlap = intersection(workspace, observation)
    envelope = hull(
        workspace,
        FloatDPExactBox([
            {exact(1.0): exact(3.0)},
            {exact(-0.25): exact(0.5)},
        ]),
    )
    lower, upper = observation.split()

    show("intersection", overlap)
    show("hull", envelope)
    show("lower half after split", lower)
    show("upper half after split", upper)

    section("Implicit regions: geometry defined by constraints")

    q0 = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    q1 = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    radius_squared = EffectiveVectorMultivariateFunction([
        sqr(q0) + sqr(q1),
    ])

    radial_codomain = RealBox(BoxDomainType([{exact(0): exact(1)}]))
    unit_disk = ConstraintSet(radius_squared, radial_codomain)

    search_window = RealBox(BoxDomainType([
        {exact(-1.25): exact(1.25)},
        {exact(-1.25): exact(1.25)},
    ]))
    bounded_disk = intersection(search_window, unit_disk)

    show("implicit unit disk", unit_disk)
    show("constraint count", unit_disk.number_of_constraints())
    show("bounded disk", bounded_disk)
    show("bounded disk domain", bounded_disk.domain())
    show("bounded disk bounding box", bounded_disk.bounding_box())
    show("disk overlaps observation", bounded_disk.overlaps(observation))

    shear = EffectiveVectorMultivariateFunction([
        q0 + q1 / 2,
        q1,
    ])
    sheared_disk = image(bounded_disk, shear)

    show("sheared disk", sheared_disk)
    show("sheared disk bounding box", sheared_disk.bounding_box())

    first_shear, second_shear = sheared_disk.split(0)
    show("first split image", first_shear.bounding_box())
    show("second split image", second_shear.bounding_box())

    section("Validated parametrised sets: dependence plus constraints")

    u0 = ValidatedScalarMultivariateFunction.coordinate(2, 0)
    u1 = ValidatedScalarMultivariateFunction.coordinate(2, 1)
    parametrisation = ValidatedVectorMultivariateFunction([
        u0 + u1,
        u0 * u1,
    ])
    parameter_domain = FloatDPExactBox([
        {exact(-1): exact(1)},
        {exact(-1): exact(1)},
    ])
    image_set = ValidatedConstrainedImageSet(parameter_domain, parametrisation)

    show("parameter domain", image_set.domain())
    show("parametrisation", image_set.function())
    show("image dimension", image_set.dimension())
    show("number of parameters", image_set.number_of_parameters())
    show("initial bounding box", image_set.bounding_box())

    image_set.new_parameter_constraint(u0 + u1 <= 1.25)

    x0 = ValidatedScalarMultivariateFunction.coordinate(2, 0)
    x1 = ValidatedScalarMultivariateFunction.coordinate(2, 1)
    state_constraint = x0 - x1 <= 1.5
    image_set.new_space_constraint(state_constraint)

    show("constrained image set", image_set)
    show("constraint count", image_set.number_of_constraints())
    show("constraint bounds", image_set.constraint_bounds())
    show("bounding box before reduction", image_set.bounding_box())

    image_set.reduce()
    show("original parameter domain", image_set.domain())
    show("reduced parameter domain", image_set.reduced_domain())
    show("bounding box after reduction", image_set.bounding_box())
    show("is empty", image_set.is_empty())

    first_image, second_image = image_set.split(0)
    show("first validated split", first_image.bounding_box())
    show("second validated split", second_image.bounding_box())

    left_parameter_domain = FloatDPExactBox([
        {exact(-1): exact(0)},
        {exact(-1): exact(1)},
    ])
    restricted = image_set.restriction(left_parameter_domain)
    show("restricted parameter domain", restricted.domain())
    show("restricted image bounds", restricted.bounding_box())

    query_box = FloatDPExactBox([
        {exact(-0.5): exact(0.5)},
        {exact(-0.5): exact(0.5)},
    ])
    show("overlaps query box", image_set.overlaps(query_box))
    show("separated from query box", image_set.separated(query_box))
    show("inside workspace", image_set.inside(workspace))

    section("Affine over-approximation: cheaper validated geometry")

    affine_set = image_set.affine_over_approximation()
    show("affine over-approximation", affine_set)
    show("affine dimension", affine_set.dimension())
    show("affine is bounded", affine_set.is_bounded())
    show("affine is empty", affine_set.is_empty())
    show("affine bounding box", affine_set.bounding_box())

    section("Interpolated curves: geometry from sampled trajectories")

    p0 = FloatDPPoint([FloatDP(0.0, dp), FloatDP(0.0, dp)])
    p1 = FloatDPPoint([FloatDP(0.5, dp), FloatDP(0.25, dp)])
    p2 = FloatDPPoint([FloatDP(1.0, dp), FloatDP(1.0, dp)])

    trajectory = InterpolatedCurve(FloatDP(0.0, dp), p0)
    trajectory.insert(FloatDP(0.5, dp), p1)
    trajectory.insert(FloatDP(1.0, dp), p2)

    show("trajectory", trajectory)
    show("stored samples", trajectory.size())
    show("trajectory dimension", trajectory.dimension())
    show("trajectory bounding box", trajectory.bounding_box())

    section("Grid cells and paving: finite set approximations")

    grid = Grid(2)
    show("grid", grid)
    show("grid origin", grid.origin())
    show("grid lengths", grid.lengths())

    primary = smallest_enclosing_primary_cell(image_set.bounding_box(), grid)
    show("smallest enclosing primary cell", primary)
    show("primary cell box", primary.box())

    coarse = image_set.outer_approximation(grid, 3)
    coarse.recombine()

    fine = image_set.outer_approximation(grid, 5)
    fine.recombine()

    show("coarse paving cells", coarse.size())
    show("coarse paving measure", coarse.measure())
    show("coarse paving bounds", coarse.bounding_box())
    show("fine paving cells", fine.size())
    show("fine paving measure", fine.measure())
    show("fine paving bounds", fine.bounding_box())

    observation_paving = GridTreePaving(grid)
    observation_paving.adjoin_over_approximation(observation, 4)
    observation_paving.recombine()

    common = intersection(fine, observation_paving)
    outside_observation = difference(fine, observation_paving)
    combined = union(fine, observation_paving)

    show("intersection cells", common.size())
    show("difference cells", outside_observation.size())
    show("union cells", combined.size())

    section("Choosing a geometry representation")
    print("Point / Box: exact locations and cheap coordinate-aligned enclosures.")
    print("ConstraintSet / BoundedConstraintSet: implicit feasible regions.")
    print("ConstrainedImageSet: effective images of constrained regions.")
    print("ValidatedConstrainedImageSet: nonlinear dependence with rigorous predicates.")
    print("ValidatedAffineConstrainedImageSet: cheaper correlated enclosure.")
    print("InterpolatedCurve: sampled paths and trajectories.")
    print("GridTreePaving: finite cell representation for subdivision and set algebra.")


if __name__ == "__main__":
    tutorial_geometry()
