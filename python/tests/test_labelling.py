from pyariadne import *


def test_labelling():
    x = RealVariable("x")
    y = RealVariable("y")

    variables = Variables2d(x, y)
    assert str(variables.x()) == str(x)
    assert str(variables.y()) == str(y)

    axes = Axes2d(-2.0, x, 2.0, -1.5, y, 1.5)
    figure = LabelledFigure(axes)

    box = LabelledRealBox(
        RealSpace([x, y]),
        RealBox([
            RealInterval(Real(Decimal("-1.0")), Real(Decimal("1.0"))),
            RealInterval(Real(Decimal("-0.5")), Real(Decimal("0.5"))),
        ]),
    )

    assert box.dimension() == 2
    figure.draw(box)
