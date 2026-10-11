from pyariadne import *


def test_io():
    point = Point2d(1.0, 2.0)
    assert point.x == 1.0
    assert point.y == 2.0

    projection = Projection2d(2, 0, 1)
    assert projection.n == 2
    assert projection.i == 0
    assert projection.j == 1

    properties = GraphicsProperties()
    properties.set_line_width(1.0)
    properties.set_fill_opacity(0.5)
