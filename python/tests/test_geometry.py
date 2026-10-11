from pyariadne import *


def test_geometry():
    dp = DoublePrecision()
    x = FloatDP(1, dp)
    box = FloatDPExactBox([(-x, x), (-x, x)])

    assert Box[FloatDP] == FloatDPExactBox
    assert box.dimension() == 2

    left, right = box.split(0)
    assert left.dimension() == 2
    assert right.dimension() == 2
