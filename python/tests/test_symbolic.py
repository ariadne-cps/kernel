from pyariadne import *


def test_symbolic():
    x = RealVariable("x")
    y = RealVariable("y")
    expression = x * x + y

    derivative(expression, x)
    make_function([x, y], expression)
    str((0 <= x) & (x <= 1))
