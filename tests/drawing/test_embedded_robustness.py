#!/usr/bin/env python3
"""Check robust embedded-only drawing without a reference renderer."""
import sys
from png_reader import read_png


def main(path):
    width, height, pixels = read_png(path)
    assert (width, height) == (234, 98), (width, height)

    def pixel(x, y):
        # Plot: x in [160,224), y in [10,58); world axes [-1,1].
        px = round(160 + (x + 1) * 32)
        py = round(10 + (1 - y) * 24)
        return pixels[py * width + px]

    background = (247, 247, 247)
    assert pixel(-0.73, 0.0) == background, "zero-opacity fill changed background"
    assert pixel(-0.25, 0.0) == (124, 251, 124), "half-opacity green fill incorrect"
    assert pixel(0.30, 0.0) == (255, 0, 0), "opaque red fill incorrect"
    assert pixel(0.75, 0.0) == (255, 255, 0), "clamped opacity incorrect"
    assert pixel(-0.9, -0.9) == (0, 0, 255), "clipped blue rectangle incorrect"
    assert pixels[25 * width + 150] == (255, 255, 255), "plot clipping changed margin"
    assert pixel(0.0, 0.0) == background, "invalid geometry changed background"
    print("Embedded robustness pixel checks passed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_embedded_robustness.py embedded-robustness.png")
    main(sys.argv[1])
