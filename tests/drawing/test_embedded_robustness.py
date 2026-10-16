#!/usr/bin/env python3
"""Embedded-only pixel and clipping checks for robustness rendering."""
import sys
from png_reader import read_png


def main(path):
    width, height, pixels = read_png(path)
    assert (width, height) == (234, 98)

    def at(x, y):
        px = round(160 + (x + 1) * 32)
        py = round(10 + (1 - y) * 24)
        return pixels[py * width + px]

    def close(actual, expected, label, tolerance=2):
        assert all(abs(a - b) <= tolerance for a, b in zip(actual, expected)), (
            f"{label}: found {actual}, expected {expected}"
        )

    # Background, transparent fill, 50% alpha, opaque fill and clamped alpha.
    background = (247, 247, 247)
    close(at(-0.7, 0.0), background, "transparent rectangle")
    close(at(-0.25, 0.0), (124, 251, 124), "half-opacity green")
    close(at(0.3, 0.0), (255, 0, 0), "opaque red")
    close(at(0.75, 0.0), (255, 255, 0), "clamped opacity")
    close(at(-0.55, 0.0), background, "untouched plot background")
    close(at(-0.9, -0.9), background, "offscreen rectangle outside plot")
    close(at(-0.9, -0.87), (0, 0, 255), "clipped blue rectangle")

    # The fixed white area outside the plot is never painted by geometry.
    close(pixels[60 * width + 100], (255, 255, 255), "left margin")
    close(pixels[30 * width + 230], (255, 255, 255), "right margin")
    print("Embedded independent opacity and clipping checks passed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_embedded_robustness.py embedded-robustness.png")
    main(sys.argv[1])
