#!/usr/bin/env python3
"""Regression guard for the axis-label stroke rendering in EmbeddedCanvas.

The previous renderer drew only horizontal glyph runs; 'x', '1', and other
characters disappeared into isolated dots although geometry comparisons passed.
"""
import sys
from png_reader import read_png


def ink(pixels, width, box):
    x0, y0, x1, y1 = box
    return sum(all(c < 90 for c in pixels[y * width + x])
               for y in range(y0, y1) for x in range(x0, x1))


def main(path):
    width, height, pixels = read_png(path)
    plot_width = width - 160 - 10
    plot_height = height - 10 - 40
    assert plot_width > 0 and plot_height > 0

    # The fixture draws x/y axis names and values -1 and 1. These regions
    # contain text only; avoid plot anti-aliasing and drawing geometry.
    x_center = 160 + plot_width // 2
    regions = {
        "horizontal axis name": ((x_center - 25, height - 37,
                                  x_center + 15, height - 2), 45),
        "vertical axis name": ((125, 10 + plot_height // 2 - 4,
                                159, 10 + plot_height // 2 + 36), 65),
        "upper numeric label": ((127, 8, 158, 45), 45),
    }
    for name, (box, minimum) in regions.items():
        count = ink(pixels, width, box)
        print(f"{name}: {count} dark pixels, minimum {minimum}")
        assert count >= minimum, f"Embedded text regression: {name}"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_embedded_labels.py comparison-embedded.png")
    main(sys.argv[1])
