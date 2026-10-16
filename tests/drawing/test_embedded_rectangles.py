#!/usr/bin/env python3
"""Guard nonzero-winding semantics for nested rectangular contours.

The fixtures contain rectangles with equal/opposite orientation, reversed
submission order, fractional coordinates and a contour touching its parent.
Checks require only the embedded renderer and the standard-library PNG reader.
"""
import sys

from png_reader import read_png


def get_pixel(image, x, y):
    width, height, pixels = image
    # The 320 x 240 test plot is offset by (160, 10); the y axis is inverted.
    plot_width = width - 170
    plot_height = height - 50
    px = round(160 + (x + 1) * plot_width / 2)
    py = round(10 + (1 - y) * plot_height / 2)
    return pixels[py * width + px]


def expect(image, x, y, filled):
    rgb = get_pixel(image, x, y)
    # The known fill colour has far less green/blue than the near-white plot.
    actual_filled = rgb[1] < 210 and rgb[2] < 210
    assert actual_filled == filled, (
        f"wrong rectangular winding at ({x}, {y}): "
        f"pixel={rgb}, expected {'filled' if filled else 'hole'}"
    )


def main(hole_scene, mixed_scene):
    holes = read_png(hole_scene)
    mixed = read_png(mixed_scene)
    assert holes[:2] == mixed[:2] == (490, 290)

    # Every nested contour has opposite orientation, irrespective of the
    # order in which the inner and outer contour were appended.
    for x in (-0.66, -0.03, 0.62):
        expect(holes, x, 0.0, filled=False)
    expect(holes, -0.91, 0.0, filled=True)
    expect(holes, -0.28, 0.0, filled=True)

    # A same-direction contour is not a hole. An opposite-direction one is.
    expect(mixed, -0.60, 0.0, filled=True)
    expect(mixed, 0.64, 0.0, filled=False)
    expect(mixed, 0.94, 0.0, filled=True)

    # An inner rectangle sharing a boundary must take the general winding
    # path; it still subtracts the inner region.
    expect(mixed, 0.0, -0.82, filled=False)
    expect(mixed, -0.45, -0.82, filled=True)
    print("Embedded nested-rectangle winding checks passed")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: test_embedded_rectangles.py scene7-embedded.png scene8-embedded.png"
        )
    main(sys.argv[1], sys.argv[2])
