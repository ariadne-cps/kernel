#!/usr/bin/env python3
"""Compare Cairo and Embedded PNG output using only the Python standard library.

Both images must have identical dimensions and approximate pixel parity. This
is deliberately tolerant of font rasterization and antialiasing differences.
"""
import argparse
import sys

from png_reader import read_png


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cairo")
    parser.add_argument("embedded")
    parser.add_argument("--max-plot-error", type=float, default=0.015)
    parser.add_argument("--max-plot-outliers", type=float, default=0.035)
    parser.add_argument("--diff", help="Diagnostic difference image (PPM)")
    args = parser.parse_args()
    a = read_png(args.cairo)
    b = read_png(args.embedded)
    assert a[:2] == b[:2], f"Canvas dimensions differ: {a[:2]} != {b[:2]}"
    width,height = a[:2]
    # Match the plot rectangle, excluding independently rasterized axis labels.
    # See LEFT/TOP/RIGHT/BOTTOM in src/drawing/{cairo,embedded}.cpp.
    x0,y0,x1,y1 = 160,10,width-10,height-40
    assert x1>x0 and y1>y0
    total=0
    error=0
    outliers=0
    max_difference=0
    for y in range(y0,y1):
        for x in range(x0,x1):
            index=y*width+x
            diffs=[abs(u-v) for u,v in zip(a[2][index],b[2][index])]
            error+=sum(diffs)
            outliers+=int(max(diffs)>32)
            max_difference=max(max_difference,max(diffs))
            total+=1
    mean=error/(3*255*total)
    fraction=outliers/total
    print(f"plot: mean RGB error={mean:.6f}, >32 outliers={fraction:.6f}, "
          f"max channel difference={max_difference}; plot pixels={total}")
    if args.diff:
        with open(args.diff,"wb") as file:
            file.write(f"P6\n{width} {height}\n255\n".encode("ascii"))
            file.write(bytes(min(255,abs(u-v)*4)
                for pa,pb in zip(a[2],b[2]) for u,v in zip(pa,pb)))
    assert mean<=args.max_plot_error, f"Plot mean error {mean:.5f} > {args.max_plot_error}"
    assert fraction<=args.max_plot_outliers, (
        f"Plot outlier fraction {fraction:.5f} > {args.max_plot_outliers}")


if __name__ == "__main__":
    main()
