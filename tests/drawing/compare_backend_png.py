#!/usr/bin/env python3
"""Compare Cairo and Embedded PNG output using only the Python standard library.

Both images must have identical dimensions and approximate pixel parity. This
is deliberately tolerant of font rasterization and antialiasing differences.
"""
import argparse
import pathlib
import struct
import sys
import zlib


def read_png(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", f"{path}: invalid PNG signature"
    offset, compressed, width, height, channels = 8, bytearray(), None, None, None
    while offset < len(data):
        length = struct.unpack_from(">I", data, offset)[0]
        tag = data[offset + 4:offset + 8]
        payload = data[offset + 8:offset + 8 + length]
        expected = struct.unpack_from(">I", data, offset + 8 + length)[0]
        actual = zlib.crc32(tag + payload) & 0xffffffff
        assert expected == actual, f"{path}: CRC mismatch in {tag!r}"
        offset += 12 + length
        if tag == b"IHDR":
            width, height, depth, mode, compression, filtr, interlace = struct.unpack(">IIBBBBB", payload)
            assert depth == 8 and mode in (2, 6) and not (compression or filtr or interlace)
            channels = 3 if mode == 2 else 4
        elif tag == b"IDAT":
            compressed.extend(payload)
        elif tag == b"IEND":
            break
    assert width and height and channels
    raw = zlib.decompress(compressed)
    stride = width * channels
    rows = []
    offset, previous = 0, bytearray(stride)
    for _ in range(height):
        mode = raw[offset]
        row = bytearray(raw[offset + 1:offset + 1 + stride])
        assert len(row) == stride
        offset += 1 + stride
        for i in range(stride):
            a = row[i - channels] if i >= channels else 0
            b = previous[i]
            c = previous[i - channels] if i >= channels else 0
            if mode == 1:
                predictor = a
            elif mode == 2:
                predictor = b
            elif mode == 3:
                predictor = (a+b)//2
            elif mode == 4:
                p = a+b-c
                distances = (abs(p-a), abs(p-b), abs(p-c))
                predictor = (a, b, c)[distances.index(min(distances))]
            else:
                assert mode == 0, f"Unsupported PNG filter {mode}"
                predictor = 0
            row[i] = (row[i] + predictor) & 255
        previous = row
        rows.extend(row[i:i+3] for i in range(0, stride, channels))
    return width, height, rows


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
