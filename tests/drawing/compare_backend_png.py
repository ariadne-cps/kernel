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
    parser.add_argument("--max-mean-error", type=float, default=0.10)
    parser.add_argument("--diff", help="Optional diagnostic RGB difference image (PPM)")
    args = parser.parse_args()
    a = read_png(args.cairo)
    b = read_png(args.embedded)
    assert a[:2] == b[:2], f"Canvas size differs: {a[:2]} != {b[:2]}"
    count = len(a[2])
    mean_error = sum(abs(i-j) for pa,pb in zip(a[2],b[2]) for i,j in zip(pa,pb)) / (count*3*255)
    differing = sum(any(abs(i-j)>32 for i,j in zip(pa,pb)) for pa,pb in zip(a[2],b[2]))
    print(f"Cairo vs Embedded {a[0]}x{a[1]}: mean normalized RGB error {mean_error:.5f};"
          f" significant differing pixels {differing}/{count} ({100*differing/count:.2f}%)")
    if args.diff:
        # Magnify RGB difference for visual diagnosis without external packages.
        with open(args.diff, "wb") as f:
            f.write(f"P6\\n{a[0]} {a[1]}\\n255\\n".encode("ascii"))
            f.write(bytes(min(255,abs(i-j)*4) for pa,pb in zip(a[2],b[2]) for i,j in zip(pa,pb)))
    assert mean_error <= args.max_mean_error, f"Difference exceeds {args.max_mean_error}"


if __name__ == "__main__":
    main()
