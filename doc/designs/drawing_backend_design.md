# Embedded drawing backend design

This document records the implementation and engineering decisions of the
self-contained 2D drawing backend in Ariadne Kernel. It describes the **current
implementation** rather than a proposed API. Keep it synchronized with changes
to the rasterizer, output format, tests, and performance benchmark.

## 1. Goals and design principles

The backend provides a compact, portable, CPU-based way to render scientific
figures to PNG without depending on an external graphics, font, or compression
library at runtime.

The main design principles are:

1. **Stable drawing abstraction.** Keep rasterization details behind the
   existing `CanvasInterface`, so figure objects and geometric algorithms do
   not manipulate pixels directly.
2. **Device-space semantics.** Transform coordinates on path insertion;
   measure stroke widths and antialiasing in output pixels.
3. **Correct compound paths.** Compute fill coverage across all subpaths
   together, with a nonzero winding rule and support for holes.
4. **Coverage before colour.** Union the coverage contributed by one drawing
   operation, then blend its colour only once per pixel.
5. **Specialize common workloads.** Retain a general polygon path but use
   simpler rasterizers for rectangles and finite-width line segments.
6. **Deterministic, inspectable implementation.** Prefer explicit numerical
   geometry and standard-library-only output encoding; benchmark the rendering
   code independently from upstream geometric computations.

## 2. Components and ownership

| File or interface | Responsibility |
| --- | --- |
| [`include/drawing/embedded.hpp`](../../include/drawing/embedded.hpp) | `EmbeddedCanvas`, `EmbeddedGraphicsBackend`, path/paint types and raster state. |
| [`src/drawing/embedded.cpp`](../../src/drawing/embedded.cpp) | Coordinate transformation, rasterization, blending, axis labels and PNG writer. |
| [`include/drawing/graphics_interface.hpp`](../../include/drawing/graphics_interface.hpp) | Backend-independent `CanvasInterface` contract. |
| [`include/drawing/graphics_base.hpp`](../../include/drawing/graphics_base.hpp) | Shared `CanvasBase::fill_boundary()` implementation. |
| [`src/drawing/graphics_manager.cpp`](../../src/drawing/graphics_manager.cpp) | Registration of the embedded backend as the default. |
| [`benchmarks/drawing/`](../../benchmarks/drawing/) | Independent absolute-performance measurements. |
| [`tests/drawing/`](../../tests/drawing/) | Image regression fixtures and label-legibility checks. |

`EmbeddedGraphicsBackend::make_canvas()` creates a canvas of the requested
drawing-area size. The `EmbeddedCanvas` owns the raster image and all temporary
path and coverage buffers. Higher-level `Figure` and `Drawable2dInterface`
objects access it only through the canvas abstraction.

## 3. Supported drawing functionality

- **Canvas initialization:** configurable x/y ranges and axis names, numeric
  endpoint labels, opaque white margins, light-coloured plot background and
  rectangular clipping to the plot.
- **Paths:** `move_to()`, `line_to()`, multiple subpaths, implicit closing
  edges for filling, and full-circle arcs through `circle()`.
- **Points:** `dot()` as a small circular path, with configurable dot radius.
- **Stroke:** configurable RGB colour and line width in device pixels;
  butt-ended segments and miter joins with a limit of 10, falling back to a
  bevel at sharp turns.
- **Fill:** configurable RGB colour and opacity; nonzero winding across the
  entire path, including oppositely oriented nested contours.
- **Boundaries:** `fill_boundary()` via the shared canvas base class.
- **Layering:** successive drawing operations are composited into the same
  opaque image in submission order.
- **Output:** `write()` produces an RGBA8 PNG file; `finalise()` explicitly
  draws the outer plot frame.

The API semantics are important: `stroke()` consumes the current path.
`fill()` fills the complete path, **then strokes its outline**, and finally
consumes it. A circle is appended to the current path rather than always
creating an independent subpath; a preceding current point is connected to
the circle's starting point.

Line and fill colours are separate. Fill opacity does not change stroke
opacity.

## 4. Canvas lifecycle and coordinate mapping

1. **Construction.** Allocate an RGBA byte buffer and a 16-bit coverage
   buffer. The raster includes a left margin of 160 px, right margin of
   10 px, top margin of 10 px and bottom margin of 40 px in addition to the
   requested plotting area.
2. **`initialise()`.** Clear pixels to opaque white; render axis labels;
   paint the plotting rectangle; compute the world-to-device affine transform;
   enable clipping for subsequent geometry.
3. **Draw.** Transform input points on insertion, accumulate subpaths, and
   rasterize/composite at `fill()` or `stroke()`.
4. **`finalise()` (optional).** Disable plotting-area clipping temporarily,
   draw the plot border and restore the drawing settings. It is not
   automatically invoked by the current `Figure::write()` code path.
5. **`write()`.** Encode the current pixels to PNG. If the supplied filename
   has no dot, append `.png`.

The transformation uses `x_device = tx + sx*x` and
`y_device = ty + sy*y`. The y scale is negative because image rows grow
downwards. Supplied axis bounds must be finite and strictly increasing.

For robustness, construction rejects zero drawing dimensions, drawing
dimensions above 16,384 in either axis, or a total allocated raster above
100 million pixels (including margins). Invalid non-finite drawing inputs
are ignored by relevant path operations.

## 5. Raster pipeline

### 5.1 Paths and subpixel masks

The implementation stores transformed device-space points in a vector of
subpaths. Each output pixel has **16 coverage bits**, representing a fixed
4 x 4 grid of subpixel sample positions.

All geometric fragments belonging to one operation OR their covered sample
bits into the same mask. The rasterizer records indices of newly touched
pixels in a dirty-pixel list. At the end of the operation, only those pixels
are visited for compositing; their masks are reset for reuse.

For a pixel mask `m`, the fractional coverage is:

```text
coverage = popcount(m) / 16
alpha    = clamp(fill_opacity, 0, 1) * coverage
result   = (1 - alpha) * existing_rgb + alpha * source_rgb
```

The stroke uses opaque source alpha and fractional edge coverage. Colour
channels are quantized to 8 bits after blending. For repeated compositions with the same paint, the backend caches channel
lookup tables indexed by coverage count and existing channel value,
preserving the same quantization. A four-entry least-recently-used cache
avoids recreating the tables when fill and outline colours alternate. This is an opaque RGB
composition model stored in an RGBA buffer with alpha 255, not a general
transparency-preserving image compositor.

**Why coverage is merged first:** drawing overlapping edges or subpaths
should not apply the same semi-transparent paint repeatedly within a
single `fill()` or `stroke()` call.

### 5.2 General polygon fill: active-edge scanlines

The generic fill rasterizer:

1. Converts path contours to non-horizontal, directed edges and assigns a
   signed winding contribution to each edge.
2. Sorts edges by their minimum y coordinate.
3. Advances through subpixel scanlines, maintaining only active edges.
4. Computes and sorts active x crossings; accumulates the winding number.
5. Marks horizontal intervals with nonzero winding directly in coverage
   masks.

Vertical edge activation uses half-open ranges to avoid double-counting
shared vertices. Filling is evaluated **across all subpaths together**;
oppositely oriented inner contours can therefore create genuine holes.

The active-edge approach avoids checking every edge for every subpixel,
which is particularly important for plots containing many polygons.

### 5.3 Fast path: axis-aligned rectangles

A contour with exactly four distinct axis-aligned corners (optionally
followed by its starting point) is recognized as a rectangle. The backend
validates both edge directions and all four unique corners before taking
this path.

The specialized fill supports a single rectangle and **two strictly nested
rectangles**. For a nested pair, contours with opposite winding leave a hole;
contours with the same winding leave the entire outer rectangle filled. The
order in which the contours were submitted does not matter. Touching,
overlapping or disjoint rectangles, as well as arbitrary quadrilaterals,
concave polygons and other compound paths, retain the general scanline
algorithm.

The rectangle renderer clamps its bounds to the plot region and generates
horizontal and vertical 4x4 sample masks directly, including fractional
pixel boundaries. For a nested pair it subtracts the inner subpixel mask
only when the windings oppose. Since the complete per-pixel coverage is
known immediately, it looks up the resulting RGB values in the paint cache
and **composites directly**, without constructing dirty-pixel lists or a
temporary coverage buffer. This avoids sorting edges, managing active-edge
lists and traversing each rectangle's covered pixels a second time.

### 5.4 Stroke: half-plane intersections

A finite-width butt-ended segment from `a` to `b` can be represented by
two pairs of linear inequalities. If `d=b-a`, its covered region satisfies:

```text
0 <= dot(p-a, d) <= dot(d, d)
abs(cross(p-a, d)) <= (line_width / 2) * length(d)
```

For fixed subpixel y, these inequalities constrain x to an interval.
The implementation intersects those bounds **once per subpixel scanline**
and ORs the corresponding range of coverage bits into the mask. This
replaces 16 independent distance/projection tests for each candidate pixel.

For **exactly horizontal or vertical segments**, the stroke is a rectangular
region in device coordinates. A separate fast path computes horizontal and
vertical 4x4 sample masks once per pixel row, including fractional boundaries
and plot clipping, then merges the masks into the existing union buffer.
This avoids the per-subpixel half-plane intersections for the four sides of
axis-aligned rectangles. It preserves butt-end caps and leaves oblique
segments on the general interval algorithm.

Joining strokes use geometric triangles to fill the inner wedge and extend
the outside to a miter tip, subject to the miter-length limit. Sharp joins
fall back to a bevel. Small triangles and round dots still use localized
sample-based rasterization.

### 5.5 Circles and tessellation

`circle()` transforms the centre and radius to device space, then builds
a full polygonal arc; nonuniform x/y scales yield elliptical contours.

The number of segments is selected from a **0.1 device-pixel chord-sagitta
target**, with a minimum of 12 and maximum of 16,384 segments. The cap
means 0.1 px is a design target, not a universal bound for arbitrarily
large radii. Both filled circles and stroked circular paths reuse the
general path rasterizers.

### 5.6 Axis-label font

The backend incorporates a limited 5 x 7 character skeleton. It renders
short horizontal, vertical and selected diagonal connections between lit
cells, plus small rounded endpoint disks, using the same antialiasing mask
as ordinary geometry.

This keeps label drawing self-contained and portable. The present font is
a simple fixed-size internal font: it does **not** provide Unicode shaping,
font discovery, kerning, scalable typography or general text layout.

## 6. PNG writer

The encoder is implemented in the rasterizer translation unit and writes:

- the PNG signature plus `IHDR`, `IDAT`, and `IEND` chunks;
- 8-bit-per-channel RGBA image rows with PNG filter type 0;
- a zlib stream using **stored (uncompressed) DEFLATE blocks**;
- Adler-32 on the data stream and CRC-32 on PNG chunks.

Using stored blocks removes the need for a compression dependency and keeps
the implementation compact. It trades off compressed-file size and some
encoding efficiency. PNG encoding and filesystem output are distinct from
path rasterization and can be measured separately.

## 7. Performance choices and measurement

Implemented performance techniques include:

1. Clipping and bounding-region checks to avoid unnecessary off-screen work.
2. Active-edge scanline filling for arbitrary compound polygons.
3. Direct composition of single rectangles and strictly nested rectangular
   frames, with correct nonzero-winding semantics for holes.
4. Interval-based stroke rasterization instead of per-sample distance tests,
   and a row-mask fast path for horizontal and vertical strokes.
5. Per-pixel bit masks and dirty-pixel tracking for single-pass compositing.
6. Device-space, curvature-aware arc tessellation with a segment-count cap.
7. Direct opaque background writes instead of per-pixel floating-point blending.
8. Four-entry least-recently-used paint cache with lookup tables for the
   17 possible 4x4 coverage counts and 256 destination values per RGB
   channel. Tables are only rebuilt for a cache miss, avoiding the cost
   of alternating fill/stroke paint. The cache belongs to each canvas.

The executable [`benchmark_drawing`](../../benchmarks/drawing/benchmark_drawing.cpp)
under `benchmarks/drawing` is the **absolute-performance reference**. It
drives `CanvasInterface` directly with deterministic command buffers whose
geometry is prepared *before* measurements. No `Figure` construction,
constrained-set approximation or other geometric computation is included
in the timed rendering interval.

Its scenarios exercise rectangles, polylines, circles, concave polygons,
compound nonzero-winding paths and a mixed workload. An aligned terminal
report separates canvas construction, initialization (including labels), pure
rasterization, and optional PNG encoding/file I/O. The benchmark uses the
Kernel command-line profile extended with benchmark-specific selectors,
and the standard logger for its report and optional progress messages. Rendering measurements
include median and p95 times; throughput is reported in primitives/second.
An optional `--breakdown` diagnostic runs full, geometry-only, fill-only
and stroke-only passes with identical precomputed commands and selective
paint suppression. Each pass is independently timed: results are not
exclusive phase durations and must not be summed.

Example Release run:

```sh
cmake --build build --target benchmark_drawing --parallel
./build/benchmarks/drawing/benchmark_drawing --backend embedded --iterations 10 --warmup 2
```

The benchmark is independent of other backend implementations; its normal
operation and absolute timing metrics remain valid as the drawing module
evolves. See [benchmark documentation](../../benchmarks/drawing/README.md)
for available scenarios and options.

Avoid recording a single machine's benchmark numbers as permanent design
properties. Compare results only with the same resolution, command mix,
compiler, optimization configuration and hardware.

## 8. Verification and maintenance constraints

Drawing regressions should preserve both **geometric semantics** and
**performance**:

- The drawing fixtures cover compound paths with holes, joins, clipping,
  circular contours, fractional-coordinate rectangles, variable-width
  strokes and concave polygons. Two additional scenes exercise nested
  rectangular paths of both winding directions and submission orders,
  including touching contours that must use general filling.
  A separate fixture covers horizontal/vertical strokes with fractional
  widths and coordinates, short segments and clipped end caps.
- A standalone embedded-only regression checks interior pixels of the
  nested-rectangle scenes; it does not depend on a reference backend.
- The embedded label regression test checks that axis text remains visible
  and connected.
- Image-comparison tests use a CTest fixture. The `test_backend_comparison`
  executable produces the PNG inputs in the drawing test's binary directory
  (normally `build/tests/drawing/`), and CTest runs it automatically when
  a dependent comparison or label test is selected. The executable must
  already have been built; CTest does not compile targets.
- PNG checks and image comparisons must evaluate the plotting area
  separately from the axis-label regions; large uniform backgrounds must
  not conceal localized errors.
- The benchmark's scenario definitions should remain stable between
  optimization iterations to allow useful comparisons.
- Optimizations must preserve the coverage-mask and single-compositing
  semantics, including on fractional pixel boundaries and clipped paths.

Important current limitations and trade-offs:

- Fixed 4x4 sampling produces quantized fractional coverage.
- Curves are tessellated rather than analytically rasterized.
- The built-in font is deliberately small and non-general-purpose.
- The PNG stream is valid but intentionally uncompressed.
- The interface targets static PNG output; animated drawing is not
  implemented by this canvas.
- The raster buffers are allocated per canvas; there is no GPU acceleration
  or incremental scene graph.
- Outer-frame rendering requires an explicit `finalise()` call.

When changing the implementation, update the relevant section here, add
focused visual/regression coverage, and re-run the isolated Release
benchmark before claiming a performance improvement.
