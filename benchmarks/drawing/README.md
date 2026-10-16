# Drawing backend benchmark

`benchmark_drawing` measures **absolute rendering performance** independently
of the Geometry and Function modules. It calls `CanvasInterface` directly,
using deterministically precomputed drawing commands, rather than constructing
`Figure` objects or constrained image sets.

**Embedded is the default and only required backend.** Cairo is accepted as
an optional transition-time baseline when the library was built with
`ARIADNE_WITH_CAIRO=ON`. Removing Cairo later requires no changes to the
benchmark's default invocation or the build target.

## Build and run

Build a Release configuration, since Debug benchmarks are misleading:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DARIADNE_WITH_CAIRO=OFF
cmake --build build --target benchmark_drawing --parallel
./build/benchmarks/drawing/benchmark_drawing --iterations 7 --warmup 2
```

The executable's location depends on the generator and build configuration.
The existing project-wide `benchmarks` target also builds this benchmark.

With Cairo available, compare identical commands without changing the
absolute numbers reported:

```sh
./build/benchmarks/drawing/benchmark_drawing --backend both --iterations 7 --warmup 2
```

Select one scenario/size (for tracking performance across commits):

```sh
./build/benchmarks/drawing/benchmark_drawing --scene rectangles --width 800 --height 800 --iterations 10 --warmup 2
```

Optionally measure PNG generation and file I/O **separately**, by writing to
`benchmark-drawing-<backend>-<scene>.png` in the current directory:

```sh
./build/benchmarks/drawing/benchmark_drawing --write-png --png-prefix /tmp/drawing-bench
```

## Scenarios

| Scenario | Primitives | Purpose |
| --- | ---: | --- |
| rectangles | 384 | Dense paving-like fills and narrow strokes |
| polylines | 72 | Long polylines, antialiasing and joins |
| circles | 120 | Curved filled paths, tessellation |
| polygons | 80 | Concave multi-edge paths |
| winding | 54 | Multiple subpaths and holes |
| mixed | 710 | Combination of all preceding scenes |

The scene command buffers, including coordinates and trigonometric values,
are prepared **outside the measured section** and are identical for Cairo
and Embedded. Every timed repetition creates a fresh canvas so no previous
frame can affect the next. Measurements use `steady_clock` and report
median and nearest-rank p95 over the measured repetitions.

Results are presented as an **aligned terminal table**, in the style of
the integration benchmarks. The report is printed through Ariadne's logging
system (default verbosity 1); no CSV files or CSV output are generated.

The benchmark extends `make_kernel_command_line_interface()` with a
`drawing-benchmark` contribution, so the usual logging options
(`--verbosity`, `--theme`, etc.) are available alongside the local
`--backend`, `--scene`, `--width`, `--height`, `--iterations`,
`--warmup`, `--write-png` and `--png-prefix` arguments.

The columns contain creation and initialization median times, pure drawing
median and nearest-rank p95 times, and shapes per second. When requested,
the report also includes PNG encoding and file-write time. Verbosity 2
logs each scene before measurement; verbosity 3 additionally logs individual
samples. Logging is outside all measured time intervals.

These measurements are not a correctness test. Continue to use the tests
under `tests/drawing` for Cairo/Embedded image comparisons and regression
checks. Compare measurements only with the same build type, machine,
compiler, size, scenario and iteration settings.
