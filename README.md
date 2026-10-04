# Ariadne Kernel

Ariadne Kernel is the standalone C++20 kernel layer extracted from Ariadne. It contains the Function, Geometry, IO, Symbolic and Solving modules and consumes the standalone Algebra and Threading repositories as direct dependencies.

## Layout

Public headers are installed from:

- `include/function`
- `include/geometry`
- `include/io`
- `include/symbolic`
- `include/solving`

Implementations live under the matching directories in `src/`.

The umbrella header is `include/ariadne-kernel.hpp`.

## Dependencies

The repository has three Git submodules:

- `submodules/configuration`
- `submodules/algebra`
- `submodules/threading`

CMake verifies that the direct Configuration revisions agree and that the Utility revision reached through Algebra matches the Utility revision reached through Threading.

The pinned revisions used for the extraction are:

- Configuration: `5286138155b55b57392d532adbd742dd542899bb`
- Algebra: `83966f292749152082b7d200dca36c31e496b183`
- Threading: `d3473103822e91ea7cde65a1226c93c1dd9b9347`
- Utility on both dependency paths: `9194f4dd7c6a89fba422bed81382d4135cddc004`

Cairo is a required external dependency of the Kernel IO graphics backend.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/kernel.git
cd kernel
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target everything
ctest --test-dir build --output-on-failure
```

## Python bindings

When a compatible Python development environment is available, the standalone build produces `pyariadne` from the lower Algebra bindings plus the Kernel bindings for Function, Geometry, IO, Solving and Symbolic. Sweeper bindings are supplied by Algebra; IO exports the CLI and graphics facilities but not Logging.

Dynamics, Evolution and Hybrid bindings are intentionally not part of this repository.

## Benchmarks

Function and Solving benchmarks are included. The full Barr3 binary benchmark is enabled only when `benchmarks/solving/data/smt_barr3_full64.bin` is present; the source/header-based regression test is included independently.

## Extraction notes

The obsolete `function/calculus_base.hpp` header was not imported because it references a missing `calculus_interface.hpp` and belongs to the legacy dynamical-calculus surface rather than the standalone kernel boundary.

Two stale forward declarations of Hybrid types were removed from Geometry during extraction.
