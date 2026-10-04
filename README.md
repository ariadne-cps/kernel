# Ariadne Kernel

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/kernel/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/kernel/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/kernel/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/kernel/actions/workflows/win.yml)
[![Coverage Status](https://github.com/ariadne-cps/kernel/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/kernel/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/kernel/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/kernel)

Ariadne Kernel is the standalone C++20 kernel layer of Ariadne. It contains the Function, Geometry, IO, Symbolic and Solving modules and consumes the standalone Algebra and Threading repositories as direct dependencies.

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

Cairo is a required external dependency of Ariadne Kernel.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/kernel.git
cd kernel
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target everything
ctest --test-dir build --output-on-failure
```

## Python bindings

When a compatible Python development environment is available, the standalone build produces `pyariadne` from the lower Algebra bindings plus the Kernel bindings for Function, Geometry, IO, Solving and Symbolic. Sweeper bindings are supplied by Algebra, while IO provides the CLI and graphics facilities.

## Contribution guidelines

If you would like to contribute to Ariadne Kernel, please contact the developer:

- Luca Geretti <luca.geretti@univr.it>

## License

Ariadne Kernel is released under the GNU General Public License v3.0.
