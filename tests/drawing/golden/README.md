# Drawing golden images

These 11 PNGs were supplied from a successful local build of `DefaultGraphicsBackend`
on 16 October 2026, using a 320 x 240 drawing area (490 x 290 pixels including
plot margins). They are versioned, read-only pixel-exact regression references
for `test_default_graphics_regression`.

The C++ test generates the corresponding images again and uses the existing
self-contained PNG reader on both generated and reference files. All RGB
pixels are compared exactly, including labels, backgrounds, margins, fill and
stroke. A mismatch fails the test and writes `*-golden-diff.ppm` in the
build's drawing-test directory.

The golden images use the stored-block PNG encoding emitted by the renderer.
Tests require no extra image library or Python scripts, and never overwrite
the golden files. For an intentional rasterization change, inspect the
difference images and verify semantic correctness before **explicitly**
replacing and reviewing the references in Git.
