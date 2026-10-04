#include <iostream>

#include "ariadne-kernel.hpp"
#include "io/cairo.hpp"

using namespace Ariadne;

int main() {
    DoublePrecision dp;
    SweeperDP sweeper(ThresholdSweeper<FloatDP>(dp,ApproximateDouble(1e-8)));

    auto domain=BoxDomainType({{-1,+1}});
    auto x=EffectiveScalarMultivariateFunction::coordinate(1,0);
    auto model=ValidatedScalarMultivariateTaylorFunctionModelDP(domain,x,sweeper);

    CairoGraphicsBackend backend;
    auto canvas=backend.make_canvas("kernel-install-check",16u,16u,false);

    std::cout << model << "\n";
    return canvas ? 0 : 1;
}
