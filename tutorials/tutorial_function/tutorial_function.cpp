#include <iostream>

#include "ariadne-kernel.hpp"

using namespace Ariadne;

template<class T>
void show(const char* label, const T& value) {
    std::cout << label << ": " << value << "\n";
}

int main() {
    DoublePrecision dp;

    std::cout << "=== Effective functions ===\n";

    auto x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
    auto x1=EffectiveScalarMultivariateFunction::coordinate(2,1);
    auto f=sqrt(sqr(x0)+sqr(x1));
    auto g=EffectiveVectorMultivariateFunction({x0+x1,x0*x1});

    show("f",f);
    show("g",g);

    auto approximate_point=FloatDPApproximationVector({3,4},dp);
    auto bounded_point=FloatDPBoundsVector({3,4},dp);

    show("f(approximate_point)",f(approximate_point));
    show("f(bounded_point)",f(bounded_point));
    show("g(approximate_point)",g(approximate_point));
    show("g(bounded_point)",evaluate(g,bounded_point));

    auto df_dx0=derivative(f,0);
    show("derivative(f,0)",df_dx0);

    auto identity2=EffectiveVectorMultivariateFunction::identity(2);
    show("compose(f,identity2)",compose(f,identity2));
    show("join(f,g)",join(f,g));

    std::cout << "\n=== Taylor function models ===\n";

    auto domain=BoxDomainType({{-1,+1},{-1,+1}});
    SweeperDP sweeper(ThresholdSweeper<FloatDP>(dp,ApproximateDouble(1e-8)));

    auto tx=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(domain,0,sweeper);
    auto ty=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(domain,1,sweeper);
    auto tidentity=ValidatedVectorMultivariateTaylorFunctionModelDP::identity(domain,sweeper);

    auto tf=1+2*tx+3*tx*ty+sqr(tx);
    auto tv=ValidatedVectorMultivariateTaylorFunctionModelDP({tx+ty,tx*ty});

    show("domain",tf.domain());
    show("codomain",tf.codomain());
    show("range",tf.range());
    show("tf",tf);
    show("tv",tv);

    show("derivative(tf,0)",derivative(tf,0));
    show("antiderivative(tf,0)",antiderivative(tf,0));

    auto subdomain=BoxDomainType({{-0.5_x,0.5_x},{-0.25_x,0.25_x}});
    show("restriction(tf,subdomain)",restriction(tf,subdomain));

    auto outer_domain=tidentity.codomain();
    auto outer_function=EffectiveScalarMultivariateFunction::coordinate(2,0);
    auto outer_model=ValidatedScalarMultivariateTaylorFunctionModelDP(
        outer_domain,outer_function,sweeper);

    show("compose(f,tidentity)",compose(f,tidentity));
    show("compose(outer_model,tidentity)",compose(outer_model,tidentity));

    auto coarse=tx+FloatDPBounds(-1,1,dp);
    show("inconsistent(tx,coarse)",inconsistent(tx,coarse));
    show("refines(tx,coarse)",refines(tx,coarse));
    show("refinement(tx,coarse)",refinement(tx,coarse));

    return 0;
}
