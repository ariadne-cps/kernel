#include <iostream>
#include <string>

#include "ariadne-kernel.hpp"

using namespace Ariadne;
using std::cout;

namespace {

void section(const std::string& title) {
    cout << "\n== " << title << " ==\n";
}

template<class T>
void show(const std::string& label, const T& value) {
    cout << label << ": " << value << "\n";
}

} // namespace

int main() {
    section("Named variables and spaces");

    RealVariable x("x");
    RealVariable y("y");
    RealVariable z("z");
    RealVariables state("state",3);
    IntegerVariable mode("mode");
    BooleanVariable enabled("enabled");

    RealSpace xyz({x,y,z});

    show("space", xyz);
    show("space dimension", xyz.dimension());
    show("index of y", xyz.index(y));
    show("generated state variables", state);
    show("integer variable", mode);
    show("boolean variable", enabled);

    section("Building symbolic real expressions");

    // Symbolic expressions keep variable names rather than committing to a
    // coordinate order.  They can therefore be manipulated before deciding
    // how they should be embedded in Euclidean space.
    RealExpression expression =
        sqr(x) + Real(2)*y + sin(z);

    RealExpression redundant =
        expression + (x-x) + log(exp(y));
    RealExpression simplified = simplify(redundant);

    show("expression", expression);
    show("redundant expression", redundant);
    show("simplified expression", simplified);

    section("Differentiation and structural analysis");

    RealExpression dx = simplify(derivative(expression,x));
    RealExpression dy = simplify(derivative(expression,y));
    RealExpression dz = simplify(derivative(expression,z));

    show("d/dx", dx);
    show("d/dy", dy);
    show("d/dz", dz);

    RealExpression repeated =
        sqr(x+y) + sin(x+y) + (x+y);
    show("expression with repeated subexpressions", repeated);
    show("nodes before CSE", count_nodes(repeated));
    show("distinct nodes", count_distinct_nodes(repeated));
    show("distinct node pointers before CSE",
         count_distinct_node_pointers(repeated));

    eliminate_common_subexpressions(repeated);

    show("distinct node pointers after CSE",
         count_distinct_node_pointers(repeated));

    section("Substitution and valuation");

    List<RealAssignment> substitutions({
        let(y)=x+Real(1),
        let(z)=Real(0)
    });
    RealExpression substituted =
        substitute(expression,substitutions);

    show("substituted expression", substituted);

    // A Valuation gives concrete values to named variables without changing
    // the expression itself.
    RealValuation valuation;
    valuation.set(x,Real(Rational(1,2)));
    valuation.set(y,Real(1));
    valuation.set(z,Real(0));

    show("valuation", valuation);
    show("evaluated expression", evaluate(expression,valuation));

    section("Assignments describe algebraic, discrete and differential rules");

    RealAssignment algebraic_rule =
        let(x)=y+Real(1);
    PrimedRealAssignment discrete_update =
        prime(x)=x+y;
    DottedRealAssignment differential_rule =
        dot(x)=-x+y;

    show("algebraic assignment", algebraic_rule);
    show("discrete update", discrete_update);
    show("differential equation", differential_rule);

    section("From symbolic expressions to callable functions");

    EffectiveScalarMultivariateFunction scalar_function =
        make_function(xyz,expression);

    Vector<RealExpression> vector_expression({
        x+y,
        sqr(z)+x
    });
    EffectiveVectorMultivariateFunction vector_function =
        make_function(xyz,vector_expression);

    show("scalar function", scalar_function);
    show("vector function", vector_function);

    // RealSpace fixes the argument ordering.  The expression remains
    // name-based, but the resulting functions use coordinates [x,y,z].
    show("x coordinate index", xyz.index(x));
    show("z coordinate index", xyz.index(z));

    section("Symbolic boxes and predicates");

    RealVariableInterval x_range(Real(0),x,Real(1));
    RealVariableInterval y_range(Real(-1),y,Real(1));

    RealVariablesBox domain({x_range,y_range});

    RealExpressionConstraintSet constraints({
        x+y <= Real(Rational(3,2)),
        sqr(x)+sqr(y) <= Real(2)
    });

    show("symbolic box", domain);
    show("symbolic constraints", constraints);

    section("Combining bounds and constraints");

    RealExpressionBoundedConstraintSet feasible_region =
        intersection(domain,constraints);

    show("bounded symbolic set", feasible_region);

    RealSpace xy({x,y});
    BoundedConstraintSet euclidean_region =
        feasible_region.euclidean_set(xy);

    show("Euclidean bounded constraint set", euclidean_region);

    // approximate_euclidean_set builds a validated computational
    // representation from the symbolic description.  This is the bridge
    // from readable model specification to the geometry/solving modules.
    ValidatedConstrainedImageSet validated_region =
        approximate_euclidean_set(feasible_region,xy);

    show("validated approximation", validated_region);
    show("validated bounding box", validated_region.bounding_box());

    section("Choosing a symbolic abstraction");
    cout
        << "Variable / RealSpace: give names to quantities and define coordinate order.\n"
        << "Expression: build and manipulate formulas before numerical evaluation.\n"
        << "Valuation: evaluate named expressions at concrete values.\n"
        << "let / prime / dot assignments: describe algebraic, discrete and differential rules.\n"
        << "make_function: compile named expressions into coordinate-based functions.\n"
        << "RealVariablesBox and expression constraint sets: define labelled feasible regions.\n"
        << "approximate_euclidean_set: turn symbolic regions into validated computational sets.\n";

    return 0;
}
