#include "macro/testing.hpp"
#include "numeric/numeric.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "symbolic/expression_set.hpp"
#include "symbolic/space.hpp"
#include "labelling/set.hpp"

using namespace Ariadne;

class TestLabelledSet {
  private:
    RealVariable x{"x"}, y{"y"}, z{"z"};

  public:
    Void test() {
        ARIADNE_TEST_CALL(test_set());
        ARIADNE_TEST_CALL(test_point_adapter());
        ARIADNE_TEST_CALL(test_box_adapter());
    }

    Void test_set() const {
        LabelledSet<RealBox> set(RealSpace({x,y,z}),RealBox({{0.1_dec,0.5_dec},{-2.3_dec,1},{0,3}}));
        ARIADNE_TEST_EQUALS(set.dimension(),3);
        ARIADNE_TEST_EQUAL(set.space(),RealSpace({x,y,z}));
        ARIADNE_TEST_EQUAL(set.variables(),RealSpace({x,y,z}).variables());
        RealVariable u("u");
        ARIADNE_TEST_EXECUTE(set.euclidean_set(RealSpace({y,z})));
        ARIADNE_TEST_FAIL(set.euclidean_set(RealSpace({y,u})));
    }

    Void test_point_adapter() const {
        Map<RealVariable,Real> bindings;
        bindings.insert(x,0.1_dec);
        bindings.insert(y,-1);
        VariablesPoint<Real> point(bindings);
        LabelledSet<RealPoint> labelled(point);
        ARIADNE_TEST_EQUALS(labelled.dimension(),2);
    }

    Void test_box_adapter() const {
        RealVariablesBox box(RealSpace({x,y}),RealBox({{0.1_dec,0.5_dec},{-2.3_dec,1}}));
        LabelledSet<RealBox> labelled(box);
        ARIADNE_TEST_EQUALS(labelled.dimension(),2);
        ARIADNE_TEST_EQUAL(labelled[x],box[x]);
    }
};

Int main() {
    TestLabelledSet().test();
    return ARIADNE_TEST_FAILURES;
}
