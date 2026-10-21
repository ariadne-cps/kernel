/***************************************************************************
 *            test_point.cpp
 *
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 ****************************************************************************/

#include "geometry/point.hpp"
#include "utility/string.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

// These explicit specializations provide the stable human-readable names
// used by diagnostics and other type-aware clients.
Void test_point_class_names() {
    ARIADNE_TEST_EQUALS(class_name<Point<Real>>(), String("Point<Real>"));
    ARIADNE_TEST_EQUALS(class_name<Point<FloatDP>>(), String("Point<FloatDP>"));
    ARIADNE_TEST_EQUALS(class_name<Point<FloatDPBounds>>(), String("Point<FloatDPBounds>"));
    ARIADNE_TEST_EQUALS(class_name<Point<FloatDPApproximation>>(), String("Point<FloatDPApproximation>"));
}

Int main() {
    ARIADNE_TEST_CALL(test_point_class_names());
    return ARIADNE_TEST_FAILURES;
}
