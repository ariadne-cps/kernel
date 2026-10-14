#ifndef ARIADNE_LABELLING_LABELLING_HPP
#define ARIADNE_LABELLING_LABELLING_HPP

#include "symbolic/space.hpp"

namespace Ariadne {

class Labelling {
  public:
    explicit Labelling(RealSpace const& space) : _space(space) { }

    RealSpace const& space() const { return _space; }
    Set<RealVariable> variables() const { return _space.variables(); }

  protected:
    RealSpace _space;
};

} // namespace Ariadne

#endif
