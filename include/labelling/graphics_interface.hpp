#ifndef ARIADNE_LABELLING_GRAPHICS_INTERFACE_HPP
#define ARIADNE_LABELLING_GRAPHICS_INTERFACE_HPP

#include "io/graphics_interface.hpp"

namespace Ariadne {

struct Variables2d;

class LabelledDrawable2dInterface {
  public:
    virtual ~LabelledDrawable2dInterface() = default;
    virtual LabelledDrawable2dInterface* labelled_clone() const = 0;
    virtual Void draw(CanvasInterface& c, const Variables2d& p) const = 0;
};

} // namespace Ariadne

#endif
