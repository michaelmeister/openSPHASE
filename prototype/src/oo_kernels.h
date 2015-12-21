#ifndef OO_KERNELS_H
#define OO_KERNELS_H

#include <cmath>

template<typename number>
class Spiky {
    Spiky(number h) : h(h) {
        kernel_actor = 10.0/(h*h*h*h*h*M_PI);
        deriv_factor = -3*kernel_factor; //-30/(h^5*Pi)
    }

    number kernel(number r) const {
        number hr = h - r;
        return kernel_factor * hr * hr * hr;
    }

    number derivative(number r) const {
        number hr = h - r;
        return deriv_factor * hr * hr;
    }

    number cutoff() const {
        return h;
    }

    number h, h5;
    number kernel_factor, deriv_factor;
};

#endif // OO_KERNELS_H
