#ifndef KERNELCUBIC_H
#define KERNELCUBIC_H

#include "vector2.h"

template<class number>
class KernelCubic {
public:
    KernelCubic(number h) : h(h) {

    }

    number operator()(Vector2<number> r) {
        number z = fabs(r.length())/h;
        if (z >= 2.0) {
            return 0.0;
        }
        if (z < 1.0) {
            number zsq = z * z;
            return 10.0/(7.0*PI*h*h)*(1 - 3.0/2.0 * zsq + 3.0/4.0 * zsq*z);
        }
        return 10.0/(4.0*7.0*PI*h*h)*pow(1.0 - z, 3);
    }

    Vector2<number> grad(Vector2<number> r) {
        number z = fabs(r.length())/h;
        if (z >= 2.0) {
            return Vector2<number>(0.0, 0.0);
        }
        r = r / r.length();
        if (z < 1.0) {
            number zsq = z * z;
            return r*(30.0/(7.0*PI*h*h*h)*(-z + 3.0/4.0 * zsq));
        }
        return r*(-30.0/(4.0*7.0*PI*h*h*h)*pow(2.0 - z, 2.0));
    }

private:
    number h;
};

#endif // KERNELCUBIC_H
