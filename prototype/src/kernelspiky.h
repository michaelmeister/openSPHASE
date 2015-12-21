#ifndef KERNELSPIKY_H
#define KERNELSPIKY_H

#include "vector2.h"
#include <cmath>

template<typename number>
class KernelSpiky {
public:
    KernelSpiky(number h) : h(h) {
    }

    number operator()(const Vector2<number> &r) {
        number d = r.length();
        if (d > h) {
            return 0;
        }
        float x = 1.f - d/h;
        number result = 10.f / (M_PI*h*h) * x * x * x;
        assert(!isnan(result));
        return result;
    }

    Vector2<number> grad(Vector2<number> &r) {
        number _r = r.length();
        number h5 = h*h*h*h*h;
        number hr = h-_r;
        return -r*10/(M_PI*h5)*hr*hr*hr;
        /*number d = r.length();
        if (d > h) {
            return Vector2<number>();
        }
        number x = 1.f - d/h;
        number first_deriv = -30.f / (M_PI*h*h*h) * x * x;
        Vector2<number> result = -first_deriv * r.normalized();
        assert(!isnan(result.x[0]));
        assert(!isnan(result.x[1]));

        return result;*/
    }

    number h;
};

#endif // KERNELSPIKY_H
