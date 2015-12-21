#ifndef KERNELPRESSURE_H
#define KERNELPRESSURE_H

#include "vector2.h"
#include <cmath>
#include <assert.h>

template<typename number>
class KernelPressure
{
public:
    KernelPressure(number h) : h(h) {
    }

    inline number operator()(Vector2<number> &r) {
        number l = r.length();
        if (l >= h) {
            return 0.0;
        }
        number x = 1.0 - l/h;
        number result = 10.0 / (M_PI * h*h) * x * x * x;
        assert(!isnan(result));
        return result;
    }

    inline Vector2<number> grad(Vector2<number> r) const {
        number l = r.length();
        if (l >= h) {
            return Vector2<number>();
        }
        number x = 1 - l/h;
        number first_deriv = -30.0f / (M_PI * h*h) * x * x * x;
        Vector2<number> result = -first_deriv * r.normalized();
        assert(!isnan(result.x));
        assert(!isnan(result.y));
        return result;
    }
private:
    number h;
};

#endif // KERNELPRESSURE_H
