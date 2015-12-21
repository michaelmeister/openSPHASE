#ifndef KERNELVISCOSITY_H
#define KERNELVISCOSITY_H

#include "vector2.h"

template<class number>
class KernelViscosity
{
public:
    KernelViscosity(number h) : h(h){
        factor_laplacian = 45.0/(PI*pow(h, 6));
    }

    number operator()(const Vector2<number> &r) {
        return 1.0;
    }

    Vector2<number> grad(const Vector2<number> &r) {
        return Vector2<number>();
    }

    number laplacian(const Vector2<number> &r) {
        number l = r.length();
        if (l < h)
            return 0.0;
        return factor_laplacian * (h - l);
    }
    number h, factor_laplacian;
};

#endif // KERNELVISCOSITY_H
