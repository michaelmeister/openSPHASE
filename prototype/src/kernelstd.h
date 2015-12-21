#ifndef KERNELSTD_H
#define KERNELSTD_H

#include "vector2.h"

template <typename number>
class KernelStd {
public:
    KernelStd(number h) : h(h) {

    }

    number operator()(const Vector2<number> &r) {
        number d = r.length();
        if (d >= h)
          return 0.f;
        else {
          float x = 1.f - d/h;
          return 4.f/(M_PI*h*h) * x * x * x;
        }
    }

    number laplacian(const Vector2<number> &r) {
        number d = r.length();
        if (d >= h)
          return 0.f;
        else {
          number x = d/h;
          return 24.f / (M_PI*h*h*h*h) * (1 - x) * (5 * x - 1);
        }
    }

    number h;
};

#endif // KERNELSTD_H
