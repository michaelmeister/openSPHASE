#ifndef KERNELPOLY6_H
#define KERNELPOLY6_H

#include "vector2.h"
#include <assert.h>

template<typename number = float>
class KernelPoly6 {
public:
    KernelPoly6(number h): h(h) {
        hh = h*h;
        hhh = hh*h;
    }

    number operator()(const Vector2<number> &r) const {
        number l = r.length();
        if (l > h)
            return 0;
        number ll = l*l;
        return 4.0/(M_PI*h*h)*pow((1-ll/hh), 3);
    }

    inline Vector2<number> grad(Vector2<number> r) const {
        assert(false);
    }

    inline number laplacian(number r) const {
        assert(false);
    }

private:
    number h, hh, hhh;
};


#endif // KERNELPOLY6_H
