#ifndef KERNELS_H
#define KERNELS_H

#include "assert.h"
#include <string>

struct Kernel {
    virtual double f(double r) const = 0;
    virtual double df(double r) const = 0;

    static Kernel *get(const std::string &name, double h);
};
#endif // KERNELS_H
