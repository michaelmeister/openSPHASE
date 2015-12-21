#ifndef NEIGHBOURDATA_H
#define NEIGHBOURDATA_H

#include "stable.h"

#include <assert.h>
#include <signal.h>
#include "kernels.h"

// TODO templatize!
struct NeighbourData {

    NeighbourData(const Vector2 &dir,
                  double dist, size_t j, const Kernel *kernel, const Vector2 &displ) {

        this->j = j;
        this->dir = dir;
        assert(!std::isnan(this->dir.x));
        assert(!std::isnan(this->dir.y));
        this->dist = dist;
        assert(!std::isnan(this->dist));
        this->displ = displ;

        //quintic or wendland
        this->kernel = kernel->f(dist);
        this->kernel_deriv = kernel->df(dist);
        this->kernel_gradient = dir / (dist)* kernel_deriv;
    }

    void update_kernel_f_only(const Kernel *kernel, const Vector2& dir)
    {
        this->dir = dir;
        this->dist = dir.length();
        this->kernel = kernel->f(dist);
#ifndef NDEBUG
        this->kernel_deriv = 0;
        this->kernel_gradient.x = 0;
        this->kernel_gradient.y = 0;
#endif
    }

    void update_kernel_df_only(const Kernel *kernel, const Vector2& dir)
    {
        this->dir = dir;
        this->dist = dir.length();
        this->kernel_deriv = kernel->df(dist);
        this->kernel_gradient = dir / (dist)* kernel_deriv;
    }

    void update_kernel(const Kernel *kernel, const Vector2& dir) {
        this->dir = dir;
        this->dist = dir.length();
        this->kernel = kernel->f(dist);
        this->kernel_deriv = kernel->df(dist);
        this->kernel_gradient = dir / (dist)* kernel_deriv;
    }

    NeighbourData() {
    }

    NeighbourData(size_t j, const Vector2 &displ) :
        j(j), displ(displ)
  #ifndef NDEBUG
      ,dist(0), kernel(0), kernel_deriv(0), kernel_gradient(0,0)
  #endif
    {}

    Vector2 dir, kernel_gradient, displ;
    double dist, kernel, kernel_deriv;
    size_t j;
};

typedef std::vector<std::vector<NeighbourData> > Neighbours;

#endif // NEIGHBOURDATA_H
