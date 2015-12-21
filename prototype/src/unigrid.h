#ifndef UNIGRID_H
#define UNIGRID_H

#include "stable.h"

#include <vector>
#include <boost/bind.hpp>
#include <assert.h>
#include <limits>
#include <algorithm>
#include <thrust/sort.h>
#include <thrust/system/omp/memory.h>
#include <thrust/iterator/retag.h>
#include <thrust/fill.h>

#include "neighbourdata.h"
#include "kernels.h"

class UniGrid {
    typedef std::pair<size_t, size_t> IntPair;
public:
    UniGrid(std::pair<double, double> scene_dim, const Kernel *kernel, double h);
    void update(Particles *ps);
    void neighbours(size_t p_i, std::vector<NeighbourData> &nb,
                    const Vector2 &displ, double distLimit) const;

    Sceneparser::Vector2<size_t> getDim() const;
    Vector2 getOffset() const;
private:
    void insert_from_cell(size_t ci, size_t p_i,
                          std::vector<NeighbourData> &nb, const Vector2 displ, double distLimit) const;
    size_t idx(const Vector2 &v) const;
    void calc_cell_indices();
    void set_cell_starts();

private:
    Particles*          particles;
    std::vector<size_t> cell_indices;
    std::vector<size_t>	cell_content;
    std::vector<size_t> grid;
    const Kernel*		kernel;
    double				h;
    double				hh;
    IntPair				grid_dim;
    Vector2             posOffset;
};




#endif // UNIGRID_H
