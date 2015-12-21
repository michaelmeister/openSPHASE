
#include "unigrid.h"

#define std_min std::numeric_limits<double>::min()
#define std_max std::numeric_limits<double>::max()

struct minmax_part_position {
    inline bool operator ()(const Particle2 *i, const Particle2 *j) const {
        return i->pos < j->pos;
    }
};
/*
struct max_dist {
    inline bool operator()(const NeighbourData &n_i, const NeighbourData &n_j) const {
        return n_i.dist < n_j.dist;
    }
};*/


UniGrid::UniGrid(std::pair<double, double> scene_dim,
                 const Kernel *kernel,
                 double h)
    : kernel(kernel), h(h), hh(h*h) {

    /*auto minmax =
        thrust::minmax_element(particles->begin(), particles->end(), minmax_part_position());

        vmin = (*minmax.first)->pos;
        vmax = (*minmax.second)->pos;*/

    grid_dim = IntPair((size_t)std::ceil(scene_dim.first / h), (size_t)std::ceil(scene_dim.second / h));
    grid.resize(grid_dim.first * grid_dim.second);
}

void UniGrid::update(Particles *ps) {
    particles = ps;
    if (cell_indices.size() != particles->size())
        cell_indices = std::vector<size_t>(particles->size(), 0);
    calc_cell_indices();
    thrust::stable_sort_by_key(thrust::retag<thrust::omp::tag>(cell_indices.begin()),
                               thrust::retag<thrust::omp::tag>(cell_indices.end()),
                               thrust::retag<thrust::omp::tag>(particles->begin()));
    set_cell_starts();
}

void UniGrid::neighbours(   size_t p_i, std::vector<NeighbourData > &nb,
                            const Vector2 &displ, double distLimit) const
{
    Particle2 *p = particles->at(p_i);

    size_t i = cell_indices[p_i];
    assert(grid[i] >= 0 && grid[i] < cell_indices.size());

    size_t center_x = i % grid_dim.first;
    size_t center_y = i / grid_dim.first;

    for (int y = -1; y < 2; y++) {
        for (int x = -1; x < 2; x++) {
            size_t idx_x = center_x + x;
            size_t idx_y = center_y + y;

            /* Check if position is outside of scene. */
            if (idx_x >= grid_dim.first || idx_y >= grid_dim.second || idx_x < 0 || idx_y < 0) {
                continue;
            }
            size_t id = idx_x + idx_y * grid_dim.first;
            if (cell_indices[grid[id]] == id)
                insert_from_cell(id, p_i, nb, Vector2(0, 0), distLimit);
            if (displ.length_squared() == 0)
                continue;

            /* periodic continuation at periodic boundaries */
            idx_x = (idx(p->pos - displ) % grid_dim.first) + x;
            idx_y = (idx(p->pos - displ) / grid_dim.first) + y;

            /* Check if position is outside of scene. */
            if (idx_x >= grid_dim.first || idx_y >= grid_dim.second || idx_x < 0 || idx_y < 0) {
                continue;
            }

            id = idx_x + idx_y * grid_dim.first;
            if (cell_indices[grid[id]] == id)
                insert_from_cell(id, p_i, nb, displ, distLimit);
        }
    }
}

/*
 * This method stores all indices of the neighbours of the particle at position p_i to nb.
 */
void UniGrid::insert_from_cell(size_t ci, size_t p_i, std::vector<NeighbourData > &nb,
                               const Vector2 displ, double distLimit) const
{
    Particle2 *p = particles->at(p_i);
    assert(ci >= 0 && ci < grid.size());
    size_t start = grid[ci];

    assert(start >= 0);
    assert(start < cell_indices.size());

    const double distLimit2 = distLimit*distLimit;
    const Vector2& p_i_pos = p->pos - displ;
    Vector2 dir;
    double dist2;

    nb.reserve(10);
    size_t maxIndex = start+cell_content[ci];

    for (size_t i = start; i < maxIndex; i++) {
        Particle2 *p_tmp = particles->at(i);
        dir = p_i_pos - p_tmp->pos;
        dist2 = dir.length_squared();
        if (dist2 > distLimit2)
            continue;

        nb.push_back(NeighbourData(i, displ));
    }
}

size_t UniGrid::idx(const Vector2 &v) const {
    assert(!std::isnan(v.x));
    assert(!std::isnan(v.y));
    size_t x, y;
    //assert(v.x >= 0);
    //assert(v.y >= 0);
    x = size_t(v.x / h);
    y = size_t(v.y / h);
    //assert(x < grid_dim.first);
    //assert(y < grid_dim.second);
    size_t ci = x + y * grid_dim.first;
    //assert(ci < grid.size());
    return ci;
}

void UniGrid::calc_cell_indices() {
    // calculate min
    Vector2 minPos(0, 0);
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        minPos.x = (std::min)(minPos.x, p_i->pos.x);
        minPos.y = (std::min)(minPos.y, p_i->pos.y);
    }
    // reset content
    for (long long i = 0; i < (long long)cell_content.size(); i++) {
        cell_content[i] = 0;
    }

    cell_content.resize(grid.size());

    posOffset = minPos;
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        size_t c_i = idx(p_i->pos - posOffset);
        cell_indices[i] = c_i;
        cell_content[c_i]++;
    }
}

void UniGrid::set_cell_starts() {
    if (!particles->size()) return;
    grid[cell_indices[0]] = 0;
    for (long long i = 1; i < (long long)cell_indices.size(); i++) {
        if (cell_indices[i - 1] != cell_indices[i]) {
            unsigned long long ci = cell_indices[i];
            //assert(ci < grid.size());
            grid[ci] = i;
            continue;
        }
    }
}

Sceneparser::Vector2<size_t> UniGrid::getDim() const
{
    return Sceneparser::Vector2<size_t>(grid_dim.first, grid_dim.second);
}

Vector2 UniGrid::getOffset() const
{
    return posOffset;
}
