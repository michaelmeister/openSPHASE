#ifndef HASHGRID_H
#define HASHGRID_H

#include "vector2.h"
#include <vector>
#include <boost/unordered_map.hpp>
#include <tbb/concurrent_hash_map.h>

struct coord {
    int x, y;

    inline bool operator==(const coord &other) const {
        return (x == other.x) && (y == other.y);
    }
};

struct hasher {

    hasher() : p1(73856093), p2(19349663) {
    }

    size_t operator() (const coord &c) const {
        return c.x*p1 ^ c.y*p2;
    }

    bool operator()(const coord &x1, const coord &x2) const {
        if (x1.x < x2.x)
            return true;
        if (x1.x > x2.x)
            return false;
        if (x1.y < x1.y)
            return true;
        return false;
    }

    size_t p1, p2;
};

struct tbb_hasher {
    inline static size_t hash( const coord& c) { return c.x*p1 ^ c.y*p2; }
    inline static bool equal( const coord& a, const coord& b ) { return a.x == b.x && a.y == b.y; }
    static const size_t p1 = 73856093;
    static const size_t p2 = 19349663;
};

/*size_t tbb_hasher::p1 = 73856093;
size_t tbb_hasher::p2 = 19349663;*/



template<typename T>
class HashGrid {
public:
    typedef tbb::concurrent_hash_map<coord, std::vector<T>, tbb_hasher > HashMap;

    HashGrid(const HashGrid &other) : h(other.h) {
        entries = other.entries;
    }

    HashGrid(float h, int n)
        : h(h), entries(n) {

    }

    void insert(const Vector2<float> &v, T data) {
        coord c={int(v.x[0]/h), int(v.x[1]/h)};
        typename HashMap::accessor a;
        if (!entries.find(a, c)) {
            entries.insert(a, c);
            a->second = std::vector<T>();
        }
        a->second.push_back(data);
    }

    std::vector<T> get(const Vector2<float> &v) {
        coord c={int(v.x[0]/h), int(v.x[1]/h)};
        typename HashMap::accessor a;

        if (entries.find(a, c)) {
            return a->second;
        }
        return std::vector<T>();
    }

    std::vector<T> neighbors(const Vector2<float> &v) {
        std::vector<T> n;
        int cx = int(v.x[0]/h);
        int cy = int(v.x[1]/h);
        for (int x = -1; x < 2; x++) {
            for (int y = -1; y < 2; y++) {
                coord c={cx+x, cy+y};
                typename HashMap::accessor a;
                if (entries.find(a, c)) {
                    std::vector<T> &tmp = a->second;
                    n.insert(n.end(), tmp.begin(), tmp.end());
                }
            }
        }
        return n;
    }

private:
    HashGrid() {

    }

    float h;
    HashMap entries;

};

#endif // HASHGRID_H
