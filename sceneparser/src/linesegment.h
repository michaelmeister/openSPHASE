#ifndef LineSegment_H
#define LineSegment_H

#include "vector2.h"

namespace Sceneparser {

template <typename number>
class LineSegment
{
public:
    LineSegment() {
    }

    LineSegment(const Vector2<number> &from, const Vector2<number> &to, const Vector2<number> &velocity, const bool &moving)
        : b(from), m (to - from), wall_velocity(velocity), is_moving(moving) {
        length = Vector2<number>::dot(m, m);
    }

    inline bool cross(const Vector2<number> &p) const {
        number t;
        Vector2<number> diff = p - b;
        t = Vector2<number>::dot(m, diff);
        t = t / length;
        diff = diff - m * t;
        Vector2<number> to = b + m;
        number dx = to.x - b.x;
        number dy = to.y - b.y;
        if (diff*Vector2<number>(-dy, dx).normalized() < 0)
            return false;
        else
            return true;
    }

    inline number dist(const Vector2<number> &p) const {
        number t;
        Vector2<number> diff = p - b;
        t = Vector2<number>::dot(m, diff);
        t = t / length;
        diff = diff - m * t;
        return diff.length();
    }

    inline Vector2<number> normal_distance(const Vector2<number> &p) const {
        number t;
        Vector2<number> diff = p - b;
        t = Vector2<number>::dot(m, diff);
        t = t / length;
        diff = diff - m * t;
        return diff;
    }

    inline Vector2<number> ref(Vector2<number> &p, number t0) const {
        Vector2<number> Q = b + m*t0;
        return Q*2 - p;
    }

    inline Vector2<number> normal() const {
        Vector2<number> to = b + m;
        number dx = to.x - b.x;
        number dy = to.y - b.y;
        return Vector2<number>(-dy, dx).normalized();
    }

    Vector2<number> b, m;
    number length;
    Vector2<number> wall_velocity;
    bool is_moving;

};

}

#endif // LineSegment_H
