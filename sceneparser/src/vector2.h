#ifndef VECTOR2_H
#define VECTOR2_H

#include <cmath>
#include <ostream>

namespace Sceneparser {

template <typename number>
class Vector2 {
public:
    inline Vector2() {
        clear();
    }

    inline Vector2(const Vector2<number> &other) {
        x = other.x;
        y = other.y;
    }

    inline Vector2(number x, number y) {
        this->x = x;
        this->y = y;
    }

    inline Vector2<number> &operator=(const Vector2<number> &other) {
        x = other.x;
        y = other.y;
        return *this;
    }

    inline Vector2<number> operator-(const Vector2<number> &other) const {
        return Vector2<number>(x-other.x, y-other.y);
    }

    inline Vector2<number> operator+(const Vector2<number> &other) const {
        return Vector2<number>(x+other.x, y+other.y);
    }

    inline void operator+=(const Vector2<number> &other) {
        x+=other.x;
        y+=other.y;
    }

    inline void operator-=(const Vector2<number> &other) {
        x-=other.x;
        y-=other.y;
    }

    inline Vector2<number> operator*(number scalar) const {
        return Vector2<number>(x*scalar, y*scalar);
    }

    inline number operator*(const Vector2<number> &other) const {
        return x*other.x + y*other.y;
    }

    inline Vector2<number> operator/(number scalar) const {
        return Vector2<number>(x/scalar, y/scalar);
    }

    inline void operator/=(number scalar) {
        x /= scalar;
        y /= scalar;
    }

    inline Vector2<number> normalized() const {
        number l = length();
        return Vector2<number>(x/l, y/l);
    }

    inline Vector2<number> normalized(number c) const {
        return Vector2<number>(x/c, y/c);
    }

    inline Vector2<number> operator-() const {
        return Vector2<number>(-x, -y);
    }

    inline number length() const {
        return std::sqrt(x*x + y*y);
    }

    inline number length_squared() const {
        return x*x + y*y;
    }

    inline void clear() {
        x = y = 0.0;
    }

    static inline number dot(const Vector2<number> &v, const Vector2<number> &w) {
        return v.x*w.x + v.y*w.y;
    }


    inline Vector2<number> rot(number theta) const {
        return Vector2<number>(cos(theta)*x - sin(theta)*y, sin(theta)*x + cos(theta)*y);
    }

    static inline number angle(const Vector2<number> &v, const Vector2<number> &w) {
        number length = v.length() * w.length();
        return std::acos(Vector2<number>::dot(v, w)/length);
    }

    inline bool operator==(const Vector2<number> &other) const {
        return x == other.x && y == other.y;
    }

    inline bool operator<(const Vector2<number> &other) const {
        if (x < other.x)
            return true;
        if (x == other.x)
            return y < other.y;
        return false;
    }

    number x, y;
};

template <typename number>
Vector2<number>  operator*(number scalar, const Vector2<number> &right) {
    return Vector2<number> (scalar * right.x, scalar*right.y);
}

template <typename number>
Vector2<number> operator/(number scalar, const Vector2<number> &right) {
    return Vector2<number> (scalar / right.x, scalar / right.y);
}

template <typename number>
std::ostream &operator<<(std::ostream &c, const Vector2<number>  &v) {
    c << "Vector[" << v.x << ", " << v.y << "]";
    return c;
}

template<typename number>
class Matrix2 {
public:

    Matrix2() : m11(0.0), m12(0.0), m21(0.0), m22(0.0) {
    }

    Matrix2(number m11, number m12,
            number m21, number m22) : m11(m11), m12(m12), m21(m21), m22(m22) {
    }

    inline
    Matrix2 inverse() const {
        number scale = number(1.0)/(m11*m22 - m12*m21);
        return Matrix2(scale*m22, -scale*m12, -scale*m21, scale*m11);
    }

    inline
    Matrix2 operator*(number alpha) const {
        return Matrix2(alpha*m11, alpha*m12, alpha*m21, alpha*m22);
    }

    inline
    Vector2<number> operator*(Vector2<number> x) const {
        return Vector2<number>(m11*x.x+m12*x.y, m21*x.x+m22*x.y);
    }

    inline
    Matrix2 operator*(const Matrix2 &m) const {
        return Matrix2(m11*m.m11+m12*m.m21,
                       m11*m.m12+m12*m.m22,
                       m21*m.m11+m22*m.m21,
                       m21*m.m12+m22*m.m22);
    }

    inline
    Matrix2 operator+(const Matrix2 &m) const {
        return Matrix2(m.m11+m11, m.m12+m12, m.m21+m21, m.m22+m22);
    }

    inline
    Matrix2&
    operator+=(const Matrix2 &m) {
        m11 += m.m11;
        m12 += m.m12;
        m21 += m.m21;
        m22 += m.m22;
        return *this;
    }

    inline
    bool operator==(const Matrix2 &m) const {
        return m11 == m.m11 && m12 == m.m12 && m21 == m.m21 && m22 == m.m22;
    }

    /**
      * return identity Matrix
      */
    static Matrix2 identity() {
        return Matrix2(1.0, 0, 0,1);
    }

//#ifndef NDEBUG
//    inline
//    void assert_not_nan() const {
//        assert(!std::isnan(m11));
//        assert(!std::isnan(m12));
//        assert(!std::isnan(m21));
//        assert(!std::isnan(m22));
//    }

//#endif

    number m11, m12, m21, m22;
};

template<typename number>
std::ostream &operator<<(std::ostream &out, const Matrix2<number> &m) {
    out << "Matrix[" << m.m11 << " " << m.m12 <<
           ";" << m.m21 << " " << m.m22 << "]";
}

}

#endif // VECTOR2_H
