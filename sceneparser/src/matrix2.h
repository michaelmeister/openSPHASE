#ifndef MATRIX2_H
#define MATRIX2_H

#include <ostream>
#include "vector2.h"

namespace Sceneparser {

template <typename number>
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
    Vector2 operator*(Vector2 x) const {
        return Vector2(m11*x.x+m12*x.y, m21*x.x+m22*x.y);
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

std::ostream &operator<<(std::ostream &out, const Matrix2 &m) {
    out << "Matrix[" << m.m11 << " " << m.m12 <<
           ";" << m.m21 << " " << m.m22 << "]";
}

}

#endif // MATRIX2_H
