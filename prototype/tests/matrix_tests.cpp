#include <gtest/gtest.h>
#include "matrix2.h"

TEST(Matrix2, mult) {
    Matrix2<double> m = Matrix2<double>::identity();
    Vector2<double> v(1.0, 3.0);
    ASSERT_EQ(v, m*v);
}

TEST(Matrix2, inverse) {
    Matrix2<double> m(1,2,3,4);
    Matrix2<double> m_inv(m.inverse());
    Matrix2<double> id = Matrix2<double>::identity();
    ASSERT_EQ(m*m_inv, id);
}
