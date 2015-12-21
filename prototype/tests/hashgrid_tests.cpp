#include "hashgrid.h"
#include <gtest/gtest.h>
#include <iostream>

TEST(HashGridTests, can_insert) {
    HashGrid<int> hg(10, 4);
    hg.insert(Vector2<float>(10, 50), 1050);
    ASSERT_TRUE(hg.get(Vector2<float>(10, 50)).size() == 1);
}

TEST(HashGridTests, does_grid) {
    HashGrid<float> hg(10, 4);
    for (float f = 0.0; f < 10.0; f += 1.0) {
        hg.insert(Vector2<float>(f, f), f);
    }

    std::vector<float> results = hg.get(Vector2<float>(0.5, 0.5));
    ASSERT_EQ(results.size(), 10);
    ASSERT_EQ(results[0], 0.0);

    for (int i = 0; i < 10; i++) {
        ASSERT_EQ(results[i], i);
    }
}

TEST(HashGridTests, neighbours) {
    HashGrid<float> hg(1.0, 4);

    for (int x  = -1; x < 2; x++) {
        for (int y  = -1; y < 2; y++) {
            hg.insert(Vector2<float>(x, y), x*y);
        }
    }

    ASSERT_EQ(hg.neighbors(Vector2<float>(0.0, 0.0)).size(), 9);
}
