#include <gtest/gtest.h>
#include <unigrid.h>
#include <time.h>
#include "particle2.h"

float frandom()
{
    srand(time(0));
    return rand() / (float) RAND_MAX;
}

TEST(UniGrid, basic) {
    std::pair<float, float> dim(100, 100);

    std::vector<Particle2<number> *> particles;

    for (int i = 0; i < 10; ++i) {
        float x = frandom()*10 + 45;
        float y = frandom()*10 + 45;
        particles.push_back(new Particle2<float>(Vector2<float>(x, y), 1.0));
    }
    particles.push_back(new Particle2<number>(Vector2<number>(50, 50)));

    UniGrid grid(&particles, dim, 10.0);

    ASSERT_EQ(grid.neighbours(11).size(), 10);
}
