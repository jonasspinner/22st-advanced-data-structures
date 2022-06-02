#include <random>

#include <gtest/gtest.h>

#include "utils.h"

TEST(UtilsTest, Select) {
    ASSERT_EQ(ads::select<uint64_t>(0b1, 1), 0);
    ASSERT_EQ(ads::select64_v2(0b1, 1), 0);
    ASSERT_EQ(ads::select64_v2(0b1111, 1), 0);
    ASSERT_EQ(ads::select64_v2(0b1111, 4), 3);
    ASSERT_EQ(ads::select64_v2(0b111111111111111111, 1), 0);
    ASSERT_EQ(ads::select64_v2(0b111111111111111111, 8), 7);

    std::mt19937_64 gen;
    std::uniform_int_distribution<uint64_t> word_dist;
    for (int k = 0; k < 1000; ++k) {
        uint64_t word = word_dist(gen);
        std::uniform_int_distribution<int> index_dist(1, std::popcount(word));
        auto i = index_dist(gen);
        if (k == 0) std::cout << word << " " << i << "\n";
        ASSERT_EQ(ads::select64_v2(word, i), ads::select(word, i)) << k << " " << word;
    }
}