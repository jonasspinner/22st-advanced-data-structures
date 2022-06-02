#include <vector>

#include <gtest/gtest.h>
#include <bitset>
#include <random>

#include "NaiveDynamicBitVector.h"
#include "DynamicBitVector.h"

TEST(NaiveBVTest, test1) {
    using ads::NaiveDynamicBitVector;

    std::vector<bool> bits{false, true, true, true};

    NaiveDynamicBitVector bv(bits);

    ASSERT_EQ(bv.size(), bits.size());

    for (size_t i = 0; i < bits.size(); ++i) {
        ASSERT_EQ(bv.access(i), bits[i]);
    }

    ASSERT_EQ(bv.rank(0, false), 0);
    ASSERT_EQ(bv.rank(0, true), 0);
    ASSERT_EQ(bv.rank(1, false), 1);
    ASSERT_EQ(bv.rank(1, true), 0);
    ASSERT_EQ(bv.rank(2, false), 1);
    ASSERT_EQ(bv.rank(2, true), 1);
    ASSERT_EQ(bv.rank(3, false), 1);
    ASSERT_EQ(bv.rank(3, true), 2);
    ASSERT_EQ(bv.rank(4, false), 1);
    ASSERT_EQ(bv.rank(4, true), 3);

    ASSERT_EQ(bv.select(1, false), 0);
    ASSERT_EQ(bv.select(1, true), 1);
    ASSERT_EQ(bv.select(2, true), 2);
    ASSERT_EQ(bv.select(3, true), 3);

    bv.insert(4, false);
    ASSERT_EQ(bv.size(), 5);
    ASSERT_EQ(bv.access(4), false);

    bv.remove(3);
    ASSERT_EQ(bv.size(), 4);
    ASSERT_EQ(bv.access(1), true);
    ASSERT_EQ(bv.access(2), true);
    ASSERT_EQ(bv.access(3), false);
}

TEST(NaiveBVTest, Lecture01) {
    using ads::NaiveDynamicBitVector;

    std::vector<bool> bits{false, true, true, false, true, true, false, true, false, false};

    NaiveDynamicBitVector bv(bits);

    ASSERT_EQ(bv.size(), bits.size());

    for (size_t i = 0; i < bits.size(); ++i) {
        ASSERT_EQ(bv.access(i), bits[i]);
    }

    ASSERT_EQ(bv.rank(5, false), 2);
    ASSERT_EQ(bv.rank(5, true), 3);

    ASSERT_EQ(bv.select(5, true), 7);
}


TEST(BVTest, Leaf) {
    using BV = ads::DynamicBitVector<uint8_t, 2>;
    using Inner = BV::Inner;
    using Leaf = BV::Leaf;

    auto eq = [](const Leaf &leaf, uint64_t x, size_t size) {
        if (leaf.size() != size) return false;
        std::bitset<64> bs(x);
        for (size_t i = 0; i < leaf.size(); ++i) {
            if (bs[leaf.size() - i - 1] != leaf.access(i)) {
                std::cerr << "bs[" << leaf.size() - i - 1 << "] = " << bs[leaf.size() - i - 1] << " leaf[" << i << "] = " << leaf.access(i) << std::endl;
                return false;
            }
        }
        return true;
    };

    Leaf leaf;
    leaf.m_bits = {0, 0};
    leaf.m_size = 12;

    ASSERT_TRUE(eq(leaf, 0b0000'0000'0000, 12));

    leaf.flip(5);
    leaf.flip(7);
    leaf.flip(8);
    leaf.flip(11);

    ASSERT_TRUE(eq(leaf, 0b0000'0101'1001, 12));

    ASSERT_EQ(leaf.rank1(0), 0);
    ASSERT_EQ(leaf.rank1(5), 0);
    ASSERT_EQ(leaf.rank1(6), 1);
    ASSERT_EQ(leaf.rank1(7), 1);
    ASSERT_EQ(leaf.rank1(8), 2);
    ASSERT_EQ(leaf.rank1(9), 3);
    ASSERT_EQ(leaf.rank1(11), 3);
    ASSERT_EQ(leaf.rank1(12), 4);

    ASSERT_EQ(leaf.select1(1), 5);
    ASSERT_EQ(leaf.select1(2), 7);
    ASSERT_EQ(leaf.select1(3), 8);
    ASSERT_EQ(leaf.select1(4), 11);

    ASSERT_EQ(leaf.select0(1), 0);
    ASSERT_EQ(leaf.select0(2), 1);
    ASSERT_EQ(leaf.select0(3), 2);
    ASSERT_EQ(leaf.select0(4), 3);
    ASSERT_EQ(leaf.select0(5), 4);
    ASSERT_EQ(leaf.select0(6), 6);
    ASSERT_EQ(leaf.select0(7), 9);
    ASSERT_EQ(leaf.select0(8), 10);

    leaf.remove(3);
    ASSERT_TRUE(eq(leaf, 0b0000'1011'001, 11));
    leaf.remove(10);
    ASSERT_TRUE(eq(leaf, 0b0000'1011'00, 10));
    leaf.remove(7);
    ASSERT_TRUE(eq(leaf, 0b0000'1010'0, 9));

    bool overflow{}, b{};

    std::tie(overflow, b) = leaf.insert(9, true);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b0000'1010'01, 10));

    std::tie(overflow, b) = leaf.insert(0, true);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b1000'0101'001, 11));

    std::tie(overflow, b) = leaf.insert(8, true);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b1000'0101'1001, 12));

    std::tie(overflow, b) = leaf.insert(7, true);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b1000'0101'1100'1, 13));

    std::tie(overflow, b) = leaf.insert(13, true);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b1000'0101'1100'11, 14));

    std::tie(overflow, b) = leaf.insert(0, false);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b0100'0010'1110'011, 15));

    std::tie(overflow, b) = leaf.insert(1, true);
    ASSERT_FALSE(overflow);
    ASSERT_TRUE(eq(leaf, 0b0110'0001'0111'0011, 16));

    std::tie(overflow, b) = leaf.insert(2, false);
    ASSERT_TRUE(overflow);
    ASSERT_TRUE(b);
    ASSERT_TRUE(eq(leaf, 0b0101'0000'1011'1001, 16));

    std::tie(overflow, b) = leaf.insert(15, false);
    ASSERT_TRUE(overflow);
    ASSERT_TRUE(b);
    ASSERT_TRUE(eq(leaf, 0b0101'0000'1011'1000, 16));
}

TEST(BVTest, BitMasks) {
    using BV = ads::DynamicBitVector<uint64_t, 4>;
    using Leaf = BV::Leaf;

    ASSERT_EQ(Leaf::bit_mask(0), 0b1);
    ASSERT_EQ(Leaf::bit_mask(63), 0b1000000000000000000000000000000000000000000000000000000000000000);
    ASSERT_EQ(Leaf::lower_bit_mask(1), 0b1);
    ASSERT_EQ(Leaf::lower_bit_mask(63), 0b0111111111111111111111111111111111111111111111111111111111111111);
    ASSERT_EQ(Leaf::lower_bit_mask(64), 0b1111111111111111111111111111111111111111111111111111111111111111);
}

TEST(BVTest, Insert) {
    using BV = ads::DynamicBitVector<>;

    auto eq = [](const auto &v, uint64_t x, size_t size) {
        if (v.size() != size) return false;
        std::bitset<64> bs(x);
        for (size_t i = 0; i < v.size(); ++i) {
            if (bs[v.size() - i - 1] != v.access(i)) {
                std::cerr << "bs[" << v.size() - i - 1 << "] = " << bs[v.size() - i - 1] << " bv[" << i << "] = " << v.access(i) << std::endl;
                return false;
            }
        }
        return true;
    };

    BV bv;

    ASSERT_EQ(bv.size(), 0);
    ASSERT_TRUE(eq(bv, 0b0, 0));

    for (size_t i = 0; i < 8; ++i) {
        bv.insert(i, (i % 2) == 0);
        ASSERT_EQ(bv.size(), i + 1);
    }
    ASSERT_TRUE(eq(bv, 0b10101010, 8));

    for (size_t i = 8; i < 16; ++i) {
        bv.insert(i, true);
        ASSERT_EQ(bv.size(), i + 1);
    }
    ASSERT_TRUE(eq(bv, 0b1010101011111111, 16));

    std::cout << "+++\n" << bv << "---\n";

    bv.insert(16, false);

    ASSERT_TRUE(eq(bv, 0b10101010111111110, 17));

    std::cout << "+++\n" << bv << "---\n";

    ASSERT_EQ(bv.size(), 17);
    ASSERT_FALSE(bv.access(16));

    for (size_t i = 0; i < 100; ++i) {
        bv.insert(bv.size(), (i % 3) == 0);
    }

    std::cout << "+++\n" << bv << "---\n";

    std::cout << bv.rank(100, true) << "\n";
    std::cout << bv.rank(100, false) << "\n";
}

TEST(BVTest, Inner) {
    using BV = ads::DynamicBitVector<>;

    BV bv;
    for (int i = 0; i < 36; ++i) {
        bv.insert(0, true);
    }

    std::cout << "+++\n" << bv << "---\n";

    bv.insert(0, true);
}

TEST(BVTest, Remove) {
    using BV = ads::DynamicBitVector<>;

    BV bv;

    for (size_t i = 0; i < 64; ++i) {
        bv.insert(i, (i % 2) == 0);
    }

    ASSERT_EQ(bv.size(), 64);

    std::cout << "+++\n" << bv << "---\n";

    for (size_t i = 63; i > 0; --i) {
        if ((i % 2) == 0) {
            bv.remove(i);
        }
    }
    bv.remove(0);

    std::cout << "+++\n" << bv << "---\n";

    for (size_t i = 0; i < 32; ++i) {
        bv.remove(0);
    }

    std::cout << "+++\n" << bv << "---\n";

    ASSERT_EQ(bv.size(), 0);
    std::mt19937_64 gen;
    std::bernoulli_distribution bit_dist;
    size_t n = 1000;
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        size_t index = index_dist(gen);
        bool bit = bit_dist(gen);
        bv.insert(index, bit);
        ASSERT_EQ(bv.access(index), bit);
    }
    ASSERT_EQ(bv.size(), n);

    for (int k = 0; k < 10; ++k) {
        for (size_t i = 0; i < n; ++i) {
            std::uniform_int_distribution<size_t> index_dist(0, bv.size());
            bv.insert(index_dist(gen), bit_dist(gen));
        }
        for (size_t i = 0; i < n; ++i) {
            std::uniform_int_distribution<size_t> index_dist(0, bv.size() - 1);
            bv.remove(index_dist(gen));
        }
    }
    ASSERT_EQ(bv.size(), n);

    std::cout << bv.size() << " " << bv.rank(bv.size(), true) << std::endl;
}



TEST(BVTest, RankSelect) {
    using BV = ads::DynamicBitVector<>;

    BV bv;

    std::mt19937_64 gen;
    std::bernoulli_distribution bit_dist;
    size_t n = 1000;
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        bv.insert(index_dist(gen), bit_dist(gen));
    }

    ASSERT_EQ(bv.rank(bv.size(), true) + bv.rank(bv.size(), false), bv.size());

    auto num_ones = bv.rank(bv.size(), true);
    auto num_zeros = bv.rank(bv.size(), false);

    ASSERT_EQ(bv.rank(bv.select(num_zeros, false) + 1, false), num_zeros);

    std::uniform_int_distribution<size_t> index_dist(0, bv.size() - 1);
    std::uniform_int_distribution<size_t> rank0_dist(1, num_zeros);
    std::uniform_int_distribution<size_t> rank1_dist(1, num_ones);
    for (int k = 0; k < 100; ++k) {
        auto rank0 = rank0_dist(gen);
        ASSERT_EQ(bv.rank(bv.select(rank0, false) + 1, false), rank0);
        auto rank1 = rank1_dist(gen);
        ASSERT_EQ(bv.rank(bv.select(rank1, true) + 1, true), rank1);

        auto idx = index_dist(gen);
        auto prev0 = bv.select(bv.rank(idx, false), false);
        auto prev1 = bv.select(bv.rank(idx, true), true);
        ASSERT_EQ(bv.select(bv.rank(prev0, false) + 1, false), prev0);
        ASSERT_EQ(bv.select(bv.rank(prev1, true) + 1, true), prev1);
    }

    std::cout << bv.size() << " " << bv.required_bits_upperbound() << std::endl;
}

