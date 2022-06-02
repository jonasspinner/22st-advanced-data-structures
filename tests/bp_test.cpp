#include <gtest/gtest.h>

#include "BP.h"
#include "BP2.h"

namespace ads {
    TEST(BPTest, Lecture03) {
        using ads::BP;

        BP bp;
        bp.insertchild(0, 2, 1);
        bp.insertchild(0, 3, 0);
        bp.insertchild(0, 2, 3);
    }

    TEST(BPTest, ExerciseSheet) {
        using ads::BP;

        BP bp;
        auto &bv = bp.bv;
        ASSERT_EQ(bv.access(0), false);
        ASSERT_EQ(bv.access(1), true);
        bv.remove(1);
        bv.remove(0);
        for (bool b: {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1}) {
            bv.insert(bp.bv.size(), b);
        }
        ASSERT_EQ(bv.rank(bv.size(), false), bv.rank(bv.size(), true));

        ASSERT_EQ(bp.open_parenthesis_for_node(0), 0);
        ASSERT_EQ(bp.open_parenthesis_for_node(1), 1);
        ASSERT_EQ(bp.open_parenthesis_for_node(2), 3);
        ASSERT_EQ(bp.open_parenthesis_for_node(3), 4);
        ASSERT_EQ(bp.open_parenthesis_for_node(4), 6);
        ASSERT_EQ(bp.open_parenthesis_for_node(5), 7);
        ASSERT_EQ(bp.open_parenthesis_for_node(6), 9);
        ASSERT_EQ(bp.open_parenthesis_for_node(7), 13);
        ASSERT_EQ(bp.open_parenthesis_for_node(8), 15);
        ASSERT_EQ(bp.open_parenthesis_for_node(9), 16);
        ASSERT_EQ(bp.open_parenthesis_for_node(10), 18);

        /*
        ASSERT_EQ(bp.close_parenthesis_for_node(0), 21);
        ASSERT_EQ(bp.close_parenthesis_for_node(1), 2);
        ASSERT_EQ(bp.close_parenthesis_for_node(2), 12);
        ASSERT_EQ(bp.close_parenthesis_for_node(3), 5);
        ASSERT_EQ(bp.close_parenthesis_for_node(4), 11);
        ASSERT_EQ(bp.close_parenthesis_for_node(5), 8);
        ASSERT_EQ(bp.close_parenthesis_for_node(6), 10);
        ASSERT_EQ(bp.close_parenthesis_for_node(7), 14);
        ASSERT_EQ(bp.close_parenthesis_for_node(8), 20);
        ASSERT_EQ(bp.close_parenthesis_for_node(9), 17);
        ASSERT_EQ(bp.close_parenthesis_for_node(10), 19);
         */
    }

    TEST(BP2Test, Leaf) {
        using BP2 = ads::BP2<uint8_t, 2>;
        using Leaf = BP2::Leaf;

        Leaf leaf;
        for (bool b : {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1}) {
            auto [overflow, _] = leaf.insert(leaf.size(), b);
            ASSERT_FALSE(overflow);
        }

        std::cout << leaf.calculate_total_excess() << " " << leaf.calculate_min_excess() << std::endl;

        ASSERT_EQ(leaf.total_excess, 0);
        ASSERT_GE(leaf.min_excess, 1);
    }
}