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
        using size_type = BP2::size_type;

        {
            Leaf leaf;
            // ab cd ef g   h
            // (()(()(()()))())
            for (bool b: {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1}) {
                auto [overflow, _] = leaf.insert(leaf.size(), b);
                ASSERT_FALSE(overflow);
            }

            //  v   i  b[i]  rank0(i)  rank1(i)  excess(i)  total   min
            //  a   0     0         0         0          0      -     -
            //  b   1     0         1         0          1      |     |
            //      2     1         2         0          2      |     |
            //  c   3     0         2         1          1      4     0
            //  d   4     0         3         1          2      |     |
            //      5     1         4         1          3      |     |
            //  e   6     0         4         2          2      |     |
            //  f   7     0         5         2          3      -     -
            //      8     1         6         2          4      -     -
            //  g   9     0         6         3          3      |     |
            //     10     1         7         3          4      |     |
            //     11     1         7         4          3     -4     -4
            //     12     1         7         5          2      |     |
            //  h  13     0         7         6          1      |     |
            //     14     1         8         6          2      |     |
            //     15     1         8         7          1      -     -
            //     16               8         8          0

            ASSERT_EQ(leaf.total_excess, 0);
            ASSERT_EQ(leaf.min_excess, 0);


            // first '(' with that level
            ASSERT_EQ(leaf.forward_search(-1, 0), std::pair(-1, 0));
            ASSERT_EQ(leaf.forward_search(-1, 1), std::pair(0, 0));
            ASSERT_EQ(leaf.forward_search(-1, 2), std::pair(1, 0));
            ASSERT_EQ(leaf.forward_search(-1, 3), std::pair(4, 0));
            ASSERT_EQ(leaf.forward_search(-1, 4), std::pair(7, 0));
            ASSERT_EQ(leaf.forward_search(-1, -1), std::pair(16, -1));

            // last '(' with that -level
            ASSERT_EQ(leaf.backward_search(16, 0), std::pair(16, 0));
            ASSERT_EQ(leaf.backward_search(16, -1), std::pair(15, 0));
            ASSERT_EQ(leaf.backward_search(16, -2), std::pair(14, 0));
            ASSERT_EQ(leaf.backward_search(16, -3), std::pair(11, 0));
            ASSERT_EQ(leaf.backward_search(16, -4), std::pair(10, 0));
            ASSERT_EQ(leaf.backward_search(16, 1), std::pair(-1, 1));


            // closing parenthesis
            ASSERT_EQ(leaf.forward_search(0, 0), std::pair(15, 0));
            ASSERT_EQ(leaf.forward_search(1, 0), std::pair(2, 0));
            ASSERT_EQ(leaf.forward_search(3, 0), std::pair(12, 0));
            ASSERT_EQ(leaf.forward_search(4, 0), std::pair(5, 0));
            ASSERT_EQ(leaf.forward_search(6, 0), std::pair(11, 0));
            ASSERT_EQ(leaf.forward_search(7, 0), std::pair(8, 0));
            ASSERT_EQ(leaf.forward_search(9, 0), std::pair(10, 0));
            ASSERT_EQ(leaf.forward_search(13, 0), std::pair(14, 0));

            // opening parenthesis
            ASSERT_EQ(leaf.backward_search(2, 0), std::pair(1, 0));
            ASSERT_EQ(leaf.backward_search(5, 0), std::pair(4, 0));
            ASSERT_EQ(leaf.backward_search(8, 0), std::pair(7, 0));
            ASSERT_EQ(leaf.backward_search(10, 0), std::pair(9, 0));
            ASSERT_EQ(leaf.backward_search(11, 0), std::pair(6, 0));
            ASSERT_EQ(leaf.backward_search(12, 0), std::pair(3, 0));
            ASSERT_EQ(leaf.backward_search(14, 0), std::pair(13, 0));
            ASSERT_EQ(leaf.backward_search(15, 0), std::pair(0, 0));

            // enclosing parenthesis
            ASSERT_EQ(leaf.backward_search(1, 2), std::pair(0, 0));
            ASSERT_EQ(leaf.backward_search(3, 2), std::pair(0, 0));
            ASSERT_EQ(leaf.backward_search(4, 2), std::pair(3, 0));
            ASSERT_EQ(leaf.backward_search(6, 2), std::pair(3, 0));
            ASSERT_EQ(leaf.backward_search(7, 2), std::pair(6, 0));
            ASSERT_EQ(leaf.backward_search(9, 2), std::pair(6, 0));
            ASSERT_EQ(leaf.backward_search(13, 2), std::pair(0, 0));
        }
        {
            Leaf left, right;
            // ab cd ef
            // (()(()((
            for (bool b: {0, 0, 1, 0, 0, 1, 0, 0}) {
                auto [overflow, _] = left.insert(left.size(), b);
                ASSERT_FALSE(overflow);
            }
            //  g   h
            // )()))())
            for (bool b: {1, 0, 1, 1, 1, 0, 1, 1}) {
                auto [overflow, _] = right.insert(right.size(), b);
                ASSERT_FALSE(overflow);
            }


            ASSERT_EQ(left.forward_search(-1, 0), std::pair(-1, 0));
            //ASSERT_EQ(left.forward_search(-1, -1), std::pair(0, 0));

            // excess(j+1) - excess(i) = d

            // d' = - (excess(8) - excess(i) - d)
            // excess(j+1) - excess(8) = - (excess(8) - excess(i) - d)
            // excess(j+1) - excess(8) = - excess(8) + excess(i) + d
            // excess(j+1) - excess(i) = d

            // closing parenthesis
            ASSERT_EQ(left.forward_search(0, 0), std::pair(8, -4));
            ASSERT_EQ(right.forward_search(-1, -4), std::pair(7, 0));
            ASSERT_EQ(left.forward_search(1, 0), std::pair(2, 0));
            ASSERT_EQ(left.forward_search(3, 0), std::pair(8, -3));
            ASSERT_EQ(right.forward_search(-1, -3), std::pair(12 - 8, 0));
            ASSERT_EQ(left.forward_search(4, 0), std::pair(5, 0));
            ASSERT_EQ(left.forward_search(6, 0), std::pair(8, -2));
            ASSERT_EQ(right.forward_search(-1, -2), std::pair(11 - 8, 0));
            ASSERT_EQ(left.forward_search(7, 0), std::pair(8, -1));
            ASSERT_EQ(right.forward_search(-1, -1), std::pair(8 - 8, 0));
            ASSERT_EQ(right.forward_search(1, 0), std::pair(10 - 8, 0));
            ASSERT_EQ(right.forward_search(13 - 8, 0), std::pair(14 - 8, 0));

            // opening parenthesis
            ASSERT_EQ(left.backward_search(2, 0), std::pair(1, 0));
            ASSERT_EQ(left.backward_search(5, 0), std::pair(4, 0));
            ASSERT_EQ(right.backward_search(8 - 8, 0), std::pair(-1, 1));
            ASSERT_EQ(left.backward_search(8, 1), std::pair(7, 0));
            ASSERT_EQ(right.backward_search(2, 0), std::pair(1, 0));
            ASSERT_EQ(right.backward_search(11 - 8, 0), std::pair(-1, 2));
            ASSERT_EQ(left.backward_search(8, 2), std::pair(6, 0));
            ASSERT_EQ(right.backward_search(12 - 8, 0), std::pair(-1, 3));
            ASSERT_EQ(left.backward_search(8, 3), std::pair(3, 0));
            ASSERT_EQ(right.backward_search(14 - 8, 0), std::pair(5, 0));
            ASSERT_EQ(right.backward_search(15 - 8, 0), std::pair(-1, 4));
            ASSERT_EQ(left.backward_search(8, 4), std::pair(0, 0));

            // enclosing parenthesis
            ASSERT_EQ(left.backward_search(1, 2), std::pair(0, 0));
            ASSERT_EQ(left.backward_search(3, 2), std::pair(0, 0));
            ASSERT_EQ(left.backward_search(4, 2), std::pair(3, 0));
            ASSERT_EQ(left.backward_search(6, 2), std::pair(3, 0));
            ASSERT_EQ(left.backward_search(7, 2), std::pair(6, 0));
            ASSERT_EQ(right.backward_search(9 - 8, 2), std::pair(-1, 2));
            ASSERT_EQ(left.backward_search(8, 2), std::pair(6, 0));
            ASSERT_EQ(right.backward_search(13 - 8, 2), std::pair(-1, 4));
            ASSERT_EQ(left.backward_search(8, 4), std::pair(0, 0));
        }
    }
}