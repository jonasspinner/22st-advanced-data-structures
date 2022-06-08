#include <gtest/gtest.h>
#include <random>

#include "BP.h"
#include "BP2.h"

namespace ads {
    TEST(BPTest, Lecture03) {
        using ads::BP;

        /*
        BP bp;
        bp.insertchild(0, 2, 1);
        bp.insertchild(0, 3, 0);
        bp.insertchild(0, 2, 3);
         */
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

    TEST(BP2Test, Inner) {
        using BP2 = ads::BP2<>;
        using Inner = BP2::Inner;
        using Leaf = BP2::Leaf;
        using NodeHandle = BP2::NodeHandle;
        using size_type = BP2::size_type;

        std::array<Leaf, 4> leaves;
        // ab c  d ef g     h   ij k
        // (()(  ()(()(  )))()  (()()))
        std::array<std::vector<bool>, 4> bits = {
                std::vector<bool>{0, 0, 1, 0},
                {0, 1, 0, 0, 1, 0},
                {1, 1, 1, 0, 1},
                {0, 0, 1, 0, 1, 1, 1}
        };
        for (size_t i = 0; i < leaves.size(); ++i) {
            for (bool b: bits[i]) {
                [[maybe_unused]] auto r = leaves[i].insert(leaves[i].size(), b);
            }
            leaves[i].total_excess = leaves[i].calculate_total_excess();
            leaves[i].min_excess = leaves[i].calculate_min_excess();
        }

        Inner b_{NodeHandle::leaf(&leaves[0]), NodeHandle::leaf(&leaves[1]), 4, 1, 4, 0, 1};
        auto b = NodeHandle::inner(&b_);
        Inner c_{NodeHandle::leaf(&leaves[2]), NodeHandle::leaf(&leaves[3]), 5, 4, -4, -4, 1};
        auto c = NodeHandle::inner(&c_);
        Inner a_{b, c, 10, 3, 0, 0, 2};
        auto a = NodeHandle::inner(&a_);

        ASSERT_EQ(BP2::access(a, 3), false);
        ASSERT_EQ(BP2::access(a, 4), false);
        ASSERT_EQ(BP2::access(a, 9), false);
        ASSERT_EQ(BP2::access(a, 10), true);
        ASSERT_EQ(BP2::access(a, 14), true);
        ASSERT_EQ(BP2::access(a, 15), false);
        ASSERT_EQ(BP2::access(a, 21), true);

        //  |               0,0                |
        //  |      4,0     |  |    -4,-4       |
        //  | 2,0|  |  2,0 |  |-3,-3|  | -1,-1 |
        //  |ab c|  |d ef g|  |   h |  |ij k   |
        //   (()(    ()(()(    )))()    (()()))

        //  |    |  |   +-+|  |-    |  |       |
        //  |    |  |+-+   |  | -   |  | +-+-  |
        //  | +-+|  |      |  |  -+-|  |+    - |
        //  |+   |  |      |  |     |  |      -|

        ASSERT_EQ(leaves[0].forward_search(0, 0), std::pair(4, -2));
        ASSERT_EQ(leaves[1].forward_search(-1, -2), std::pair(6, -4));
        ASSERT_EQ(leaves[2].forward_search(-1, -4), std::pair(5, -1));
        ASSERT_EQ(leaves[3].forward_search(-1, -1), std::pair(6, 0));

        ASSERT_EQ(BP2::forward_search(b, 0, 0), std::pair(BP2::right_end, -4));
        ASSERT_EQ(BP2::forward_search(b, 1, 0), std::pair(2, 0));
        ASSERT_EQ(BP2::forward_search(b, 3, 0), std::pair(BP2::right_end, -3));
        ASSERT_EQ(BP2::forward_search(b, 4, 0), std::pair(5, 0));
        ASSERT_EQ(BP2::forward_search(b, 6, 0), std::pair(BP2::right_end, -2));
        ASSERT_EQ(BP2::forward_search(b, 7, 0), std::pair(8, 0));
        ASSERT_EQ(BP2::forward_search(b, 9, 0), std::pair(BP2::right_end, -1));

        ASSERT_EQ(BP2::forward_search(c, 3, 0), std::pair(4, 0));
        ASSERT_EQ(BP2::forward_search(c, 5, 0), std::pair(10, 0));
        ASSERT_EQ(BP2::forward_search(c, 6, 0), std::pair(7, 0));
        ASSERT_EQ(BP2::forward_search(c, 8, 0), std::pair(9, 0));

        ASSERT_EQ(BP2::forward_search(a, 0, 0), std::pair(21, 0));
        ASSERT_EQ(BP2::forward_search(a, 1, 0), std::pair(2, 0));
        ASSERT_EQ(BP2::forward_search(a, 3, 0), std::pair(12, 0));
        ASSERT_EQ(BP2::forward_search(a, 4, 0), std::pair(5, 0));
        ASSERT_EQ(BP2::forward_search(a, 6, 0), std::pair(11, 0));
        ASSERT_EQ(BP2::forward_search(a, 7, 0), std::pair(8, 0));
        ASSERT_EQ(BP2::forward_search(a, 9, 0), std::pair(10, 0));
        ASSERT_EQ(BP2::forward_search(a, 13, 0), std::pair(14, 0));
        ASSERT_EQ(BP2::forward_search(a, 15, 0), std::pair(20, 0));
        ASSERT_EQ(BP2::forward_search(a, 16, 0), std::pair(17, 0));
        ASSERT_EQ(BP2::forward_search(a, 18, 0), std::pair(19, 0));


        //   (()(    ()(()(    )))()    (()()))
        //   0123    456789    01234    5678901
        //   0123    012345    01234    0123456
        ASSERT_EQ(leaves[3].backward_search(6, 0), std::pair(-1, 1));
        ASSERT_EQ(leaves[2].backward_search(5, 1), std::pair(-1, 4));
        ASSERT_EQ(leaves[1].backward_search(6, 4), std::pair(-1, 2));
        ASSERT_EQ(leaves[0].backward_search(4, 2), std::pair(0, 0));

        ASSERT_EQ(BP2::backward_search(b, 2, 0), std::pair(1, 0));
        ASSERT_EQ(BP2::backward_search(b, 5, 0), std::pair(4, 0));
        ASSERT_EQ(BP2::backward_search(b, 8, 0), std::pair(7, 0));
        ASSERT_EQ(BP2::backward_search(b, BP2::right_end, 4), std::pair(0, 0));
        ASSERT_EQ(BP2::backward_search(b, BP2::right_end, 3), std::pair(3, 0));
        ASSERT_EQ(BP2::backward_search(b, BP2::right_end, 2), std::pair(6, 0));
        ASSERT_EQ(BP2::backward_search(b, BP2::right_end, 1), std::pair(9, 0));

        ASSERT_EQ(BP2::backward_search(c, 0, 0), std::pair(BP2::left_end, 1));
        ASSERT_EQ(BP2::backward_search(c, 1, 0), std::pair(BP2::left_end, 2));
        ASSERT_EQ(BP2::backward_search(c, 2, 0), std::pair(BP2::left_end, 3));
        ASSERT_EQ(BP2::backward_search(c, 4, 0), std::pair(3, 0));
        ASSERT_EQ(BP2::backward_search(c, 7, 0), std::pair(6, 0));
        ASSERT_EQ(BP2::backward_search(c, 9, 0), std::pair(8, 0));
        ASSERT_EQ(BP2::backward_search(c, 10, 0), std::pair(5, 0));
        ASSERT_EQ(BP2::backward_search(c, 11, 0), std::pair(BP2::left_end, 4));


        ASSERT_EQ(BP2::backward_search(a, 2, 0), std::pair(1, 0));
        ASSERT_EQ(BP2::backward_search(a, 5, 0), std::pair(4, 0));
        ASSERT_EQ(BP2::backward_search(a, 8, 0), std::pair(7, 0));
        ASSERT_EQ(BP2::backward_search(a, 10, 0), std::pair(9, 0));
        ASSERT_EQ(BP2::backward_search(a, 11, 0), std::pair(6, 0));
        ASSERT_EQ(BP2::backward_search(a, 12, 0), std::pair(3, 0));
        ASSERT_EQ(BP2::backward_search(a, 14, 0), std::pair(13, 0));
        ASSERT_EQ(BP2::backward_search(a, 17, 0), std::pair(16, 0));
        ASSERT_EQ(BP2::backward_search(a, 19, 0), std::pair(18, 0));
        ASSERT_EQ(BP2::backward_search(a, 20, 0), std::pair(15, 0));
        ASSERT_EQ(BP2::backward_search(a, 21, 0), std::pair(0, 0));

        ASSERT_EQ(BP2::backward_search(a, 1, 2), std::pair(0, 0));
        ASSERT_EQ(BP2::backward_search(a, 3, 2), std::pair(0, 0));
        ASSERT_EQ(BP2::backward_search(a, 4, 2), std::pair(3, 0));
        ASSERT_EQ(BP2::backward_search(a, 6, 2), std::pair(3, 0));
        ASSERT_EQ(BP2::backward_search(a, 7, 2), std::pair(6, 0));
        ASSERT_EQ(BP2::backward_search(a, 9, 2), std::pair(6, 0));
        ASSERT_EQ(BP2::backward_search(a, 13, 2), std::pair(0, 0));
        ASSERT_EQ(BP2::backward_search(a, 15, 2), std::pair(0, 0));
        ASSERT_EQ(BP2::backward_search(a, 16, 2), std::pair(15, 0));
        ASSERT_EQ(BP2::backward_search(a, 18, 2), std::pair(15, 0));

        {
            Leaf leaf;
            for (bool b: bits[0]) {
                leaf.insert(leaf.size(), b);
            }
            for (bool b: bits[1]) {
                leaf.insert(leaf.size(), b);
            }

            ASSERT_EQ(leaf.forward_search(0, 0), std::pair(10, -4));
        }
    }

    TEST(BP2Test, BVInsert) {
        using BP2 = ads::BP2<>;
        using size_type = BP2::size_type;

        BP2 bp;

        std::mt19937_64 gen;
        std::bernoulli_distribution bit_dist;
        size_t n = 1000;
        for (size_t i = 0; i < n; ++i) {
            std::uniform_int_distribution<size_type> index_dist(0, bp.size());
            bp.insert(index_dist(gen), bit_dist(gen));
        }
        for (size_t i = 0; i < n; ++i) {
            bp.insert(0, bit_dist(gen));
        }
        for (size_t i = 0; i < n; ++i) {
            bp.insert(bp.size() / 2, bit_dist(gen));
        }
        for (size_t i = 0; i < n; ++i) {
            bp.insert(bp.size(), bit_dist(gen));
        }
        ASSERT_EQ(bp.size(), 4 * n + 2);
    }

}