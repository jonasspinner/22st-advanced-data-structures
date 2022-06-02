#ifndef ADVANCED_DATA_STRUCTURES_BP_H
#define ADVANCED_DATA_STRUCTURES_BP_H

#include <cstdint>

#include "DynamicBitVector.h"

namespace ads {
    class BP {
    public:
        using size_type = int32_t;
        constexpr static bool open = false;
        constexpr static bool close = true;
    private:

        size_type open_parenthesis_for_node(size_type v) {
            return bv.select(v + 1, false);
        }

        size_type close_parenthesis_for_node(size_type v) {
            // TODO
            assert(false);
        }

        size_type find_close(size_type i) {
            // TODO
            assert(false);
        }

        size_type find_open(size_type i) {
            // TODO
            assert(false);
        }

        size_type excess(size_type i) {
            return bv.rank(i, open) - bv.rank(i, close);
        }

        size_type enclose(size_type i) {
            // TODO
            assert(false);
        }

    public:
        BP() : bv({false, true}) {}

        void insertchild(size_type v, size_type i, size_type k) {
            // TODO
            assert(false);
        }

        void deletenode(size_type v) {
            assert(v != 0);
            // TODO
            assert(false);
        }

        size_type ith_child(size_type v, size_type i) {
            // TODO
            assert(false);
        }

        size_type parent(size_type v) {
            return enclose(v);
        }

        size_type subtree_size(size_type v) {
            return (find_close(v) - v + 1) / 2;
        }

        template<class F>
        void preorder_degree_sequence(F f) {

        }

    private:
        FRIEND_TEST(BPTest, Lecture03);

        FRIEND_TEST(BPTest, ExerciseSheet);

        DynamicBitVector<> bv;
    };
}

#endif //ADVANCED_DATA_STRUCTURES_BP_H
