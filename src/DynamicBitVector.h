#ifndef ADVANCED_DATA_STRUCTURES_DYNAMICBITVECTOR_H
#define ADVANCED_DATA_STRUCTURES_DYNAMICBITVECTOR_H

#include <cstdint>
#include <array>
#include <cassert>
#include <bit>
#include <limits>
#include <tuple>

#include <gtest/gtest.h>

#include "TaggedPointer.h"
#include "utils.h"
#include "SmallStaticBitVector.h"
#include "SmallDynamicBitVector.h"

namespace ads {
    template<class BlockType = uint64_t, size_t NumBlocks = 64>
    class DynamicBitVector {
        static_assert(std::is_unsigned_v<BlockType>);
    public:
        using size_type = int32_t;
        using height_type = int32_t;
    public:
        struct Inner;
        struct Leaf;

        // ## Inner and Leaf Node ##

        class NodeHandle {
            explicit NodeHandle(TaggedPointer<Inner, Leaf> ptr) : m_ptr(ptr) {}

        public:
            [[nodiscard]] constexpr static auto new_leaf() { return leaf(new Leaf); }

            [[nodiscard]] constexpr static auto inner(Inner *ptr) {
                return NodeHandle(TaggedPointer<Inner, Leaf>::first(ptr));
            }

            [[nodiscard]] constexpr static auto leaf(Leaf *ptr) {
                return NodeHandle(TaggedPointer<Inner, Leaf>::second(ptr));
            }

            [[nodiscard]] constexpr std::pair<Inner *, Leaf *> cast() { return m_ptr.cast(); }

            [[nodiscard]] constexpr std::pair<const Inner *, const Leaf *> cast() const { return m_ptr.cast(); }

            [[nodiscard]] constexpr size_type height() const {
                auto [inner, leaf] = m_ptr.cast();
                return inner ? inner->height : 0;
            }

        private:
            TaggedPointer<Inner, Leaf> m_ptr;
        };

        struct Inner {
            NodeHandle left;
            NodeHandle right;
            size_type left_size{};
            size_type left_ones{};
            height_type height{};
        };

        struct Leaf {
            using block_type = BlockType;
            constexpr static size_type max_num_blocks = NumBlocks;

            static_assert(max_num_blocks > 1);
            using Bits = SmallDynamicBitVector<block_type, max_num_blocks>;
            constexpr static size_type block_width = std::numeric_limits<block_type>::digits; // Number of bits per block.

            // Split leave if it has more that `max_num_bits` and merge two leaves if both have `min_num_bits`.
            constexpr static size_type max_num_bits = block_width * max_num_blocks;
            constexpr static size_type min_num_bits = max_num_bits / 4;

            Bits m_bits{};

            [[nodiscard]] size_type size() const { return m_bits.size(); }

            [[nodiscard]] size_type required_bits_upperbound() const {
                return m_bits.required_bits_upperbound();
            }

            [[nodiscard]] bool access(size_type i) const { return m_bits.access(i); }

            void flip(size_type i) { m_bits.flip(i); }

            [[nodiscard]] size_type rank1(size_type i) const { return m_bits.rank1(i); }

            [[nodiscard]] size_type select1(size_type i) const { return m_bits.select1(i); }

            [[nodiscard]] size_type select0(size_type i) const { return m_bits.select0(i); }

            bool remove(size_type i) { return m_bits.remove(i); }

            bool pop_front() { return m_bits.pop_front(); }

            bool pop_back() { return m_bits.pop_back(); }

            [[nodiscard]] std::pair<bool, bool> insert(size_type i, bool b) { return m_bits.insert(i, b); }

            [[nodiscard]] std::pair<bool, bool> push_front(bool b) { return m_bits.push_front(b); }

            [[nodiscard]] bool push_back(bool b) { return m_bits.push_back(b); }

            void split_into_empty(Leaf &other) {
                m_bits.split_into_empty(other.m_bits);
                assert(std::get<0>(check_integrity(*this)));
                assert(std::get<0>(check_integrity(other)));
            }

            friend std::ostream &operator<<(std::ostream &os, const Leaf &leaf) {
                os << "Leaf(size=" << leaf.size() << ", bits=";
                for (size_type i = 0; i < leaf.size(); ++i) {
                    os << (int) leaf.access(i);
                }
                return os << ")";
            }
        };

    private:
        // ## Query operations on nodes ##
        // bool access(node, i)
        // bool access(inner, i)
        // void flip(node, i)
        // void flip(inner, i)
        // void rank1(node, i, rank&)
        // void rank1(inner, i, rank&)
        // void select1(node, i, result&)
        // void select1(inner, i, result&)
        // void select0(node, i, result&)
        // void select0(inner, i, result&)

        static bool access(const NodeHandle &node, size_type i) {
            // Dispatch dependent on node type.
            auto [inner, leaf] = node.cast();
            return inner ? access(*inner, i) : leaf->access(i);
        }

        static bool access(const Inner &inner, size_type i) {
            // If index i is the left subtree, go there, otherwise adjust i and go to the right subtree.
            return i < inner.left_size ? access(inner.left, i) : access(inner.right, i - inner.left_size);
        }

        static void flip(NodeHandle node, size_type i) {
            // Dispatch dependent on node type.
            auto [inner, leaf] = node.cast();
            inner ? flip(*inner, i) : leaf->flip(i);
        }

        static void flip(Inner &inner, size_type i) {
            // If index i is the left subtree, go there, otherwise adjust i and go to the right subtree.
            i < inner.left_size ? flip(inner.left, i) : flip(inner.right, i - inner.left_size);
        }

        static void rank1(const NodeHandle &node, size_type i, size_type &rank) {
            // Dispatch dependent on node type.
            // Update rank via reference to allow for tail call optimization
            auto [inner, leaf] = node.cast();
            if (inner) {
                return rank1(*inner, i, rank);
            } else {
                rank += leaf->rank1(i);
            }
        }

        static void rank1(const Inner &inner, size_type i, size_type &rank) {
            // If index i is the left subtree, go there.
            if (i < inner.left_size) return rank1(inner.left, i, rank);
            // Otherwise, inner.left_ones ones are in the left subtree. Adjust i and go to the right subtree.
            rank += inner.left_ones;
            return rank1(inner.right, i - inner.left_size, rank);
        }

        static void select1(const NodeHandle &node, size_type i, size_type &result) {
            assert(i > 0);
            // Dispatch dependent on node type.
            // Update rank via reference to allow for tail call optimization
            auto [inner, leaf] = node.cast();
            if (inner) {
                return select1(*inner, i, result);
            } else {
                result += leaf->select1(i);
            }
        }

        static void select1(const Inner &inner, size_type i, size_type &result) {
            // If `i`th one bit is the left subtree, go there.
            if (i <= inner.left_ones) return select1(inner.left, i, result);
            // Otherwise, inner.left_size bits are in the left subtree. Adjust i and go to the right subtree.
            result += inner.left_size;
            return select1(inner.right, i - inner.left_ones, result);
        }

        static void select0(const NodeHandle &node, size_type i, size_type &result) {
            assert(i > 0);
            auto [inner, leaf] = node.cast();
            if (inner) {
                return select0(*inner, i, result);
            } else {
                result += leaf->select0(i);
            }
        }

        static void select0(const Inner &inner, size_type i, size_type &result) {
            auto num_zeros_on_left = inner.left_size - inner.left_ones;
            if (i <= num_zeros_on_left) return select0(inner.left, i, result);
            result += inner.left_size;
            return select0(inner.right, i - num_zeros_on_left, result);
        }

        // ## Tree modifying operations ##
        // Inner* split(leaf, high_bit)
        // {node, changed} insert(node, i, b)
        // {node, deleted_bit} remove(node, i)
        // {node, deleted_bit} remove_second_level(inner, i)
        // {node, deleted_bit} remove_first_level(inner, i)
        // Leaf* merge(inner)

        /**
         * Split up one full leaf into two leaves of equal size and append `high_bit` to the right one.
         * @param leaf
         * @param high_bit
         * @return the new inner node that is parent to the two leaves
         */
        static Inner *split(Leaf *leaf, bool high_bit) {
            assert(leaf->size() == Leaf::max_num_bits);
            static_assert(Leaf::max_num_blocks > 1);

            Leaf *left = leaf;
            Leaf *right = new Leaf;

            left->split_into_empty(*right);
            [[maybe_unused]] auto overflow_ = right->push_back(high_bit);

            assert(!overflow_);
            assert(left->size() == Leaf::max_num_bits / 2);
            assert(right->size() == Leaf::max_num_bits / 2 + 1);
            assert(left->size() >= Leaf::min_num_bits);
            assert(right->size() >= Leaf::min_num_bits);

            auto left_size = left->size();
            auto left_ones = left->m_bits.num_ones();
            assert(left_ones == left->rank1(left->size()));

            auto *new_inner = new Inner{
                    NodeHandle::leaf(left), NodeHandle::leaf(right),
                    left_size, left_ones,
                    1};
            assert(std::get<0>(check_integrity(*new_inner)));
            return new_inner;
        }

        /**
         *
         * @param node the subtree root
         * @param i
         * @param b
         * @return {root, changed} the (maybe updated) subtree root and whether the subtree height has changed
         */
        static std::pair<NodeHandle, bool> insert(NodeHandle node, size_type i, bool b) {
            auto [inner, leaf] = node.cast();
            if (leaf) {
                assert(i <= leaf->size());
                // NOTE: i can be the size of the leaf, if it is the right most leaf.

                // Insert b at position i
                auto [overflow, high_bit] = leaf->insert(i, b);
                if (overflow) {
                    // Split leaf if needed
                    node = NodeHandle::inner(split(leaf, high_bit));
                    assert(std::get<0>(check_integrity(node)));
                    return {node, true};
                } else {
                    assert(std::get<0>(check_integrity(node)));
                    return {node, false};
                }
            } else if (i < inner->left_size) {
                bool changed{};
                std::tie(inner->left, changed) = insert(inner->left, i, b);

                // Update statistics
                inner->left_size++;
                inner->left_ones += b;

                // If the height of the subtree changed, recalculate height and if it changed, balance the tree if needed.
                if (changed) {
                    auto h = std::max<int>(inner->height, height(inner->left) + 1);
                    if (inner->height != h) {
                        inner->height = h;

                        assert(std::get<0>(check_integrity(*inner)));
                        return {NodeHandle::inner(balance(inner)), true};
                    }
                }
                return {node, false};
            } else {
                i -= inner->left_size;

                bool changed{};
                std::tie(inner->right, changed) = insert(inner->right, i, b);

                // If the height of the subtree changed, recalculate height and if it changed, balance the tree if needed.
                if (changed) {
                    auto h = std::max(inner->height, height(inner->right) + 1);
                    if (inner->height != h) {
                        inner->height = h;
                        assert(std::get<0>(check_integrity(*inner)));
                        return {NodeHandle::inner(balance(inner)), true};
                    }
                }
                return {node, false};
            }
        }


        /**
         *
         * @param node
         * @param i
         * @return {root, deleted_bit} the (maybe updated) subtree root and the deleted bit
         */
        static std::pair<NodeHandle, bool> remove(NodeHandle node, size_type i) {
            auto [inner, leaf] = node.cast();
            if (leaf) {
                // NOTE: This only applies if the whole tree is only one leaf node.
                //       Otherwise, the cases with height == 1 or height == 2 apply.
                auto deleted_bit = leaf->remove(i);
                return {node, deleted_bit};
            }
            // Special cases for height 1 and 2 inner nodes.
            if (inner->height == 1) { return remove_first_level(inner, i); }
            if (inner->height == 2) { return remove_second_level(inner, i); }

            bool deleted_bit{};
            if (i < inner->left_size) {
                std::tie(inner->left, deleted_bit) = remove(inner->left, i);

                // Update statistics.
                inner->left_size--;
                inner->left_ones -= deleted_bit;
            } else {
                std::tie(inner->right, deleted_bit) = remove(inner->right, i - inner->left_size);
            }

            // The height might have decreased by one.
            inner->height = std::max(height(inner->left), height(inner->right)) + 1;
            return {NodeHandle::inner(balance(inner)), deleted_bit};
        }

        /**
         *
         * @param a
         * @param i
         * @return {root, deleted_bit} the (maybe updated) subtree root and the deleted bit
         */
        static std::pair<NodeHandle, bool> remove_second_level(Inner *a, size_type i) {
            assert(std::get<0>(check_integrity(*a)));
            assert(a->height == 2);

            // Following cases are distinguished:
            // 1.    i is in the left subtree
            // 1.1     and the left child is an inner node.
            // 1.2     and the left child is a leaf
            // 1.2.1     and the leaf only has the minimum number of bits left.
            // 1.2.2     and the leaf has more than the minimum number of bits.
            // 2.    i is in the right subtree
            // 2.1     and the right child is an inner node.
            // 2.2     and the right child is a leaf
            // 2.2.1     and the leaf only has the minimum number of bits left.
            // 2.2.2     and the leaf has more than the minimum number of bits.

            if (i < a->left_size) {
                // 1.    i is in the left subtree
                //   a
                // b
                auto [b_inner, b_leaf] = a->left.cast();
                if (b_inner) {
                    // 1.1     and the left child is an inner node.
                    // Strategy: Use level 1 algorithm on left child.
                    bool deleted_bit{};
                    std::tie(a->left, deleted_bit) = remove_first_level(b_inner, i);

                    a->left_size--;
                    a->left_ones -= deleted_bit;
                    a->height = std::max(height(a->left), height(a->right)) + 1;

                    assert(std::get<0>(check_integrity(*a)));
                    return {NodeHandle::inner(balance(a)), deleted_bit};
                } else {
                    // 1.2     and the left child is a leaf
                    assert(b_leaf);
                    if (b_leaf->size() == Leaf::min_num_bits) {
                        // 1.2.1     and the leaf only has the minimum number of bits left.
                        //    a
                        //  b   c
                        //    d
                        // Strategy: Delete bit from left child
                        //           and append the first bit of the right subtree to the left subtree.
                        auto c = a->right.cast().first;
                        assert(c);
                        auto deleted_bit = b_leaf->remove(i);
                        bool moved_bit{};
                        std::tie(a->right, moved_bit) = remove_first_level(c, 0);

                        // Append moved bit to b. This will not overflow, as b has min_num_bits - 1 bits.
                        [[maybe_unused]] auto overflow_ = b_leaf->push_back(moved_bit);
                        assert(!overflow_);

                        a->left_ones -= deleted_bit;
                        a->left_ones += moved_bit;
                        a->height = std::max(height(a->left), height(a->right)) + 1;

                        assert(std::get<0>(check_integrity(*a)));
                        assert(b_leaf->size() >= Leaf::min_num_bits);
                        return {NodeHandle::inner(balance(a)), deleted_bit};
                    } else {
                        // 1.2.2     and the leaf has more than the minimum number of bits.
                        // Strategy: Directly delete bit from the left child.
                        assert(b_leaf->size() > Leaf::min_num_bits);
                        auto deleted_bit = b_leaf->remove(i);

                        a->left_size--;
                        a->left_ones -= deleted_bit;

                        assert(std::get<0>(check_integrity(*a)));
                        assert(b_leaf->size() >= Leaf::min_num_bits);
                        return {NodeHandle::inner(a), deleted_bit};
                    }
                }
            } else {
                // 2.    i is in the right subtree
                // a
                //   b

                i -= a->left_size;

                auto [b_inner, b_leaf] = a->right.cast();
                if (b_inner) {
                    // 2.1     and the right child is an inner node.
                    // Strategy: Use level 1 algorithm on right child.
                    bool deleted_bit{};
                    std::tie(a->right, deleted_bit) = remove_first_level(b_inner, i);

                    a->height = std::max<int>(height(a->left), height(a->right)) + 1;

                    assert(std::get<0>(check_integrity(*a)));
                    return {NodeHandle::inner(balance(a)), deleted_bit};
                } else {
                    // 2.2     and the right child is a leaf

                    assert(b_leaf);
                    assert(a->left_size == (a->left.cast().first->left.cast().second->size() +
                                            a->left.cast().first->right.cast().second->size()));
                    if (b_leaf->size() == Leaf::min_num_bits) {
                        // 2.2.1     and the leaf only has the minimum number of bits left.
                        //    a
                        //  c   b
                        //    d
                        // Strategy: Delete bit from right child
                        //           and prepend the last bit of the left subtree to the right subtree.
                        auto c = a->left.cast().first;
                        assert(c);
                        auto deleted_bit = b_leaf->remove(i);
                        bool moved_bit{};
                        std::tie(a->left, moved_bit) = remove_first_level(c, a->left_size - 1);
                        [[maybe_unused]] auto [overflow_, high_bit_] = b_leaf->push_front(moved_bit);
                        assert(!overflow_);

                        a->left_size--;
                        a->left_ones -= moved_bit;
                        a->height = std::max(height(a->left), height(a->right)) + 1;

                        assert(std::get<0>(check_integrity(*a)));
                        assert(b_leaf->size() >= Leaf::min_num_bits);
                        return {NodeHandle::inner(balance(a)), deleted_bit};
                    } else {
                        // 2.2.2     and the leaf has more than the minimum number of bits.
                        // Strategy: Directly delete bit from the right child.
                        assert(b_leaf->size() > Leaf::min_num_bits);

                        auto deleted_bit = b_leaf->remove(i);

                        assert(std::get<0>(check_integrity(*a)));
                        assert(b_leaf->size() >= Leaf::min_num_bits);
                        return {NodeHandle::inner(a), deleted_bit};
                    }
                }
            }
            assert(false);
        }

        static std::pair<NodeHandle, bool> remove_first_level(Inner *a, size_type i) {
            assert(std::get<0>(check_integrity(*a)));
            assert(a->height == 1);
            //   a
            // b   c
            auto b = a->left.cast().second;
            auto c = a->right.cast().second;
            assert(b && c);

            // Following cases are distinguished:
            // 1.    i is in the left child
            // 1.1     and the left child only has the minimum number of bits left
            // 1.1.1     and the right child only has the minimum number of bits left.
            // 1.1.2     and the right child has more than the minimum number of bits left.
            // 1.2     and the left child has more than the minimum number of bits left.
            // 2.    i is in the right child
            // 2.1     and the right child only has the minimum number of bits left
            // 2.1.1     and the left child only has the minimum number of bits left.
            // 2.1.2     and the left child has more than the minimum number of bits left.
            // 2.2     and the right child has more than the minimum number of bits left.

            if (i < a->left_size) {
                // 1.    i is in the left child

                if (b->size() <= Leaf::min_num_bits) {
                    // 1.1     and the left child only has the minimum number of bits left

                    if (c->size() <= Leaf::min_num_bits) {
                        // 1.1.1     and the right child only has the minimum number of bits left.
                        // Strategy: Merge both children and remove the bit in the resulting leaf.

                        auto *leaf = merge(a);
                        auto deleted_bit = leaf->remove(i);
                        assert(std::get<0>(check_integrity(*leaf)));
                        return {NodeHandle::leaf(leaf), deleted_bit};
                    } else {
                        // 1.1.2     and the right child has more than the minimum number of bits left.
                        // Strategy: Delete bit from the left child
                        //           and append the first bit of the right child to the left child.

                        if constexpr(false && Leaf::max_num_bits >= Leaf::min_num_bits + 2 * Leaf::block_width) {
                            if (c->size() >= Leaf::min_num_bits + 2 * Leaf::block_width) {
                                size_type num_moved_blocks{};
                                while (c->size() >= b->size() + Leaf::block_width) {
                                    auto block = c->m_bits.pop_front_block();
                                    a->left_size += Leaf::block_width;
                                    a->left_ones += std::popcount(block);
                                    b->m_bits.push_back_block_aligned(block);
                                    num_moved_blocks++;
                                }
                                auto deleted_bit = b->remove(i);
                                a->left_size--;
                                a->left_ones -= deleted_bit;

                                assert(std::get<0>(check_integrity(*a)));
                                return {NodeHandle::inner(a), deleted_bit};
                            }
                        }

                        auto deleted_bit = b->remove(i);
                        auto moved_bit = c->pop_front();
                        [[maybe_unused]] auto overflow_ = b->push_back(moved_bit);
                        assert(!overflow_);

                        a->left_ones -= deleted_bit;
                        a->left_ones += moved_bit;

                        assert(std::get<0>(check_integrity(*a)));
                        // assert(b->size() >= Leaf::min_num_bits && c->size() >= Leaf::min_num_bits);
                        return {NodeHandle::inner(a), deleted_bit};
                    }
                } else {
                    // 1.2     and the left child has more than the minimum number of bits left.
                    // Strategy: Directly delete bit from left child.

                    assert(b->size() > Leaf::min_num_bits);

                    auto deleted_bit = b->remove(i);

                    a->left_size--;
                    a->left_ones -= deleted_bit;

                    assert(std::get<0>(check_integrity(*a)));
                    return {NodeHandle::inner(a), deleted_bit};
                }
            } else {
                // 2.    i is in the right child

                if (c->size() <= Leaf::min_num_bits) {
                    // 2.1     and the right child only has the minimum number of bits left

                    if (b->size() <= Leaf::min_num_bits) {
                        // 2.1.1     and the left child only has the minimum number of bits left.
                        // Strategy: Merge both children and remove the bit in the resulting leaf.

                        auto *leaf = merge(a);
                        auto deleted_bit = leaf->remove(i);

                        assert(std::get<0>(check_integrity(*leaf)));
                        return {NodeHandle::leaf(leaf), deleted_bit};
                    } else {
                        // 2.1.2     and the left child has more than the minimum number of bits left.
                        // Strategy: Delete bit from the right child
                        //           and prepend the last bit of the left child to the right child.

                        assert(b->size() > Leaf::min_num_bits);
                        if constexpr (false && Leaf::max_num_bits >= Leaf::min_num_bits + 2 * Leaf::block_width) {
                            if (b->size() >= Leaf::min_num_bits + 2 * Leaf::block_width &&
                                b->size() % Leaf::block_width == 0) {
                                //std::cout << "could move block: " << b->size() << " " << c->size() << std::endl;
                                i -= a->left_size;
                                size_type num_moved_bits{};
                                while (b->size() % Leaf::block_width != 0) {
                                    auto moved_bit = b->pop_back();
                                    [[maybe_unused]] auto [overflow_, high_bit_] = c->push_front(moved_bit);
                                    a->left_size--;
                                    a->left_ones -= moved_bit;
                                    ++i;
                                    num_moved_bits++;
                                }
                                assert(b->size() % Leaf::block_width == 0);

                                size_type num_moved_blocks{};
                                while (b->size() >= c->size() + Leaf::block_width) {
                                    auto block = b->m_bits.pop_back_block_aligned();
                                    a->left_size -= Leaf::block_width;
                                    a->left_ones -= std::popcount(block);
                                    c->m_bits.push_front_block(block);
                                    i += Leaf::block_width;
                                    num_moved_blocks++;
                                }
                                auto deleted_bit = c->remove(i);

                                //std::cout << "moved " << num_moved_bits << " bits and " << num_moved_blocks << " blocks"
                                //          << std::endl;
                                //std::cout << "      " << b->size() << " " << c->size() << std::endl;

                                assert(std::get<0>(check_integrity(*a)));
                                assert(b->size() >= Leaf::min_num_bits && c->size() >= Leaf::min_num_bits);
                                return {NodeHandle::inner(a), deleted_bit};

                            }
                        }

                        auto deleted_bit = c->remove(i - a->left_size);
                        auto moved_bit = b->pop_back();
                        [[maybe_unused]] auto [overflow_, high_bit_] = c->push_front(moved_bit);
                        assert(!overflow_);

                        a->left_size--;
                        a->left_ones -= moved_bit;

                        assert(std::get<0>(check_integrity(*a)));
                        return {NodeHandle::inner(a), deleted_bit};
                    }
                } else {
                    // 2.2     and the right child has more than the minimum number of bits left.
                    // Strategy: Directly delete bit from right child.

                    auto deleted_bit = c->remove(i - a->left_size);

                    assert(std::get<0>(check_integrity(*a)));
                    return {NodeHandle::inner(a), deleted_bit};
                }
            }
            assert(false);
        }

        /**
         * Merge the two children of a height 1 inner node into one leaf node. Deletes the inner node and one of the leaf nodes.
         * Both leaf nodes must have the minimum number of bits.
         * @param a an inner node with two leaves as children
         * @return a leaf node with the bits of both leaves combined
         */
        static Leaf *merge(Inner *a) {
            assert(a->height == 1);

            Leaf *b = a->left.cast().second;
            Leaf *c = a->right.cast().second;
            assert(b && c);

            if (b->size() == Leaf::min_num_bits && c->size() == Leaf::min_num_bits) {
                if constexpr (Leaf::min_num_bits % Leaf::block_width == 0) {
                    // If the minimum number of bits aligns with block boundaries, we can copy whole blocks.
                    b->m_bits.concat_block_aligned(c->m_bits);
                } else {
                    for (size_type j = 0; j < Leaf::min_num_bits; ++j) {
                        [[maybe_unused]] auto overflow_ = b->push_back(c->access(j));
                        assert(!overflow_);
                    }
                }
                assert(b->size() == 2 * Leaf::min_num_bits);
            } else {
                for (int i = 0; i < c->size(); ++i) {
                    [[maybe_unused]] auto overflow_ = b->push_back(c->access(i));
                    assert(!overflow_);
                }
            }
            delete a;
            delete c;

            assert(std::get<0>(check_integrity(*b)));
            return b;
        }

        // ### Node info dispatch ###
        // int height(node)

        static int height(const NodeHandle &node) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                assert(inner->height == std::max(height(inner->left), height(inner->right)) + 1);
                return inner->height;
            } else { return 0; }
        }

#ifndef NDEBUG

        static int size_(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? size_(inner->left) + size_(inner->right) : leaf->size();
        }

        static int ones_(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? ones_(inner->left) + ones_(inner->right) : leaf->rank1(leaf->size());
        }

#endif

        // ## AVL Tree balance operations ##
        // int balance_factor(node)
        // int balance_factor(inner)
        // Inner* balance(inner)
        // Inner* rotate_left_left(inner)
        // Inner* rotate_left_right(inner)
        // Inner* rotate_right_left(inner)
        // Inner* rotate_right_right(inner)

        static int balance_factor(const NodeHandle &node) {
            auto [inner, leaf] = node.cast();
            return inner ? balance_factor(*inner) : 0;
        }

        static int balance_factor(const Inner &node) {
            return height(node.left) - height(node.right);
        }

        static Inner *balance(Inner *node) {
            auto factor = balance_factor(*node);
            if (factor > 1) {
                if (balance_factor(node->left) > 0) {
                    return rotate_left_left(node);
                } else {
                    return rotate_left_right(node);
                }
            } else if (factor < -1) {
                if (balance_factor(node->right) > 0) {
                    return rotate_right_left(node);
                } else {
                    return rotate_right_right(node);
                }
            }
            return node;
        }

        static Inner *rotate_left_left(Inner *a) {
            //     a            b
            //   b       ->   c   a
            // c
            assert(a);
            Inner *b = a->left.cast().first;
            assert(b);

            a->left = b->right;
            b->right = NodeHandle::inner(a);

            a->left_size -= b->left_size;
            a->left_ones -= b->left_ones;

            a->height = std::max(height(a->left), height(a->right)) + 1;
            b->height = std::max(height(b->left), a->height) + 1;

            assert(std::get<0>(check_integrity(NodeHandle::inner(b))));
            return b;
        }

        static Inner *rotate_left_right(Inner *a) {
            //   a         c
            // b    ->   b   a
            //   c
            assert(a);
            Inner *b = a->left.cast().first;
            assert(b);
            Inner *c = b->right.cast().first;
            assert(c);

            b->right = c->left;
            a->left = c->right;
            c->left = NodeHandle::inner(b);
            c->right = NodeHandle::inner(a);

            a->left_size -= b->left_size + c->left_size;
            a->left_ones -= b->left_ones + c->left_ones;
            c->left_size += b->left_size;
            c->left_ones += b->left_ones;
            a->height = std::max(height(a->left), height(a->right)) + 1;
            b->height = std::max(height(b->left), height(b->right)) + 1;
            c->height = std::max(a->height, b->height) + 1;

            assert(std::get<0>(check_integrity(*c)));
            return c;
        }

        static Inner *rotate_right_left(Inner *a) {
            //  a          c
            //    b  ->  a   b
            //  c
            assert(a);
            Inner *b = a->right.cast().first;
            assert(b);
            Inner *c = b->left.cast().first;
            assert(c);
            a->right = c->left;
            b->left = c->right;
            c->left = NodeHandle::inner(a);
            c->right = NodeHandle::inner(b);

            b->left_size -= c->left_size;
            b->left_ones -= c->left_ones;
            c->left_size += a->left_size;
            c->left_ones += a->left_ones;
            a->height = std::max(height(a->left), height(a->right)) + 1;
            b->height = std::max(height(b->left), height(b->right)) + 1;
            c->height = std::max(a->height, b->height) + 1;

            assert(std::get<0>(check_integrity(*c)));
            return c;
        }

        static Inner *rotate_right_right(Inner *a) {
            //  a             b
            //    b    ->   a   c
            //      c
            assert(a);
            Inner *b = a->right.cast().first;
            assert(b);
            a->right = b->left;
            b->left = NodeHandle::inner(a);

            b->left_size += a->left_size;
            b->left_ones += a->left_ones;
            a->height = std::max(height(a->left), height(a->right)) + 1;
            b->height = std::max<int>(a->height, height(b->right)) + 1;

            assert(std::get<0>(check_integrity(*b)));
            return b;
        }

        static void delete_node(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? delete_node(inner) : delete_node(leaf);
        }

        static void delete_node(Inner *inner) {
            delete_node(inner->left);
            delete_node(inner->right);
            delete inner;
        }

        static void delete_node(Leaf *leaf) { delete leaf; }

        // ## Debug checks ##
#ifndef NDEBUG

        static std::tuple<bool, size_type, size_type, size_type>
        check_integrity(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            if (inner) { return check_integrity(*inner); }
            if (leaf) { return check_integrity(*leaf); }
            return {false, 0, 0, 0};
        }

        static std::tuple<bool, size_type, size_type, size_type>
        check_integrity(const Inner &inner) {
            auto [left_ok, left_size, left_ones, left_height] =
                    check_integrity(inner.left);
            auto [right_ok, right_size, right_ones, right_height] =
                    check_integrity(inner.right);
            auto height = std::max(left_height, right_height) + 1;
            bool ok = left_ok && right_ok
                      && (inner.left_size == left_size)
                      && (inner.left_ones == left_ones)
                      && (inner.height == height);
            assert(ok);
            return {ok, left_size + right_size, left_ones + right_ones, height};
        }

        static std::tuple<bool, size_type, size_type, size_type>
        check_integrity(const Leaf &leaf) {
            auto size = leaf.size();
            auto ones = leaf.rank1(leaf.size());
            bool ok = true;
            return {ok, size, ones, 0};
        }

#endif

        // ## Required bits upperbound ##

        static size_type required_bits_upperbound(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                return required_bits_upperbound(inner->left) + required_bits_upperbound(inner->right) +
                       sizeof(Inner) * 8;
            } else {
                return leaf->required_bits_upperbound();
            }
        }

    public:
        // ## Constructor/Destructor ##

        DynamicBitVector() : m_root(NodeHandle::leaf(new Leaf)), m_size(0) {}

        explicit DynamicBitVector(const std::vector<bool> &bits) : m_root(NodeHandle::leaf(new Leaf)), m_size(0) {
            for (auto bit: bits) {
                insert(size(), bit);
            }
        }

        ~DynamicBitVector() { delete_node(m_root); }

        // ## BV operations ##
        // ### BV modifying operations ##
        // void insert(i, b)
        // void remove(i)
        // void flip(i)

        void insert(size_type i, bool b) {
            assert(i <= size());
            bool changed{};
            std::tie(m_root, changed) = insert(m_root, i, b);
            m_size++;

#ifndef NDEBUG
            auto [ok, size, ones, height] = check_integrity(m_root);
            assert(ok && (m_size == size));
#endif
        }

        void remove(size_type i) {
            assert(i < size());
            bool deleted_bit{};
            std::tie(m_root, deleted_bit) = remove(m_root, i);
            m_size--;

#ifndef NDEBUG
            auto [ok, size, ones, height] = check_integrity(m_root);
            assert(ok && (m_size == size));
#endif
        }

        void flip(size_type i) const {
            assert(i < size());
            flip(m_root, i);
        }

        // ### BV query operations ##
        // bool access(i)
        // size_type rank(i, b)
        // size_type select(i, b)
        // size_type size()

        [[nodiscard]] bool access(size_type i) const { return access(m_root, i); }

        [[nodiscard]] size_type rank(size_type i, bool b) const {
            assert(i <= size());
            size_type rank = 0;
            rank1(m_root, i, rank);
            return b ? rank : i - rank;
        }

        [[nodiscard]] size_type select(size_type i, bool b) const {
            size_type result = 0;
            if (b) {
                select1(m_root, i, result);
                return result;
            } else {
                select0(m_root, i, result);
                return result;
            }
        }

        [[nodiscard]] size_type size() const {
            return m_size;
        }

        // ## Misc ##
        // clear()
        // std::ostream << overload
        // size_type required_bits_upperbound()

        void clear() {
            delete_node(m_root);
            m_root = new Leaf;
            m_size = 0;
        }

        friend std::ostream &operator<<(std::ostream &os, const DynamicBitVector &bv) {
            print(os, bv.m_root, 0);
            return os;
        }

        [[nodiscard]] size_type required_bits_upperbound() const {
            return required_bits_upperbound(m_root) + (sizeof m_root) * 8;
        }

    private:
        static void print(std::ostream &os, NodeHandle node, size_type indent) {
            auto [inner, leaf] = node.cast();
            for (size_type i = 0; i < indent; ++i) os << ' ';
            if (inner) {
                os << "Inner(num=" << inner->left_size << ", ones=" << inner->left_ones << ", height=" << inner->height
                   << ", balance_factor=" << balance_factor(*inner) << ")\n";
                print(os, inner->left, indent + 1);
                print(os, inner->right, indent + 1);
            } else {
                os << *leaf << "\n";
            }
        }

        FRIEND_TEST(BVTest, Leaf);

        FRIEND_TEST(BVTest, Inner);

        FRIEND_TEST(BVTest, Insert);

        FRIEND_TEST(BVTest, Remove);

        NodeHandle m_root;
        size_type m_size;
    };
}

#endif //ADVANCED_DATA_STRUCTURES_DYNAMICBITVECTOR_H
