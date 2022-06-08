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

namespace ads {
    template<class BlockType = uint64_t, size_t NumBlocks = 8>
    class DynamicBitVector {
        static_assert(std::is_unsigned_v<BlockType>);
    public:
        using size_type = int32_t;
        using height_type = int32_t;
    public:
        struct Inner;
        struct Leaf;

        class NodeHandle {
            explicit NodeHandle(TaggedPointer<Inner, Leaf> ptr) : m_ptr(ptr) {}

        public:
            static auto new_leaf() { return leaf(new Leaf); }

            static auto inner(Inner *ptr) { return NodeHandle(TaggedPointer<Inner, Leaf>::first(ptr)); }

            static auto leaf(Leaf *ptr) { return NodeHandle(TaggedPointer<Inner, Leaf>::second(ptr)); }

            std::pair<Inner *, Leaf *> cast() { return m_ptr.cast(); }

            std::pair<const Inner *, const Leaf *> cast() const { return m_ptr.cast(); }

            [[nodiscard]] size_type height() const {
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
            constexpr static size_type num_blocks = NumBlocks;

            static_assert(num_blocks > 1);
            using Bits = std::array<block_type, num_blocks>;
            constexpr static size_type block_width = std::numeric_limits<block_type>::digits; // Number of bits per block.

            // Split leave if it has more that `max_num_bits` and merge two leaves if both have `min_num_bits`.
            constexpr static size_type max_num_bits = block_width * num_blocks;
            constexpr static size_type min_num_bits = max_num_bits / 4;

            Bits m_bits{};
            size_type m_size{};

            /**
             * @return the number of bits stored.
             */
            [[nodiscard]] size_type size() const { return m_size; }

            /**
             * @return the maximum number of bits that can be stored.
             */
            [[nodiscard]] constexpr size_type capacity() const { return m_bits.size() * block_width; }

            [[nodiscard]] bool access(size_type i) const {
                assert(0 <= i);
                assert(i < size());
                size_type j = i / block_width;
                size_type k = i % block_width;
                return (m_bits[j] >> k) & bit_mask(0);
            }

            void flip(size_type i) {
                assert(i < size());
                size_type j = i / block_width;
                size_type k = i % block_width;
                m_bits[j] ^= bit_mask(k);
            }

            [[nodiscard]] size_type rank1(size_type i) const {
                assert(0 <= i);
                assert(i <= size());
                size_type j = i / block_width;
                size_type rank = 0;
                // Sum up number of ones in all blocks before block `j`.
                for (size_type block_idx = 0; block_idx < num_blocks && block_idx < j; ++block_idx) {
                    rank += std::popcount(m_bits[block_idx]);
                    i -= block_width;
                }
                // This covers the case that `i` equals `size` and prevents memory access after the last block
                if (i == 0) return rank;
                assert(i < block_width);
                // Only count number of ones in the bits lower than `i`.
                block_type masked_block = m_bits[j] & lower_bit_mask(i);
                rank += std::popcount(masked_block);
                return rank;
            }

            [[nodiscard]] size_type select1(size_type i) const {
                assert(i > 0);
                size_type j = 0;
                size_type result = 0;
                // Skip blocks until the `i`th ones bit is in block `j`.
                // `result` is the number of bits in the previous blocks.
                // `i` is updated such that it corresponds to the `i`th one bit in block `j`.
                while (j < num_blocks - 1) {
                    size_type count = std::popcount(m_bits[j]);
                    if (i > count) {
                        i -= count;
                        result += block_width;
                        j++;
                    } else { break; }
                }
                if (i == 0) return result;
                assert(j < num_blocks);
                block_type block = m_bits[j];
                return result + ads::select(block, i);
            }

            [[nodiscard]] size_type select0(size_type i) const {
                assert(i > 0);
                size_type j = 0;
                size_type result = 0;
                // Skip blocks until the `i`th zero bit is in block `j`.
                // `result` is the number of bits in the previous blocks.
                // `i` is updated such that it corresponds to the `i`th zero bit in block `j`.
                while (j < num_blocks - 1) {
                    size_type count = block_width - std::popcount(m_bits[j]);
                    if (i > count) {
                        i -= count;
                        result += block_width;
                        j++;
                    } else { break; }
                }
                if (i == 0) return result;
                assert(j < num_blocks);

                // Use select_1(i) on inverted block to calculate select_0(i)
                block_type block = ~m_bits[j];
                return result + ads::select(block, i);
            }

            /**
             * Remove the bit at position `i` \in [0..size-1].
             * @param i
             * @return The bit that has been deleted.
             */
            bool remove(size_type i) {
                assert(i < size());
                assert(size() != 0);
                size_type j = i / block_width;
                size_type k = i % block_width;

                bool deleted_bit = (m_bits[j] >> k) & bit_mask(0);

                auto lo_mask = lower_bit_mask(k);
                auto hi_mask = ~(lower_bit_mask(k + 1));
                assert((lo_mask | hi_mask) == ~(bit_mask(k)));

                // [x_0,...,x_{k-1}][  x_k  ][x_{k+1},...,x_{B-1}]
                // [x_0,...,x_{k-1}][x_{k+1},...,x_{B-1}][       ]
                m_bits[j] = (m_bits[j] & lo_mask) | ((m_bits[j] & hi_mask) >> 1);

                // Shift down the higher blocks by one position and place lowest bit at the highest position of previous block.
                while (++j < num_blocks) {
                    auto m = ((m_bits[j] & static_cast<block_type>(1)) << (block_width - 1));
                    m_bits[j - 1] |= m;
                    m_bits[j] >>= 1;
                }
                m_size--;
                return deleted_bit;
            }

            /**
             * Inserts a bit at position `i` \in [0..size].
             * @param i
             * @param b
             * @return {overflow, high_bit} Whether or not the leaf was previously full and if yes, the high bit that has been pushed out.
             */
            [[nodiscard]] std::pair<bool, bool> insert(size_type i, bool b) {
                assert(i <= size());
                if (i == capacity()) return {true, b};
                assert(i < capacity());
                size_type j = i / block_width;
                size_type k = i % block_width;

                bool highest_bit_is_set = m_bits[j] & bit_mask(block_width - 1);
                auto lo_mask = lower_bit_mask(k);
                auto hi_mask = ~lo_mask;
                assert((lo_mask | hi_mask) == ~(static_cast<block_type>(0)));

                // Leave lower bits and shift up higher bits. x_{B-1} was saved to `highest_bit_is_set`.
                // [x_0,...,x_{k-1}][x_k,   ...  ,x_{B-1}]
                // [x_0,...,x_{k-1}][ b ][x_k,...,x_{B-2}]
                m_bits[j] = (m_bits[j] & lo_mask) | ((m_bits[j] & hi_mask) << 1);
                if (b) m_bits[j] |= bit_mask(k);

                // Shift up the higher blocks by one position and insert the last bit of previous block at the lowest position.
                while (++j < num_blocks) {
                    bool prev_highest_bit_is_set = highest_bit_is_set;
                    highest_bit_is_set = (m_bits[j] >> (block_width - 1)) & bit_mask(0);
                    m_bits[j] <<= 1;
                    if (prev_highest_bit_is_set) m_bits[j] |= bit_mask(0);
                }

                if (size() == capacity()) {
                    return {true, highest_bit_is_set};
                }
                assert(size() < capacity());
                m_size++;
                return {false, false};
            }

            friend std::ostream &operator<<(std::ostream &os, const Leaf &leaf) {
                os << "Leaf(size=" << leaf.size() << ", bits=";
                for (size_type i = 0; i < leaf.size(); ++i) {
                    os << (int) leaf.access(i);
                }
                return os << ")";
            }

            [[nodiscard]] static constexpr block_type bit_mask(size_type i) {
                // Returns block with bit at position `i` set.
                assert(i < block_width);
                return static_cast<block_type>(1) << i;
            }

            [[nodiscard]] static constexpr block_type lower_bit_mask(size_type i) {
                // Returns block with all bits at positions less than `i` set.
                assert(i <= block_width);
                if (i == block_width)
                    return ~static_cast<block_type>(0);
                return (static_cast<block_type>(1) << i) - static_cast<block_type>(1);
            }
        };

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

        /**
         * Split up one full leaf into two leaves of equal size and append `high_bit` to the right one.
         * @param leaf
         * @param high_bit
         * @return the new inner node that is parent to the two leaves
         */
        static Inner *split(Leaf *leaf, bool high_bit) {
            assert(leaf->size() == leaf->capacity());
            static_assert(Leaf::num_blocks > 1);

            Leaf *left = leaf;
            Leaf *right = new Leaf;

            // Copy blocks from left to right, starting at the middle of left and the first block of right.
            size_type k = 0;
            for (size_type j = Leaf::num_blocks / 2; j < Leaf::num_blocks; ++j, ++k) {
                right->m_bits[k] = left->m_bits[j];
                left->m_bits[j] = static_cast<typename Leaf::block_type>(0);
            }
            // Adjust the size of both leaves.
            left->m_size -= k * Leaf::block_width;
            right->m_size += k * Leaf::block_width;
            assert(k < Leaf::num_blocks);

            // Append `high_bit` to the right leaf.
            if (high_bit)
                right->m_bits[k] |= Leaf::bit_mask(0);
            right->m_size++;

            auto left_size = left->size();
            auto left_ones = left->rank1(left_size);

            auto *new_inner = new Inner{
                    NodeHandle::leaf(left), NodeHandle::leaf(right),
                    left_size, left_ones,
                    1};
            assert(std::get<0>(check_integrity(NodeHandle::inner(new_inner))));
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
                    return {NodeHandle::inner(split(leaf, high_bit)), true};
                } else {
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
                        return {NodeHandle::inner(balance(inner)), true};
                    }
                }
                return {node, false};
            }
        }

        void insert_iter_(size_type i, bool b) {
            //std::cout << "insert(" << i << ", " << b << ")\n";
            constexpr static bool left = true;
            constexpr static bool right = false;

            auto node = m_root;
            auto &path = m_path;
            path.clear();

            auto [inner, leaf] = node.cast();
            while (inner) {
                if (i < inner->left_size) {
                    inner->left_size++;
                    inner->left_ones += b;

                    path.emplace_back(inner, left);
                    std::tie(inner, leaf) = inner->left.cast();
                } else {
                    i -= inner->left_size;
                    path.emplace_back(inner, right);
                    std::tie(inner, leaf) = inner->right.cast();
                }
            }

            auto [overflow, high_bit] = leaf->insert(i, b);
            if (overflow) {
                //std::cout << "overflow\n";
                node = NodeHandle::inner(split(leaf, high_bit));
                if (path.empty()) {
                    m_root = node;
                    return;
                }
                bool was_left_child{};

                size_type h = 1;
                while (!path.empty()) {
                    ++h;

                    std::tie(inner, was_left_child) = path.back();
                    path.pop_back();

                    //std::cout << "   " << inner << " " << (was_left_child ? "left" : "right") << " " << inner->height
                    //          << " " << height << "\n";

                    if (inner->height + 1 == h) {
                        inner->height = h;
                        if (was_left_child) {
                            inner->left = node;
                            inner = balance(inner);
                            assert(std::get<0>(check_integrity(NodeHandle::inner(inner))));
                            //assert(inner_prime == inner);
                            assert(std::abs(balance_factor(*inner)) < 2);
                            node = NodeHandle::inner(inner);
                        } else {
                            inner->right = node;
                            inner = balance(inner);
                            //assert(inner_prime == inner);
                            assert(std::abs(balance_factor(*inner)) < 2);
                            assert(std::get<0>(check_integrity(NodeHandle::inner(inner))));
                            node = NodeHandle::inner(inner);
                        }
                    } else {
                        if (was_left_child) { inner->left = node; } else { inner->right = node; }
                        assert(inner->height == std::max(height(inner->left), height(inner->right)) + 1);
                        assert(node.height() < inner->height);
                        assert(std::get<0>(check_integrity(NodeHandle::inner(inner))));
                        assert(std::abs(balance_factor(*inner)) < 2);
                        return;
                    }
                }
                m_root = node;
                assert(std::get<0>(check_integrity(m_root)));
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
            assert(std::get<0>(check_integrity(NodeHandle::inner(a))));
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

                    assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

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
                        [[maybe_unused]] auto [overflow_, high_bit_] = b_leaf->insert(b_leaf->size(), moved_bit);
                        assert(!overflow_);

                        a->left_ones -= deleted_bit;
                        a->left_ones += moved_bit;
                        a->height = std::max(height(a->left), height(a->right)) + 1;

                        assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

                        return {NodeHandle::inner(balance(a)), deleted_bit};
                    } else {
                        // 1.2.2     and the leaf has more than the minimum number of bits.
                        // Strategy: Directly delete bit from the left child.
                        assert(b_leaf->size() > Leaf::min_num_bits);
                        auto deleted_bit = b_leaf->remove(i);

                        a->left_size--;
                        a->left_ones -= deleted_bit;

                        assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

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

                    assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

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
                        [[maybe_unused]] auto [overflow_, high_bit_] = b_leaf->insert(0, moved_bit);
                        assert(!overflow_);

                        a->left_size--;
                        a->left_ones -= moved_bit;
                        a->height = std::max(height(a->left), height(a->right)) + 1;

                        assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

                        return {NodeHandle::inner(balance(a)), deleted_bit};
                    } else {
                        // 2.2.2     and the leaf has more than the minimum number of bits.
                        // Strategy: Directly delete bit from the right child.
                        assert(b_leaf->size() > Leaf::min_num_bits);
                        return {NodeHandle::inner(a), b_leaf->remove(i)};
                    }
                }
            }
            assert(false);
        }

        static std::pair<NodeHandle, bool> remove_first_level(Inner *a, size_type i) {
            assert(std::get<0>(check_integrity(NodeHandle::inner(a))));
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

                if (b->size() == Leaf::min_num_bits) {
                    // 1.1     and the left child only has the minimum number of bits left

                    if (c->size() == Leaf::min_num_bits) {
                        // 1.1.1     and the right child only has the minimum number of bits left.
                        // Strategy: Merge both children and remove the bit in the resulting leaf.

                        auto *leaf = merge(a);
                        auto deleted_bit = leaf->remove(i);
                        return {NodeHandle::leaf(leaf), deleted_bit};
                    } else {
                        // 1.1.2     and the right child has more than the minimum number of bits left.
                        // Strategy: Delete bit from the left child
                        //           and append the first bit of the right child to the left child.

                        assert(c->size() > Leaf::min_num_bits);

                        auto deleted_bit = b->remove(i);
                        auto moved_bit = c->remove(0);
                        [[maybe_unused]] auto [overflow_, high_bit_] = b->insert(b->size(), moved_bit);
                        assert(!overflow_);

                        a->left_ones -= deleted_bit;
                        a->left_ones += moved_bit;

                        return {NodeHandle::inner(a), deleted_bit};
                    }
                } else {
                    // 1.2     and the left child has more than the minimum number of bits left.
                    // Strategy: Directly delete bit from left child.

                    assert(b->size() > Leaf::min_num_bits);

                    auto deleted_bit = b->remove(i);
                    a->left_size--;
                    a->left_ones -= deleted_bit;
                    return {NodeHandle::inner(a), deleted_bit};
                }
            } else {
                // 2.    i is in the right child

                if (c->size() == Leaf::min_num_bits) {
                    // 2.1     and the right child only has the minimum number of bits left

                    if (b->size() == Leaf::min_num_bits) {
                        // 2.1.1     and the left child only has the minimum number of bits left.
                        // Strategy: Merge both children and remove the bit in the resulting leaf.

                        auto *leaf = merge(a);
                        auto deleted_bit = leaf->remove(i);
                        return {NodeHandle::leaf(leaf), deleted_bit};
                    } else {
                        // 2.1.2     and the left child has more than the minimum number of bits left.
                        // Strategy: Delete bit from the right child
                        //           and prepend the last bit of the left child to the right child.

                        assert(b->size() > Leaf::min_num_bits);

                        auto deleted_bit = c->remove(i - a->left_size);
                        auto moved_bit = b->remove(b->size() - 1);
                        [[maybe_unused]] auto [overflow_, high_bit_] = c->insert(0, moved_bit);
                        assert(!overflow_);

                        a->left_size--;
                        a->left_ones -= moved_bit;

                        return {NodeHandle::inner(a), deleted_bit};
                    }
                } else {
                    // 2.2     and the right child has more than the minimum number of bits left.
                    // Strategy: Directly delete bit from right child.

                    assert(c->size() > Leaf::min_num_bits);
                    return {NodeHandle::inner(a), c->remove(i - a->left_size)};
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
            assert(b->size() == Leaf::min_num_bits);
            assert(c->size() == Leaf::min_num_bits);

            if constexpr (Leaf::min_num_bits % Leaf::block_width == 0) {
                // If the minimum number of bits aligns with block boundaries, we can copy whole blocks.
                constexpr auto min_num_blocks = Leaf::min_num_bits / Leaf::block_width;
                size_type k = 0;
                for (size_type j = min_num_blocks; j < 2 * min_num_blocks; ++j, ++k) {
                    b->m_bits[j] = c->m_bits[k];
                }
                b->m_size += k * Leaf::block_width;
            } else {
                for (size_type j = 0; j < Leaf::min_num_bits; ++j) {
                    [[maybe_unused]] auto [overflow_, high_bit_] = b->insert(b->size(), c->access(j));
                    assert(!overflow_);
                }
            }
            assert(b->size() == 2 * Leaf::min_num_bits);
            delete a;
            delete c;
            return b;
        }

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

            assert(std::get<0>(check_integrity(NodeHandle::inner(c))));
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

            assert(std::get<0>(check_integrity(NodeHandle::inner(c))));
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

            assert(std::get<0>(check_integrity(NodeHandle::inner(b))));
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

#ifndef NDEBUG

        static std::tuple<bool, size_type, size_type, size_type> check_integrity(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                auto [left_ok, left_size, left_ones, left_height] = check_integrity(inner->left);
                auto [right_ok, right_size, right_ones, right_height] = check_integrity(inner->right);
                auto height = std::max(left_height, right_height) + 1;
                bool ok = left_ok && right_ok
                          && (inner->left_size == left_size)
                          && (inner->left_ones == left_ones)
                          && (inner->height == height);
                return {ok, left_size + right_size, left_ones + right_ones, height};
            }
            if (leaf) {
                auto size = leaf->size();
                auto ones = leaf->rank1(leaf->size());
                return {true, size, ones, height(node)};
            }
            return {false, 0, 0, 0};
        }

#endif

        static size_type required_bits(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                return required_bits(inner->left) + required_bits(inner->right) + sizeof(Inner) * 8;
            } else {
                return sizeof(Leaf) * 8;
            }
        }

    public:
        DynamicBitVector() : m_root(NodeHandle::leaf(new Leaf)), m_size(0) {}

        explicit DynamicBitVector(const std::vector<bool> &bits) : m_root(NodeHandle::leaf(new Leaf)), m_size(0) {
            for (auto bit: bits) {
                insert(size(), bit);
            }
        }

        ~DynamicBitVector() { delete_node(m_root); }

        [[nodiscard]] bool access(size_type i) const { return access(m_root, i); }

        void insert(size_type i, bool b) {
            assert(i <= size());
            bool changed{};
            std::tie(m_root, changed) = insert(m_root, i, b);
            //insert_iter_(i, b);
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
            return required_bits(m_root) + (sizeof m_root) * 8;
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
        std::vector<std::pair<Inner *, bool>> m_path;
    };
}

#endif //ADVANCED_DATA_STRUCTURES_DYNAMICBITVECTOR_H
