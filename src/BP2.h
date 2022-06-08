#ifndef ADVANCED_DATA_STRUCTURES_BP2_H
#define ADVANCED_DATA_STRUCTURES_BP2_H

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
    template<class BlockType = uint64_t, size_t NumBlocks = 64>
    class BP2 {
        static_assert(std::is_unsigned_v<BlockType>);
    public:
        using size_type = int32_t;
        using height_type = int32_t;
        constexpr static bool open = false;
        constexpr static bool close = true;
        constexpr static size_type right_end = std::numeric_limits<size_type>::max();
        constexpr static size_type left_end = std::numeric_limits<size_type>::min();
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
            size_type total_excess{};
            size_type min_excess{};
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
            size_type total_excess{};
            size_type min_excess{};

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

                total_excess += (m_bits[j] >> k) & bit_mask(0) ? +1 : -1;
                min_excess = calculate_min_excess();
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

                total_excess += deleted_bit ? 1 : -1;
                min_excess = calculate_min_excess();
                return deleted_bit;
            }

            /**
             * Inserts a bit at position `i` \in [0..size].
             * @param i
             * @param b
             * @return {overflow, high_bit} Whether or not the leaf was previously full and if yes, the high bit that has been pushed out.
             */
            [[nodiscard]] std::pair<bool, bool> insert(size_type i, bool b) {
                // TODO: update total and min excess
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
                    std::tie(total_excess, min_excess) = calculate_excess_info();
                    return {true, highest_bit_is_set};
                }
                assert(size() < capacity());
                m_size++;
                std::tie(total_excess, min_excess) = calculate_excess_info();
                return {false, false};
            }

            [[nodiscard]] std::pair<size_type, size_type> forward_search(size_type i, size_type d) const {
                // NOTE: different definition than slides
                //       forward_search(i, d) = min {j > i : excess(j+1) - excess(i) = d }
                assert(-1 <= i && i < size());
                if (i == -1 && d == 0) return {i, 0};

                auto excess_i = i >= 0 ? i - 2 * rank1(i) : 0;
                assert(i < 0 || excess_i == i - 2 * rank1(i));

                size_type j = i;
                auto excess_jp1 = excess_i + (i >= 0 ? (access(j) ? -1 : 1) : 0);
                assert(excess_jp1 == (j + 1) - 2 * rank1(j + 1));

                ++j;
                // TODO: make more efficient
                for (; j < size(); ++j) {
                    excess_jp1 += access(j) ? -1 : 1;
                    assert(excess_jp1 == (j + 1) - 2 * rank1(j + 1));
                    if (excess_jp1 - excess_i == d) {
                        return {j, 0};
                    }
                }
                return {j, -(excess_jp1 - excess_i - d)};
            }

            [[nodiscard]] std::pair<size_type, size_type> backward_search(size_type i, size_type d) const {
                // NOTE: different definition than slides
                //       forward_search(i, d) = min {j < i : excess(i+1) - excess(j) = d }
                assert(0 <= i && i <= size());
                if (i == size() && d == 0) return {i, 0};

                auto excess_ip1 = i < size() ? (i + 1) - 2 * rank1(i + 1) : (i) - 2 * rank1(i);

                size_type j = i;
                auto excess_j = excess_ip1 - (i < size() ? (access(j) ? -1 : 1) : 0);
                assert(excess_j == j - 2 * rank1(j));

                --j;
                // TODO: make more efficient
                for (; j >= 0; --j) {
                    excess_j -= access(j) ? -1 : 1;
                    assert(excess_j == j - 2 * rank1(j));
                    if (excess_ip1 - excess_j == d) {
                        return {j, 0};
                    }
                }
                return {j, -(excess_ip1 - excess_j - d)};
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

            [[nodiscard]] size_type calculate_total_excess() const {
                auto rank_1 = rank1(size());
                auto rank_0 = size() - rank_1;
                return rank_0 - rank_1;
            }

            [[nodiscard]] size_type calculate_min_excess() const {
                // TODO: make more efficient
                size_type excess = 0;
                size_type current_min_excess = 0;
                for (size_type i = 0; i < size(); ++i) {
                    if (access(i)) {
                        excess--;
                        current_min_excess = std::min(current_min_excess, excess);
                    } else {
                        excess++;
                    }
                }
                return current_min_excess;
            }

            [[nodiscard]] std::pair<size_type, size_type> calculate_excess_info() const {
                size_type excess = 0;
                size_type current_min_excess = 0;

                auto update_for_block = [&](block_type block, size_type num_bits = block_width) {
                    for (size_type k = 0; k < num_bits; ++k, block >>= 1) {
                        if (block & bit_mask(0)) {
                            excess--;
                            current_min_excess = std::min(current_min_excess, excess);
                        } else {
                            excess++;
                        }
                    }
                };

                size_type last_block_idx = (size() - 1) / block_width;

                // first block
                update_for_block(m_bits[0], std::min(size(), block_width));

                for (size_type j = 1; j < last_block_idx; ++j) {
                    size_type block_ones = std::popcount(m_bits[j]);
                    // check if the worst case (block_ones 1s followed by (block_width - block_ones) zeros) would not
                    // result in an update to the min excess.
                    if (excess - block_ones >= current_min_excess) {
                        excess += block_width - 2 * block_ones;
                    } else {
                        update_for_block(m_bits[j]);
                    }
                }
                if (last_block_idx > 0) {
                    // last block
                    size_type block_num_bits = size() - last_block_idx * block_width;
                    assert(1 <= block_num_bits && block_num_bits <= block_width);
                    size_type block_ones = std::popcount(m_bits[last_block_idx]);
                    if (excess - block_ones >= current_min_excess) {
                        excess += block_num_bits - 2 * block_ones;
                    } else {
                        update_for_block(m_bits[last_block_idx], block_num_bits);
                    }
                }


                assert(excess == size() - 2 * rank1(size()));
                assert(current_min_excess == calculate_min_excess());
                return {excess, current_min_excess};
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
            // TODO: update total and min excess
            assert(false);
            // Dispatch dependent on node type.
            auto [inner, leaf] = node.cast();
            inner ? flip(*inner, i) : leaf->flip(i);
        }

        static void flip(Inner &inner, size_type i) {
            // TODO: update total and min excess
            assert(false);
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

        static std::pair<size_type, size_type> forward_search(const NodeHandle node, size_type i, size_type d) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                auto left_size = inner->left_size;

                if (i < left_size) {
                    // descend left
                    auto [idx, d_prime] = forward_search(inner->left, i, d);
                    if (idx == right_end) {
                        if (d_prime >= min_excess(inner->right)) {
                            // descend right
                            auto [idx_r, d_prime_r] = forward_search(inner->right, left_end, d_prime);
                            assert(idx_r != right_end);
                            idx_r += left_size;
                            return {idx_r, d_prime_r};
                        }
                        // ascend
                        d_prime -= total_excess(inner->right);
                        return {idx, d_prime};
                    }
                    return {idx, d_prime};
                } else {
                    // descend right
                    auto [idx, d_prime] = forward_search(inner->right, i - left_size, d);
                    if (idx == right_end) {
                        return {idx, d_prime};
                    }
                    idx += left_size;
                    return {idx, d_prime};
                }
            } else {
                if (i == left_end) { i = -1; }
                auto [idx, d_prime] = leaf->forward_search(i, d);
                if (idx == leaf->m_size) { idx = right_end; }
                return {idx, d_prime};
            }
        }

        static std::pair<size_type, size_type> backward_search(const NodeHandle node, size_type i, size_type d) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                auto left_size = inner->left_size;

                if (i < left_size) {
                    return backward_search(inner->left, i, d);
                } else if (i == right_end){
                    if (total_excess(inner->right) - d >= min_excess(inner->right)) { // TODO: check
                        auto [idx, d_prime] = backward_search(inner->right, right_end, d);
                        assert(idx != left_end);
                        idx += left_size;
                        return {idx, d_prime};
                    }
                    d -= total_excess(inner->left);
                    return backward_search(inner->left, right_end, d);
                } else {
                    auto [idx, d_prime] = backward_search(inner->right, i - left_size, d);
                    if (idx == left_end) {
                        if (total_excess(inner->left) - d_prime >= min_excess(inner->left)) { // TODO: check
                            auto [idx_l, d_prime_l] = backward_search(inner->left, right_end, d_prime);
                            assert(idx_l != left_end);
                            return {idx_l, d_prime_l};
                        }
                        d_prime -= total_excess(inner->left);
                        return {idx, d_prime};
                    }
                    idx += left_size;
                    return {idx, d_prime};
                }
            } else {
                if (i == right_end) { i = leaf->size(); }
                auto [idx, d_prime] = leaf->backward_search(i, d);
                if (idx == -1) { idx = left_end; }
                return {idx, d_prime};
            }
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
            left->total_excess = left_size - 2 * left_ones;
            left->min_excess = left->calculate_min_excess();

            auto right_size = right->size();
            auto right_ones = right->rank1(right_size);
            right->total_excess = right_size - 2 * right_ones;
            right->min_excess = right->calculate_min_excess();

            auto total_excess = left->total_excess + right->total_excess;
            auto min_excess = std::min(left->min_excess, left->total_excess + right->min_excess);

            auto *new_inner = new Inner{
                    NodeHandle::leaf(left), NodeHandle::leaf(right),
                    left_size, left_ones,
                    total_excess, min_excess,
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

                update_excess_info_from_children_info(*inner);

                // If the height of the subtree changed, recalculate height and if it changed, balance the tree if needed.
                if (changed) {
                    auto h = std::max<int>(inner->height, height(inner->left) + 1);
                    if (inner->height != h) {
                        inner->height = h;
                        node = NodeHandle::inner(balance(inner));
                        assert(std::get<0>(check_integrity(node)));
                        return {node, true};
                    }
                }
                return {node, false};
            } else {
                i -= inner->left_size;

                bool changed{};
                std::tie(inner->right, changed) = insert(inner->right, i, b);

                update_excess_info_from_children_info(*inner);

                // If the height of the subtree changed, recalculate height and if it changed, balance the tree if needed.
                if (changed) {
                    auto h = std::max(inner->height, height(inner->right) + 1);
                    if (inner->height != h) {
                        inner->height = h;
                        assert(std::get<0>(check_integrity(NodeHandle::inner(inner))));
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

            update_excess_info_from_children_info(*inner);

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

                    update_excess_info_from_children_info(*a);

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

                        update_excess_info_from_children_info(*a);

                        assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

                        return {NodeHandle::inner(balance(a)), deleted_bit};
                    } else {
                        // 1.2.2     and the leaf has more than the minimum number of bits.
                        // Strategy: Directly delete bit from the left child.
                        assert(b_leaf->size() > Leaf::min_num_bits);
                        auto deleted_bit = b_leaf->remove(i);

                        a->left_size--;
                        a->left_ones -= deleted_bit;

                        update_excess_info_from_children_info(*a);

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

                    update_excess_info_from_children_info(*a);

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

                        update_excess_info_from_children_info(*a);

                        assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

                        return {NodeHandle::inner(balance(a)), deleted_bit};
                    } else {
                        // 2.2.2     and the leaf has more than the minimum number of bits.
                        // Strategy: Directly delete bit from the right child.
                        assert(b_leaf->size() > Leaf::min_num_bits);
                        auto deleted_bit = b_leaf->remove(i);

                        update_excess_info_from_children_info(*a);

                        assert(std::get<0>(check_integrity(NodeHandle::inner(a))));

                        return {NodeHandle::inner(a), deleted_bit};
                    }
                }
            }
            assert(false);
        }

        static std::pair<NodeHandle, bool> remove_first_level(Inner *a, size_type i) {
            // TODO: update total and min excess
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

                        a->total_excess = b->total_excess + c->total_excess;
                        a->min_excess = std::min(b->min_excess, b->total_excess + c->min_excess);

                        return {NodeHandle::inner(a), deleted_bit};
                    }
                } else {
                    // 1.2     and the left child has more than the minimum number of bits left.
                    // Strategy: Directly delete bit from left child.

                    assert(b->size() > Leaf::min_num_bits);

                    auto deleted_bit = b->remove(i);

                    a->left_size--;
                    a->left_ones -= deleted_bit;

                    a->total_excess = b->total_excess + c->total_excess;
                    a->min_excess = std::min(b->min_excess, b->total_excess + c->min_excess);

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
                        c->insert(0, moved_bit);

                        a->left_size--;
                        a->left_ones -= moved_bit;

                        a->total_excess = b->total_excess + c->total_excess;
                        a->min_excess = std::min(b->min_excess, b->total_excess + c->min_excess);

                        return {NodeHandle::inner(a), deleted_bit};
                    }
                } else {
                    // 2.2     and the right child has more than the minimum number of bits left.
                    // Strategy: Directly delete bit from right child.

                    assert(c->size() > Leaf::min_num_bits);
                    auto deleted_bit = c->remove(i - a->left_size);

                    a->total_excess = b->total_excess + c->total_excess;
                    a->min_excess = std::min(b->min_excess, b->total_excess + c->min_excess);

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

            std::tie(b->total_excess, b->min_excess) = b->calculate_excess_info();

            return b;
        }

        static int height(const NodeHandle &node) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                assert(inner->height == std::max(height(inner->left), height(inner->right)) + 1);
                return inner->height;
            } else { return 0; }
        }

        static size_type total_excess(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? inner->total_excess : leaf->total_excess;
        }

        static size_type min_excess(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? inner->min_excess : leaf->min_excess;
        }

        static void update_excess_info_from_children_info(Inner &node) {
            auto left_total_excess = total_excess(node.left);
            node.total_excess = left_total_excess + total_excess(node.right);
            node.min_excess = std::min(min_excess(node.left), left_total_excess + min_excess(node.right));
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

        static size_type total_excess_explicit(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? total_excess_explicit(inner->left) + total_excess_explicit(inner->right)
                         : leaf->calculate_total_excess();
        }

        static size_type min_excess_explicit(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            return inner ? std::min(min_excess_explicit(inner->left),
                                    total_excess(inner->left) + min_excess_explicit(inner->right))
                         : leaf->calculate_min_excess();
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

            b->total_excess = a->total_excess;
            b->min_excess = a->min_excess;
            update_excess_info_from_children_info(*a);

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

            c->total_excess = a->total_excess;
            c->min_excess = a->min_excess;
            update_excess_info_from_children_info(*a);
            update_excess_info_from_children_info(*b);

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

            c->total_excess = a->total_excess;
            c->min_excess = a->min_excess;
            update_excess_info_from_children_info(*a);
            update_excess_info_from_children_info(*b);

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

            b->total_excess = a->total_excess;
            b->min_excess = a->min_excess;
            update_excess_info_from_children_info(*a);

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

        static std::tuple<bool, size_type, size_type, size_type, size_type, size_type> check_integrity(NodeHandle node) {
            auto [inner, leaf] = node.cast();
            if (inner) {
                auto [left_ok, left_size, left_ones, left_height, left_total_excess, left_min_excess] = check_integrity(
                        inner->left);
                auto [right_ok, right_size, right_ones, right_height, right_total_excess, right_min_excess] = check_integrity(
                        inner->right);
                auto height = std::max(left_height, right_height) + 1;
                auto total_excess = left_total_excess + right_total_excess;
                auto min_excess = std::min(left_min_excess, left_total_excess + right_min_excess);
                bool ok = left_ok && right_ok
                          && (inner->left_size == left_size)
                          && (inner->left_ones == left_ones)
                          && (inner->height == height)
                          && (inner->total_excess == total_excess)
                          && (inner->min_excess == min_excess);
                return {ok, left_size + right_size, left_ones + right_ones, height, total_excess, min_excess};
            }
            if (leaf) {
                auto size = leaf->size();
                auto ones = leaf->rank1(leaf->size());
                auto total_excess = leaf->calculate_total_excess();
                auto min_excess = leaf->calculate_min_excess();
                bool ok = (leaf->total_excess == total_excess)
                          && (leaf->min_excess == min_excess);
                return {ok, size, ones, height(node), total_excess, min_excess};
            }
            return {false, 0, 0, 0, 0, 0};
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

        [[nodiscard]] bool access(size_type i) const { return access(m_root, i); }

        void insert(size_type i, bool b) {
            assert(i <= size());
            bool changed{};
            std::tie(m_root, changed) = insert(m_root, i, b);
            m_bv_size++;

#ifndef NDEBUG
            auto [ok, size, ones, height, total_excess, min_excess] = check_integrity(m_root);
            assert(ok && (m_bv_size == size)); // && (total_excess == 0) && (min_excess >= 0));
#endif
        }

        void remove(size_type i) {
            assert(i < size());
            bool deleted_bit{};
            std::tie(m_root, deleted_bit) = remove(m_root, i);
            m_bv_size--;

#ifndef NDEBUG
            auto [ok, size, ones, height, total_excess, min_excess] = check_integrity(m_root);
            assert(ok && (m_bv_size == size) && (total_excess == 0) && (min_excess >= 0));
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
            return m_bv_size;
        }

        [[nodiscard]] size_type find_close(size_type i) const {
            assert(access(i) == open);
            return forward_search(m_root, i, 0);
        }

        [[nodiscard]] size_type find_open(size_type i) const {
            assert(access(i) == close);
            return backward_search(m_root, i, 0);
        }

        [[nodiscard]] size_type excess(size_type i) const {
            return rank(i, open) - rank(i, close);
        }

        [[nodiscard]] size_type enclose(size_type i) const {
            assert(access(i) == open);
            return backward_search(m_root, i, 2);
        }


    public:
        BP2() : m_root(NodeHandle::leaf(new Leaf)), m_bv_size(0) {
            insert(size(), open);
            insert(size(), close);
        }

        explicit BP2(const std::vector<bool> &bits) : m_root(NodeHandle::leaf(new Leaf)), m_bv_size(0) {
            for (auto bit: bits) {
                insert(size(), bit);
            }
        }

        ~BP2() { delete_node(m_root); }

        void insertchild(size_type v, size_type i, size_type k) {
            // TODO
            auto v_idx = select(v+1, open);
            assert(access(v_idx) == open);

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
            return 0;
        }

        size_type parent(size_type v) {
            return enclose(v);
        }

        size_type subtree_size(size_type v) {
            return (find_close(v) - v + 1) / 2;
        }


        friend std::ostream &operator<<(std::ostream &os, const BP2 &bv) {
            print(os, bv.m_root, 0);
            return os;
        }

        [[nodiscard]] size_type required_bits_upperbound() const {
            return required_bits(m_root) + (sizeof m_root) * 8;
        }

        [[nodiscard]] size_type num_nodes() const {
            return m_bv_size / 2;
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

        FRIEND_TEST(BP2Test, Leaf);

        FRIEND_TEST(BP2Test, Inner);

        FRIEND_TEST(BP2Test, Insert);

        FRIEND_TEST(BP2Test, Remove);

        FRIEND_TEST(BP2Test, BVInsert);

        NodeHandle m_root;
        size_type m_bv_size;
    };
}

#endif //ADVANCED_DATA_STRUCTURES_BP2_H
