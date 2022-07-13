#ifndef ADVANCED_DATA_STRUCTURES_SMALLSTATICBITVECTOR_H
#define ADVANCED_DATA_STRUCTURES_SMALLSTATICBITVECTOR_H

#include <cstdint>
#include <array>
#include <limits>
#include <bit>
#include <cassert>
#include <iostream>

#include "utils.h"

namespace ads {
    template<class BlockType, std::size_t MaxNumBlocks>
    class SmallStaticBitVector {
    private:
        using size_type = int32_t;
        using block_type = BlockType;
        constexpr static size_type block_width = std::numeric_limits<block_type>::digits; // Number of bits per block.
        constexpr static size_type max_num_blocks = MaxNumBlocks;
        constexpr static size_type max_num_bits = max_num_blocks * block_width;

        static_assert(max_num_blocks > 1);
        using Blocks = std::array<block_type, max_num_blocks>;

    private:
        Blocks m_blocks{};
        size_type m_size{};
    public:
        SmallStaticBitVector() = default;

        SmallStaticBitVector(Blocks blocks, size_type size) : m_blocks(blocks), m_size(size) {
            assert(size <= max_num_blocks * block_width);
        }

        [[nodiscard]] size_type required_bits_upperbound() const {
            return (sizeof m_blocks + sizeof m_size) * 8;
        }

        /**
         * @return the number of bits stored.
         */
        [[nodiscard]] constexpr size_type size() const { return m_size; }

        [[nodiscard]] constexpr bool empty() const { return size() == 0; }

        [[nodiscard]] constexpr size_type num_blocks() const { return (size() + block_width - 1) / block_width; }

        [[nodiscard]] constexpr block_type blocks(size_type j) const {
            assert(j < num_blocks());
            return m_blocks[j];
        }

        void clear() {
            for (size_type j = 0; j < num_blocks(); ++j) {
                m_blocks[j] = 0;
            }
            m_size = 0;
        }

        [[nodiscard]] bool access(size_type i) const {
            assert(i < size());
            size_type j = i / block_width;
            size_type k = i % block_width;
            return (m_blocks[j] >> k) & bit_mask(0);
        }

        void flip(size_type i) {
            assert(i < size());
            size_type j = i / block_width;
            size_type k = i % block_width;
            m_blocks[j] ^= bit_mask(k);
        }

        [[nodiscard]] size_type rank(size_type i, bool b) const {
            return b ? rank1(i) : i - rank1(i);
        }

        [[nodiscard]] size_type rank1(size_type i) const {
            assert(i <= size());
            size_type j = i / block_width;
            size_type rank = 0;
            // Sum up number of ones in all blocks before block `j`.
            for (size_type block_idx = 0; block_idx < num_blocks() && block_idx < j; ++block_idx) {
                rank += std::popcount(m_blocks[block_idx]);
                i -= block_width;
            }
            // This covers the case that `i` equals `size` and prevents memory access after the last block
            if (i == 0) return rank;
            assert(i < block_width);
            // Only count number of ones in the bits lower than `i`.
            block_type masked_block = m_blocks[j] & lower_bit_mask(i);
            rank += std::popcount(masked_block);
            return rank;
        }

        [[nodiscard]] size_type select(size_type i, bool b) const {
            return b ? select1(i) : select0(i);
        }

        [[nodiscard]] size_type select1(size_type i) const {
            assert(i > 0);
            size_type j = 0;
            size_type result = 0;
            // Skip blocks until the `i`th ones bit is in block `j`.
            // `result` is the number of bits in the previous blocks.
            // `i` is updated such that it corresponds to the `i`th one bit in block `j`.
            while (j < num_blocks() - 1) {
                size_type count = std::popcount(m_blocks[j]);
                if (i > count) {
                    i -= count;
                    result += block_width;
                    j++;
                } else { break; }
            }
            if (i == 0) return result;
            assert(j < num_blocks());
            block_type block = m_blocks[j];
            return result + ads::select(block, i);
        }

        [[nodiscard]] size_type select0(size_type i) const {
            assert(i > 0);
            size_type j = 0;
            size_type result = 0;
            // Skip blocks until the `i`th zero bit is in block `j`.
            // `result` is the number of bits in the previous blocks.
            // `i` is updated such that it corresponds to the `i`th zero bit in block `j`.
            while (j < num_blocks() - 1) {
                size_type count = block_width - std::popcount(m_blocks[j]);
                if (i > count) {
                    i -= count;
                    result += block_width;
                    j++;
                } else { break; }
            }
            if (i == 0) return result;
            assert(j < num_blocks());

            // Use select_1(i) on inverted block to calculate select_0(i)
            block_type block = ~m_blocks[j];
            return result + ads::select(block, i);
        }

        [[nodiscard]] size_type num_ones() const {
            size_type num = 0;
            for (size_type j = 0; j < num_blocks(); ++j) {
                num += std::popcount(m_blocks[j]);
            }
            return num;
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

            bool deleted_bit = (m_blocks[j] >> k) & bit_mask(0);

            auto lo_mask = lower_bit_mask(k);
            auto hi_mask = ~(lower_bit_mask(k + 1));
            assert((lo_mask | hi_mask) == ~(bit_mask(k)));

            // [x_0,...,x_{k-1}][  x_k  ][x_{k+1},...,x_{B-1}]
            // [x_0,...,x_{k-1}][x_{k+1},...,x_{B-1}][       ]
            m_blocks[j] = (m_blocks[j] & lo_mask) | ((m_blocks[j] & hi_mask) >> 1);

            // Shift down the higher blocks by one position and place lowest bit at the highest position of previous block.
            while (++j < num_blocks()) {
                auto m = ((m_blocks[j] & static_cast<block_type>(1)) << (block_width - 1));
                m_blocks[j - 1] |= m;
                m_blocks[j] >>= 1;
            }
            m_size--;
            return deleted_bit;
        }

        /**
         * Remove the highest bit.
         * @return The bit that has been deleted
         */
        bool pop_back() {
            assert(!empty());
            m_size--;
            size_type j = size() / block_width;
            size_type k = size() % block_width;
            auto deleted_bit = (m_blocks[j] >> k) & 1;
            m_blocks[j] &= ~bit_mask(k);
            return deleted_bit;
        }

        block_type pop_back_block_aligned() {
            assert(size() >= block_width);
            assert(size() % block_width == 0);
            size_type last_block_idx = num_blocks() - 1;
            auto block = m_blocks[last_block_idx];
            m_blocks[last_block_idx] = 0;
            m_size -= block_width;
            return block;
        }

        bool pop_front() {
            assert(!empty());
            return remove(0);
        }

        block_type pop_front_block() {
            assert(size() >= block_width);
            block_type first_block = m_blocks[0];
            size_type last_block_idx = num_blocks() - 1;
            for (size_type j = 0; j < last_block_idx; ++j) {
                m_blocks[j] = m_blocks[j + 1];
            }
            m_blocks[last_block_idx] = 0;
            m_size -= block_width;
            return first_block;
        }

        /**
         * Inserts a bit at position `i` \in [0..size].
         * @param i
         * @param b
         * @return {overflow, high_bit} Whether or not the leaf was previously full and if yes, the high bit that has been pushed out.
         */
        [[nodiscard]] std::pair<bool, bool> insert(size_type i, bool b) {
            assert(i <= size());
            if (i == max_num_bits) return {true, b};
            assert(i < max_num_bits);
            size_type j = i / block_width;
            size_type k = i % block_width;

            bool highest_bit_is_set = m_blocks[j] & bit_mask(block_width - 1);
            auto lo_mask = lower_bit_mask(k);
            auto hi_mask = ~lo_mask;
            assert((lo_mask | hi_mask) == ~(static_cast<block_type>(0)));

            // Leave lower bits and shift up higher bits. x_{B-1} was saved to `highest_bit_is_set`.
            // [x_0,...,x_{k-1}][x_k,   ...  ,x_{B-1}]
            // [x_0,...,x_{k-1}][ b ][x_k,...,x_{B-2}]
            m_blocks[j] = (m_blocks[j] & lo_mask) | ((m_blocks[j] & hi_mask) << 1);
            if (b) m_blocks[j] |= bit_mask(k);

            // Shift up the higher blocks by one position and insert the last bit of previous block at the lowest position.
            while (++j < num_blocks()) {
                bool prev_highest_bit_is_set = highest_bit_is_set;
                highest_bit_is_set = (m_blocks[j] >> (block_width - 1)) & bit_mask(0);
                m_blocks[j] <<= 1;
                if (prev_highest_bit_is_set) m_blocks[j] |= bit_mask(0);
            }

            if (size() == max_num_bits) {
                return {true, highest_bit_is_set};
            }
            assert(size() < max_num_bits);
            m_size++;
            return {false, false};
        }

        bool push_back(bool b) {
            if (size() == max_num_bits) return true;
            if (b) {
                size_type j = size() / block_width;
                size_type k = size() % block_width;
                m_blocks[j] |= bit_mask(k);
            }
            m_size++;
            return false;
        }

        void push_back_block_aligned(block_type block) {
            assert(size() % block_width == 0);
            assert(size() <= max_num_bits - block_width);
            size_type j = size() / block_width;
            m_blocks[j] = block;
            m_size += block_width;
        }

        [[nodiscard]] std::pair<bool, bool> push_front(bool b) {
            return insert(0, b);
        }

        void push_front_block(block_type block) {
            assert(size() + block_width <= max_num_bits);
            for (size_type j = num_blocks() - 1; j > 0; --j) {
                m_blocks[j] = m_blocks[j - 1];
            }
            m_blocks[0] = block;
            m_size += block_width;
        }

        void concat_block_aligned(const SmallStaticBitVector &right) {
            assert(size() % block_width == 0);
            assert(right.size() % block_width == 0);
            assert(size() + right.size() <= max_num_bits);
            size_type k = 0;
            for (size_type j = size() / block_width; j < (size() + right.size()) / block_width; ++j, ++k) {
                m_blocks[j] = right.m_blocks[k];
            }
            m_size += k * block_width;
            assert(size() <= max_num_bits);
        }

        void split_into_empty(SmallStaticBitVector &right) {
            assert(size() == max_num_bits);
            assert(right.size() == 0);
            // Copy blocks from left to right, starting at the middle of left and the first block of right.
            size_type k = 0;
            static_assert(max_num_blocks % 2 == 0);
            for (size_type j = num_blocks() / 2; j < num_blocks(); ++j, ++k) {
                right.m_blocks[k] = m_blocks[j];
                m_blocks[j] = 0;
            }
            m_size -= k * block_width;
            right.m_size += k * block_width;
        }


        static std::pair<size_type, size_type> balance_sizes(SmallStaticBitVector &left, SmallStaticBitVector &right) {
            assert(false);
            // TODO
            static_assert(max_num_blocks % 2 == 0);

            size_type num_ones_moved_to_left = 0;
            size_type num_ones_moved_to_right = 0;

            const size_type original_left_num_blocks = left.m_size / block_width;
            const size_type original_right_num_blocks = right.m_size / block_width;

            if (original_left_num_blocks + original_right_num_blocks > 2 * max_num_blocks) return {};

            if (left.m_size + 2 * block_width <= right.m_size) {
                if (left.m_size % block_width == 0) {
                    // ab______ cdefgh__
                    const size_type shift = (original_right_num_blocks - original_left_num_blocks) / 2;
                    assert(original_left_num_blocks + shift <= max_num_blocks);
                    assert(original_right_num_blocks >= shift);

                    size_type idx = 0;
                    for (; idx < shift; ++idx) {
                        num_ones_moved_to_left += std::popcount(right.m_blocks[idx]);
                        left.m_blocks[original_left_num_blocks + idx] = right.m_blocks[idx];
                        right.m_blocks[idx] = right.m_blocks[idx + shift];
                    }
                    // abcd____ efefgh__
                    for (; idx < original_right_num_blocks - shift; ++idx)
                        right.m_blocks[idx] = right.m_blocks[idx + shift];
                    // abcd____ efghgh__
                    for (; idx < original_right_num_blocks; ++idx)
                        right.m_blocks[idx] = 0;
                    // abcd____ efgh____

                    left.m_size += shift * block_width;
                    right.m_size -= shift * block_width;
                } else {

                }
            } else if (left.m_size >= right.m_size + 2 * block_width) {
                if (left.m_size % block_width == 0) {
                    // abcdef__ gh______
                    const size_type shift = (original_left_num_blocks - original_right_num_blocks) / 2;
                    assert(original_left_num_blocks >= shift);
                    assert(original_right_num_blocks + shift <= max_num_blocks);

                    size_type idx = original_right_num_blocks - 1;
                    for (; idx >= shift; --idx) {
                        right.m_blocks[idx + shift] = right.m_blocks[idx];
                    }
                    // abcdef__ ghgh____
                    assert(idx = shift - 1);
                    for (; idx >= 0; --idx) {
                        const auto l = original_left_num_blocks - shift + idx;
                        num_ones_moved_to_right += std::popcount(left.m_blocks[l]);
                        right.m_blocks[idx] = left.m_blocks[l];
                        left.m_blocks[l] = 0;
                    }
                    // abcd____ efgh____

                    left.m_size -= shift * block_width;
                    right.m_size += shift * block_width;
                } else {

                }
            }

            return {num_ones_moved_to_left, num_ones_moved_to_right};
        }

        friend std::ostream &operator<<(std::ostream &os, const SmallStaticBitVector &bv) {
            os << "StaticBitVector(size=" << bv.size() << ", bits=";
            for (size_type i = 0; i < bv.size(); ++i) {
                os << (int) bv.access(i);
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
}


#endif //ADVANCED_DATA_STRUCTURES_SMALLSTATICBITVECTOR_H
