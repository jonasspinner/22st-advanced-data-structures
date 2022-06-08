#ifndef ADVANCED_DATA_STRUCTURES_STATICBITVECTOR_H
#define ADVANCED_DATA_STRUCTURES_STATICBITVECTOR_H

#include <cstdint>
#include <array>
#include <limits>
#include <bit>
#include <cassert>
#include <iostream>

#include "utils.h"

namespace ads {
    template<class BlockType, std::size_t NumBlocks>
    class StaticBitVector {
    private:
        using size_type = int32_t;
        using block_type = BlockType;
        constexpr static size_type num_blocks = NumBlocks;

        static_assert(num_blocks > 1);
        using Bits = std::array<block_type, num_blocks>;
        constexpr static size_type block_width = std::numeric_limits<block_type>::digits; // Number of bits per block.

    private:
        Bits m_bits{};
        size_type m_size{};
    public:

        /**
         * @return the number of bits stored.
         */
        [[nodiscard]] constexpr size_type size() const { return m_size; }

        /**
         * @return the maximum number of bits that can be stored.
         */
        [[nodiscard]] constexpr size_type capacity() const { return m_bits.size() * block_width; }

        [[nodiscard]] constexpr block_type blocks(size_type j) const {
            assert(j < m_bits.size());
            return m_bits[j];
        }

        void clear() {
            for (size_type j = 0; j < (size() - 1) / block_width; ++j) {
                m_bits[j] = 0;
            }
            m_size = 0;
        }

        [[nodiscard]] bool access(size_type i) const {
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

        [[nodiscard]] size_type rank(size_type i, bool b) const {
            return b ? rank1(i) : i - rank1(i);
        }

        [[nodiscard]] size_type rank1(size_type i) const {
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

        [[nodiscard]] size_type num_ones() const {
            size_type last_block_idx = (size() - 1) / block_width;
            size_type num = 0;
            for (size_type j = 0; j <= last_block_idx; ++j) {
                num += std::popcount(m_bits[j]);
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

        block_type remove_first_block() {
            assert(size() >= block_width);
            block_type first_block = m_bits[0];
            size_type last_block_idx = (size() - 1) / block_width;
            for (size_type j = 0; j < last_block_idx; ++j) {
                m_bits[j] = m_bits[j+1];
            }
            m_bits[last_block_idx] = 0;
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

        bool append(bool b) {
            if (size() == capacity()) return true;
            size_type j = size() / block_width;
            size_type k = size() % block_width;
            m_bits[j] |= bit_mask(k);
            m_size++;
            return false;
        }

        void append_block_aligned(block_type block) {
            assert(size() % block_width == 0);
            assert(size() <= capacity() - block_width);
            size_type j = size() / block_width;
            m_bits[j] = block;
            m_size += block_width;
        }

        friend std::ostream &operator<<(std::ostream &os, const StaticBitVector &bv) {
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


#endif //ADVANCED_DATA_STRUCTURES_STATICBITVECTOR_H
