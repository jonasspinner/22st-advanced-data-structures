#ifndef ADVANCED_DATA_STRUCTURES_BITITERATOR_H
#define ADVANCED_DATA_STRUCTURES_BITITERATOR_H

#include <cstdint>
#include <type_traits>
#include <limits>
#include <cassert>

namespace ads {
    template<class BlockType>
    class BlockBitRange {
    private:
        static_assert(std::is_unsigned_v<BlockType>);
        using block_type = BlockType;
        using size_type = uint32_t;
        static constexpr size_type block_width = std::numeric_limits<block_type>::digits;

        class EndIterator {
        };

        class Iterator {
        public:
            constexpr Iterator(block_type block, size_type remaining) : m_block(block), m_remaining(remaining) {}

            [[nodiscard]] constexpr bool operator*() const { return m_block & 1; }

            constexpr Iterator &operator++() {
                m_block >>= 1;
                --m_remaining;
                return *this;
            }

            [[nodiscard]] constexpr bool operator!=(EndIterator) { return m_remaining != 0; }

        private:
            block_type m_block{};
            size_type m_remaining{};
        };

    public:
        constexpr explicit BlockBitRange(block_type block, size_type start = 0, size_type end = block_width) : m_begin(
                block >> start, end - start) {
            assert(0 <= start);
            assert(start <= end);
            assert(end <= block_width);
        }

        [[nodiscard]] constexpr auto begin() const { return m_begin; }

        [[nodiscard]] constexpr auto end() const { return EndIterator{}; }

    private:
        Iterator m_begin;
    };
}

#endif //ADVANCED_DATA_STRUCTURES_BITITERATOR_H
