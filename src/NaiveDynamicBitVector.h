#ifndef ADVANCED_DATA_STRUCTURES_NAIVEDYNAMICBITVECTOR_H
#define ADVANCED_DATA_STRUCTURES_NAIVEDYNAMICBITVECTOR_H

#include <cstddef>
#include <cassert>
#include <iostream>
#include <vector>

namespace ads {
    class NaiveDynamicBitVector {
    public:
        using size_type = std::size_t;

        NaiveDynamicBitVector() = default;

        explicit NaiveDynamicBitVector(size_type n) : m_bits(n) {}

        explicit NaiveDynamicBitVector(const std::vector<bool> &bits) : m_bits(bits) {}

        explicit NaiveDynamicBitVector(std::vector<bool> &&bits) : m_bits(std::move(bits)) {}

        [[nodiscard]] bool access(size_type i) const {
            return m_bits[i];
        }

        void insert(size_type i, bool b) {
            assert(i <= size());
            m_bits.insert(m_bits.begin() + i, b);
        }

        void remove(size_type i) {
            assert(i < size());
            m_bits.erase(m_bits.begin() + i);
        }

        void flip(size_type i) {
            assert(i < size());
            m_bits[i] = !m_bits[i];
        }

        [[nodiscard]] size_type rank(size_type i, bool b) const {
            assert(i <= size());
            size_type rank = 0;
            for (size_type j = 0; j < i; ++j) {
                rank += m_bits[j] == b;
            }
            return rank;
        }

        [[nodiscard]] size_type select(size_type i, bool b) const {
            assert(i > 0);
            for (size_type j = 0; j < size(); ++j)
                if (m_bits[j] == b && --i == 0) return j;
            throw std::invalid_argument("there is no ith bit of the given type");
        }

        [[nodiscard]] size_type size() const { return m_bits.size(); }

        [[nodiscard]] bool empty() const { return m_bits.empty(); }

    private:
        friend std::ostream &operator<<(std::ostream &os, const NaiveDynamicBitVector &bv) {
            for (auto b: bv.m_bits) { os << b; }
            return os;
        }

        std::vector<bool> m_bits;
    };
}

#endif //ADVANCED_DATA_STRUCTURES_NAIVEDYNAMICBITVECTOR_H
