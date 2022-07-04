#ifndef ADVANCED_DATA_STRUCTURES_TAGGEDPOINTER_H
#define ADVANCED_DATA_STRUCTURES_TAGGEDPOINTER_H

#include <cstdint>
#include <utility>
#include <cassert>

namespace ads {
    template<class A, class B>
    class TaggedPointer {
        static_assert(sizeof(uintptr_t) == sizeof(void *));
        static_assert(!std::is_same_v<A, B>);
    private:
        constexpr explicit TaggedPointer(void *ptr) : m_ptr(ptr) {
#ifdef NDEBUG
            static_assert(sizeof(TaggedPointer) == sizeof(void *));
#endif
            static_assert(alignof(A) > 1);
            static_assert(alignof(B) > 1);
            assert(ptr);
        }

    public:
        [[nodiscard]] constexpr static TaggedPointer first(A *ptr) {
            assert(ptr);
            TaggedPointer p(ptr);
#ifndef NDEBUG
            p.m_a = ptr;
#endif
            return p;
        }

        [[nodiscard]] constexpr static TaggedPointer second(B *ptr) {
            assert(ptr);
            TaggedPointer p(as_tagged(ptr));
#ifndef NDEBUG
            p.m_b = ptr;
#endif
            return p;
        }

        constexpr explicit operator bool() const {
            return as_untagged(m_ptr) != nullptr;
        }

        [[nodiscard]] constexpr std::pair<A *, B *> cast() {
            if (is_tagged(m_ptr)) {
                return {nullptr, static_cast<B *>(as_untagged(m_ptr))};
            } else {
                return {static_cast<A *>(m_ptr), nullptr};
            }
        }

        [[nodiscard]] constexpr std::pair<const A *, const B *> cast() const {
            if (is_tagged(m_ptr)) {
                return {nullptr, static_cast<const B *>(as_untagged(m_ptr))};
            } else {
                return {static_cast<const A *>(m_ptr), nullptr};
            }
        }

    private:
        [[nodiscard]] constexpr static bool is_tagged(void *ptr) {
            return reinterpret_cast<uintptr_t>(ptr) & mask;
        }

        [[nodiscard]] constexpr static void *as_untagged(void *ptr) {
            return reinterpret_cast<void *>(reinterpret_cast<uintptr_t>(ptr) & ~mask);
        }

        [[nodiscard]] constexpr static void *as_tagged(void *ptr) {
            return reinterpret_cast<void *>(reinterpret_cast<uintptr_t>(ptr) | mask);
        }

        static constexpr uintptr_t mask = 1;

        void *m_ptr{nullptr};
#ifndef NDEBUG
        A* m_a{nullptr};
        B* m_b{nullptr};
#endif
    };
}

#endif //ADVANCED_DATA_STRUCTURES_TAGGEDPOINTER_H
