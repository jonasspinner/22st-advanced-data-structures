#ifndef ADVANCED_DATA_STRUCTURES_AVLTREE_H
#define ADVANCED_DATA_STRUCTURES_AVLTREE_H

#include <algorithm>

#include "TaggedPointer.h"

namespace ads {
    template<class InnerNodeData, class LeafNodeData, class InnerNodeInfoUpdater>
    class avl_tree_options {
    public:
        using inner_node_data_type = InnerNodeData;
        using leaf_node_data_type = LeafNodeData;
        using inner_node_data_updater = InnerNodeInfoUpdater;

    private:
        static_assert(std::is_invocable_r_v<void, inner_node_data_updater,
                inner_node_data_type &, const inner_node_data_type &, const inner_node_data_type &>);
        static_assert(std::is_invocable_r_v<void, inner_node_data_updater,
                inner_node_data_type &, const inner_node_data_type &, const leaf_node_data_type &>);
        static_assert(std::is_invocable_r_v<void, inner_node_data_updater,
                inner_node_data_type &, const leaf_node_data_type &, const inner_node_data_type &>);
        static_assert(std::is_invocable_r_v<void, inner_node_data_updater,
                inner_node_data_type &, const leaf_node_data_type &, const leaf_node_data_type &>);
    };

    enum class direction : uint8_t {
        left, right, stop
    };
    template <class InnerNodeData, class LeafNodeData, class Navigator>
    class avl_tree_query_options {
    public:
        using inner_node_data_type = InnerNodeData;
        using leaf_node_data_type = LeafNodeData;
        using navigator = Navigator;
    private:
        static_assert(std::is_invocable_r_v<direction, navigator, const inner_node_data_type &, const inner_node_data_type &>);
        static_assert(std::is_invocable_r_v<direction, navigator, const inner_node_data_type &, const leaf_node_data_type &>);
        static_assert(std::is_invocable_r_v<direction, navigator, const leaf_node_data_type &, const inner_node_data_type &>);
        static_assert(std::is_invocable_r_v<direction, navigator, const leaf_node_data_type &, const leaf_node_data_type &>);
    };

    template<class Options>
    class avl_tree {
    private:
        using height_type = int32_t;
        using balance_type = height_type;

        using inner_node_data_type = typename Options::inner_node_data_type;
        using leaf_node_data_type = typename Options::leaf_node_data_type;
        using inner_node_data_updater = typename Options::inner_node_data_updater;

        struct Inner;
        struct Leaf;
        using node_ptr = TaggedPointer<Inner, Leaf>;
        struct Inner {
            inner_node_data_type data{};
            node_ptr left, right;
            height_type height{1};
        };
        struct Leaf {
            leaf_node_data_type data{};
        };

        static height_type height(node_ptr node) {
            auto [inner, leaf] = node.cast();
            return inner ? inner->height : 0;
        }

        static balance_type balance_factor(const Inner &inner) {
            return height(inner.left) - height(inner.right);
        }

        static balance_type balance_factor(node_ptr node) {
            auto [inner, leaf] = node.cast();
            return inner ? balance_factor(*inner) : 0;
        }

        static Inner *balance(Inner *node);

        static Inner *rotate_left_left(Inner *a);

        static Inner *rotate_left_right(Inner *a);

        static Inner *rotate_right_left(Inner *a);

        static Inner *rotate_right_right(Inner *a);

        static void update(Inner &node) {
            auto [left_inner, left_leaf] = node.left.cast();
            auto [right_inner, right_leaf] = node.right.cast();
            if (left_inner) {
                if (right_inner) {
                    inner_node_data_updater{}(node.data, left_inner->data, right_inner->data);
                } else {
                    inner_node_data_updater{}(node.data, left_inner->data, right_leaf->data);
                }
            } else {
                if (right_inner) {
                    inner_node_data_updater{}(node.data, left_leaf->data, right_inner->data);
                } else {
                    inner_node_data_updater{}(node.data, left_leaf->data, right_leaf->data);
                }
            }
        }

    public:
        avl_tree() : m_root(node_ptr::second(new Leaf)) {};

        template<class Navigator>
        void query(Navigator &navigator) const {
            auto [inner, leaf] = m_root.cast();
            return inner ? query(navigator, *inner) : query(navigator, *leaf);
        }


        template<class Navigator>
        static void query(Navigator &navigator, const Inner &node) {
            auto [left_inner, left_leaf] = node.left.cast();
            auto [right_inner, right_leaf] = node.right.cast();

            if (left_inner) {
                if (right_inner) {
                    direction dir = navigator.visit_inner(left_inner->data, right_inner->data);
                    switch (dir) {
                        case direction::left: return query(navigator, *left_inner);
                        case direction::right: return query(navigator, *right_inner);
                        case direction::stop: return;
                    }
                } else {
                    direction dir = navigator.visit_inner(left_inner->data, right_leaf->data);
                    switch (dir) {
                        case direction::left: return query(navigator, *left_inner);
                        case direction::right: return query(navigator, *right_leaf);
                        case direction::stop: return;
                    }
                }
            } else {
                if (right_inner) {
                    direction dir = navigator.visit_inner(left_leaf->data, right_inner->data);
                    switch (dir) {
                        case direction::left: return query(navigator, *left_leaf);
                        case direction::right: return query(navigator, *right_inner);
                        case direction::stop: return;
                    }
                } else {
                    direction dir = navigator.visit_inner(left_leaf->data, right_leaf->data);
                    switch (dir) {
                        case direction::left: return query(navigator, *left_leaf);
                        case direction::right: return query(navigator, *right_leaf);
                        case direction::stop: return;
                    }
                }
            }
        }

        template<class Navigator>
        static void query(Navigator &navigator, const Leaf &node) {
            navigator.visit_leaf(node.data);
        }

    private:
        node_ptr m_root;
    };


    template<class Options>
    auto avl_tree<Options>::balance(Inner *node) -> Inner* {
        static_assert(std::is_signed_v<balance_type>);
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

    template<class Options>
    auto avl_tree<Options>::rotate_left_left(Inner *a) -> Inner * {
        //     a            b
        //   b       ->   c   a
        // c
        assert(a);
        Inner *b = a->left.cast().first;
        assert(b);

        a->left = b->right;
        b->right = node_ptr::first(a);

        update(*a);
        update(*b);

        a->height = std::max(height(a->left), height(a->right)) + 1;
        b->height = std::max<int>(height(b->left), a->height) + 1;
        return b;
    }

    template<class Options>
    auto avl_tree<Options>::rotate_left_right(Inner *a) -> Inner * {
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
        c->left = node_ptr::first(b);
        c->right = node_ptr::first(a);

        update(*a);
        update(*b);
        update(*c);

        a->height = std::max(height(a->left), height(a->right)) + 1;
        b->height = std::max(height(b->left), height(b->right)) + 1;
        c->height = std::max(a->height, b->height) + 1;

        assert(std::get<0>(check_integrity(node_ptr::first(c))));
        return c;
    }

    template<class Options>
    auto avl_tree<Options>::rotate_right_left(Inner *a) -> Inner * {
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
        c->left = node_ptr::first(a);
        c->right = node_ptr::first(b);

        update(*a);
        update(*b);
        update(*c);

        a->height = std::max(height(a->left), height(a->right)) + 1;
        b->height = std::max(height(b->left), height(b->right)) + 1;
        c->height = std::max(a->height, b->height) + 1;

        return c;
    }

    template<class Options>
    auto avl_tree<Options>::rotate_right_right(Inner *a) -> Inner * {
        //  a             b
        //    b    ->   a   c
        //      c
        assert(a);
        Inner *b = a->right.cast().first;
        assert(b);

        a->right = b->left;
        b->left = node_ptr::first(a);

        update(*a);
        update(*b);

        a->height = std::max(height(a->left), height(a->right)) + 1;
        b->height = std::max<int>(a->height, height(b->right)) + 1;
        return b;
    }
}


#endif //ADVANCED_DATA_STRUCTURES_AVLTREE_H
