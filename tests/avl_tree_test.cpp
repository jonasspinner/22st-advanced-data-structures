#include "AVLTree.h"

using size_type = int32_t;
struct InnerNode {
    size_type size{};
    size_type ones{};
};
struct LeafNode {
    size_type size{};
    size_type ones{};
    std::array<uint64_t, 8> bits{};
    
    bool access(size_type i) const {
        auto j = i / bits.size();
        auto k = i % bits.size();
        return bits[j] >> k & 1;
    }
    
    size_type rank(size_type i) const {
        return 0;
    }
    size_type select(size_type i) const {
        return 0;
    }
};
struct InnerNodeUpdater {
    template <class Left, class Right>
    void operator()(InnerNode &node, const Left &left, const Right &right) const {
        node.size = left.size + right.size;
        node.ones = left.ones + right.ones;
    }
};
struct LeafInserter {
};
using Options = ads::avl_tree_options<InnerNode, LeafNode, InnerNodeUpdater>;
using avl_tree = ads::avl_tree<Options>;


struct AccessNavigator {
    size_type i{};
    bool result{};
    
    template <class Left, class Right>
    ads::direction visit_inner(const Left &left, const Right &) {
        if (i < left.size) {
            return ads::direction::left;
        } else {
            i -= left.size;
            return ads::direction::right;
        }
    }

    void visit_leaf(const LeafNode &node) {
        result += node.access(i);
    }
};

struct RankNavigator {
    size_type i{};
    size_type result{};

    template <class Left, class Right>
    ads::direction visit_inner(const Left &left, const Right &) {
        if (i < left.size) {
            return ads::direction::left;
        } else {
            i -= left.size;
            result += left.ones;
            return ads::direction::right;
        }
    }

    void visit_leaf(const LeafNode &node) {
        result += node.rank(i);
    }
};

struct SelectNavigator {
    size_type i{};
    size_type result{};

    template <class Left, class Right>
    ads::direction visit_inner(const Left &left, const Right &) {
        if (i <= left.ones) {
            return ads::direction::left;
        } else {
            i -= left.ones;
            result += left.size;
            return ads::direction::right;
        }
    }

    void visit_leaf(const LeafNode &node) {
        result += node.select(i);
    }
};

static_assert(std::is_invocable_r_v<ads::direction, decltype(&SelectNavigator::visit_inner<InnerNode, LeafNode>), SelectNavigator&, InnerNode, LeafNode>);


int main() {
    avl_tree tree;

    AccessNavigator navigator{10};
    tree.query(navigator);
}
