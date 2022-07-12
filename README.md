# Fortgeschrittene Datenstrukturen

+ [Vorlesungswebseite](https://algo2.iti.kit.edu/4264.php)
+ [Aufgabenbeschreibung](https://algo2.iti.kit.edu/download/project_ss22.pdf)
+ [Wikipedia AVL tree](https://en.wikipedia.org/wiki/AVL_tree)

### Bitvektoren in Blattknoten

`SmallDynamicBitVector<>` und `SmallStaticBitVector<>`


## Aufgabe 1: Dynamischer Bitvektor `ads::DynamicBitVector<>`

    ## Inner and Leaf Node ##

    ## Query operations on nodes ##

        bool access(node, i)
        bool access(inner, i)
        void flip(node, i)
        void flip(inner, i)
        void rank1(node, i, rank&)
        void rank1(inner, i, rank&)
        void select1(node, i, result&)
        void select1(inner, i, result&)
        void select0(node, i, result&)
        void select0(inner, i, result&)

    ## Tree modifying operations ##
    Operations that insert or remove bits from the tree.

        Inner* split(leaf, high_bit)
        {node, changed} insert(node, i, b)
        {node, deleted_bit} remove(node, i)
        {node, deleted_bit} remove_second_level(inner, i)
            -- Special case for nodes with height two.
               Allows finding neighboring leaves for merging.
        {node, deleted_bit} remove_first_level(inner, i)
            -- Special case for nodes with height one.
               Allows finding neighboring leaves for merging.
        Leaf* merge(inner)

    ### Node info dispatch ###

        int height(node)

    ## AVL Tree balance operations ##

        int balance_factor(node)
        int balance_factor(inner)
        Inner* balance(inner)
            -- Checks whether balancing is needed and applies one of the four rotations
        Inner* rotate_left_left(inner)
        Inner* rotate_left_right(inner)
        Inner* rotate_right_left(inner)
        Inner* rotate_right_right(inner)

    ## Node deletion ##
    Used in destructor.

        void delete_node(NodeHandle)
        void delete_node(Inner*)
        void delete_node(Leaf*)

    ## Debug checks ##

    ## Required bits upperbound ##

    ## Constructor/Destructor ##

    ## BV operations ##

    ### BV modifying operations ##

        void insert(i, b)
        void remove(i)
        void flip(i)

    ### BV query operations ##
    
        bool access(i)
        size_type rank(i, b)
        size_type select(i, b)
        size_type size()


    ## Misc ##

        void clear()    
        std::ostream << overload
        size_type required_bits_upperbound()


## Aufgabe 2: Balanced Parantheses `ads::BP<>`

    ## Inner and Leaf Node ##

    ## Query operations on nodes ##

        bool access(node, i)
        bool access(inner, i)
        void rank1(node, i, rank&)
        void rank1(inner, i, rank&)
        void select1(node, i, result&)
        void select1(inner, i, result&)
        void select0(node, i, result&)
        void select0(inner, i, result&)

    ## Tree modifying operations ##
    Operations that insert or remove bits from the tree.

        Inner* split(leaf, high_bit)
        {node, changed} insert(node, i, b)
        {node, deleted_bit} remove(node, i)
        {node, deleted_bit} remove_second_level(inner, i)
            -- Special case for nodes with height two.
               Allows finding neighboring leaves for merging.
        {node, deleted_bit} remove_first_level(inner, i)
            -- Special case for nodes with height one.
               Allows finding neighboring leaves for merging.
        Leaf* merge(inner)

    ### Node info dispatch ###

        int height(node)
        size_type total_excess(node)
        size_type min_excess(node)
        void update_excess_info_from_children_info(inner)

    ## AVL Tree balance operations ##

        int balance_factor(node)
        int balance_factor(inner)
        Inner* balance(inner)
            -- Checks whether balancing is needed and applies one of the four rotations
        Inner* rotate_left_left(inner)
        Inner* rotate_left_right(inner)
        Inner* rotate_right_left(inner)
        Inner* rotate_right_right(inner)

    ## Node deletion ##
    Used in destructor.

        void delete_node(NodeHandle)
        void delete_node(Inner*)
        void delete_node(Leaf*)

    ## Debug checks ##

    ## Required bits upperbound ##

    ## BV operations ##

    ### BV modifying operations ##
    These operations are allowed to disregard BP invariants,
    i.e. an even and balanced number of bits.

        void insert(i, b)
        void remove(i)

    ### BV query operations ##
    
        bool access(i)
        size_type rank(i, b)
        size_type select(i, b)
        size_type size()

    ### BP query operations ###

        size_type find_close(i)
        size_type find_open(i)
        size_type excess(i)
        size_type enclose(i)

    ### BV map operation ###

        void preorder_map(f)
        void preorder_map(node, f)

    ## Constructor/Destructor ##

    ### Supported tree operations ###
            
        void insertchild(v, i, k)
        void deletenode(v)
        size_type ith_child(v, i)
        size_type parent(v)
        size_type subtree_size(v)

    ## Misc ##
    
        std::ostream << overload
        size_type required_bits_upperbound()
        size_type num_nodes()
        std::vector<size_type> preorder_out_degree_sequence()
            -- Uses preorder_map to iterate over all bits and calculate the out degree sequence
               in O(n) time.


Unterschiede von BV and BP:
```bash
diff src/DynamicBitVector.h src/BP.h
```