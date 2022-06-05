#include <random>

#include "io.h"
#include "DynamicBitVector.h"

int main(int argc, char *argv[]) {
    size_t base_size{1 << 30};
    std::stringstream ss;
    if (argc == 2) {
        ss << argv[1];
        ss >> base_size;
    }

    if (ss.bad() || base_size < 16 || base_size > (1 << 30))
        throw std::runtime_error("");

    size_t grow_size = base_size / 4;

    std::vector<ads::bv_operation> operations;

    std::mt19937_64 gen;
    std::bernoulli_distribution bit_dist;

    ads::DynamicBitVector<> bv;

    std::cout << base_size << "\n";
    for (size_t j = 0; j < base_size; ++j) {
        auto b = bit_dist(gen);
        bv.insert(bv.size(), b);
        std::cout << b << "\n";
    }

    for (size_t j = 0; j < grow_size; ++j) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        auto i = index_dist(gen);
        auto b = bit_dist(gen);
        bv.insert(i, b);
        operations.push_back(ads::bv_operation{ads::bv_operation_kind::insert, i , b});
    }

    for (size_t j = 0; j < grow_size; ++j) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size() - 1);
        auto i = index_dist(gen);
        bv.remove(i);
        operations.push_back(ads::bv_operation{ads::bv_operation_kind::remove, i , false});
    }

    for (size_t j = 0; j < grow_size; ++j) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        auto i = index_dist(gen);
        auto b = bit_dist(gen);
        operations.push_back(ads::bv_operation{ads::bv_operation_kind::rank, i , b});
    }

    auto num_ones = bv.rank(bv.size(), true);
    auto num_zeros = bv.rank(bv.size(), false);

    for (size_t j = 0; j < grow_size; ++j) {
        auto b = bit_dist(gen);
        auto num_b = b ? num_ones : num_zeros;
        std::uniform_int_distribution<size_t> index_dist(1, num_b);
        auto i = index_dist(gen);
        operations.push_back(ads::bv_operation{ads::bv_operation_kind::select, i , b});
    }

    for (auto op : operations) {
        std::cout << op << "\n";
    }
}