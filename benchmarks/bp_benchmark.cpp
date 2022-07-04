#include <benchmark/benchmark.h>
#include <random>

#include "BP2.h"

template <class BlockType, size_t NumBlocks>
static void BM_ForwardSearch(benchmark::State& state) {
    using block_type = BlockType;
    using size_type = int32_t;
    constexpr size_type num_blocks = NumBlocks;
    std::uniform_int_distribution<BlockType> block_dist(std::numeric_limits<block_type>::min(), std::numeric_limits<block_type>::max());
    std::array<block_type, num_blocks> bits{};

    using Leaf = typename ads::BP2<block_type, num_blocks>::Leaf;

    std::minstd_rand gen;

    Leaf leaf;

    for (int i = 0; i < Leaf::max_num_bits; ++i) {
        leaf.push_back(i < Leaf::max_num_bits / 2);
    }
    std::uniform_int_distribution<size_type> index_dist(0, leaf.size() - 1);
    std::uniform_int_distribution<size_type> d_dist(-64, 64);

    size_type num_fails{}, num_tries{};
    for (auto _ : state) {
        auto i = index_dist(gen) / 2;
        //auto d = d_dist(gen);
        //benchmark::DoNotOptimize(i + d);
        // benchmark::DoNotOptimize(leaf.rank1(i));
        //benchmark::DoNotOptimize(leaf.forward_search(i, d));
        //benchmark::DoNotOptimize(leaf.calculate_min_excess());
        //benchmark::DoNotOptimize(leaf.calculate_excess_info());
        auto [idx, d_prime] = leaf.forward_search(i, 0);
        benchmark::DoNotOptimize(idx);
        num_tries++;
        if (idx == -1 || idx == leaf.size()) { num_fails++; }
    }
    std::cout << num_tries << " " << num_fails << std::endl;
}
BENCHMARK(BM_ForwardSearch<uint64_t, 64>);

BENCHMARK_MAIN();