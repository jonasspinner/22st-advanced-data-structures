#include <benchmark/benchmark.h>
#include <random>

#include "BP.h"

template <class BlockType, size_t NumBlocks>
static void BM_ForwardSearch(benchmark::State& state) {
    using block_type = BlockType;
    using size_type = int32_t;
    constexpr size_type num_blocks = NumBlocks;
    std::uniform_int_distribution<BlockType> block_dist(std::numeric_limits<block_type>::min(), std::numeric_limits<block_type>::max());

    using Leaf = typename ads::BP<block_type, num_blocks>::Leaf;

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
        if (idx < 0 || leaf.size() <= idx) { num_fails++; }
    }
}
BENCHMARK(BM_ForwardSearch<uint64_t, 4>);
BENCHMARK(BM_ForwardSearch<uint64_t, 8>);
BENCHMARK(BM_ForwardSearch<uint64_t, 16>);
BENCHMARK(BM_ForwardSearch<uint64_t, 32>);
BENCHMARK(BM_ForwardSearch<uint64_t, 64>);
BENCHMARK(BM_ForwardSearch<uint64_t, 128>);


template <class BP>
BP build_random_inserts_bp_tree(int n) {
    std::minstd_rand gen;
    BP bp;
    for (int i = 0; i < n; ++i) {
        std::uniform_int_distribution<int> vertex_dist(0, bp.num_nodes() - 1);
        bp.insertchild(vertex_dist(gen), 1, 0);
    }
    return bp;
}



template <class BP>
static void BM_InsertChild(benchmark::State& state) {
    int n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BP bp = build_random_inserts_bp_tree<BP>(n);

    for (auto _ : state) {
        std::uniform_int_distribution<int> vertex_dist(0, bp.num_nodes() - 1);
        bp.insertchild(vertex_dist(gen), 1, 0);
        if (bp.num_nodes() == n + 128) {
            state.PauseTiming();
            for (int i = 0; i < 128; ++i) {
                auto v = std::clamp(vertex_dist(gen), 1, n-1);
                bp.deletenode(v);
            }
            state.ResumeTiming();
        }
    }
}
//BENCHMARK(BM_Insert<ads::BP<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
//BENCHMARK(BM_Insert<ads::BP<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_InsertChild<ads::BP<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_InsertChild<ads::BP<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_InsertChild<ads::BP<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_InsertChild<ads::BP<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<22);

template <class BP>
static void BM_DeleteNode(benchmark::State& state) {
    int n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BP bp = build_random_inserts_bp_tree<BP>(n);

    for (auto _ : state) {
        std::uniform_int_distribution<int> vertex_dist(1, bp.num_nodes() - 1);
        if (bp.num_nodes() == n) {
            state.PauseTiming();
            for (int i = 0; i < 128; ++i) {
                bp.insertchild(vertex_dist(gen), 1, 0);
            }
            state.ResumeTiming();
        }
        bp.deletenode(vertex_dist(gen));
    }
}

//BENCHMARK(BM_DeleteNode<ads::BP<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
//BENCHMARK(BM_DeleteNode<ads::BP<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_DeleteNode<ads::BP<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<20);
BENCHMARK(BM_DeleteNode<ads::BP<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<20);
BENCHMARK(BM_DeleteNode<ads::BP<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<20);
BENCHMARK(BM_DeleteNode<ads::BP<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<20);
BENCHMARK(BM_DeleteNode<ads::BP<uint64_t, 32>>)->RangeMultiplier(2)->Range(8, 1<<20);
BENCHMARK(BM_DeleteNode<ads::BP<uint64_t, 64>>)->RangeMultiplier(2)->Range(8, 1<<20);

BENCHMARK_MAIN();