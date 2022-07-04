#include <benchmark/benchmark.h>
#include <random>
#include <algorithm>

#include "NaiveDynamicBitVector.h"
#include "DynamicBitVector.h"


template <class BV>
static void BM_InsertN(benchmark::State& state) {
    size_t n = state.range(0);
    for (auto _ : state) {
        std::minstd_rand gen;
        std::bernoulli_distribution bit_dist;

        BV bv;
        for (size_t i = 0; i < n; ++i) {
            std::uniform_int_distribution<size_t> index_dist(0, bv.size());

            bv.insert(index_dist(gen), bit_dist(gen));
        }
        benchmark::DoNotOptimize(bv.rank(bv.size(), true));
    }
}

BENCHMARK(BM_InsertN<ads::NaiveDynamicBitVector>)->Arg(1 << 5)->Arg(1 << 10)->Arg(1 << 15);
BENCHMARK(BM_InsertN<ads::DynamicBitVector<uint8_t, 2>>)->Arg(1 << 5)->Arg(1 << 10)->Arg(1 << 15);
BENCHMARK(BM_InsertN<ads::DynamicBitVector<uint64_t, 2>>)->Arg(1 << 5)->Arg(1 << 10)->Arg(1 << 15);
BENCHMARK(BM_InsertN<ads::DynamicBitVector<uint64_t, 4>>)->Arg(1 << 5)->Arg(1 << 10)->Arg(1 << 15);
BENCHMARK(BM_InsertN<ads::DynamicBitVector<uint64_t, 8>>)->Arg(1 << 5)->Arg(1 << 10)->Arg(1 << 15);
BENCHMARK(BM_InsertN<ads::DynamicBitVector<uint64_t, 16>>)->Arg(1 << 5)->Arg(1 << 10)->Arg(1 << 15);


template <class BV>
static void BM_Insert(benchmark::State& state) {
    int n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BV bv;
    for (int i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        bv.insert(index_dist(gen), bit_dist(gen));
    }
    for (auto _ : state) {
        std::uniform_int_distribution<int> index_dist(0, bv.size());
        bv.insert(index_dist(gen), bit_dist(gen));
        if (bv.size() == n + 128) {
            state.PauseTiming();
            for (int i = 0; i < 128; ++i) {
                bv.remove(std::clamp(index_dist(gen), 0, n-1));
            }
            state.ResumeTiming();
        }
    }
}
BENCHMARK(BM_Insert<ads::NaiveDynamicBitVector>)->RangeMultiplier(2)->Range(8, 1<<15);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 32>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 64>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Insert<ads::DynamicBitVector<uint64_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<24);


template <class BV>
static void BM_Remove(benchmark::State& state) {
    size_t n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BV bv;
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        bv.insert(index_dist(gen), bit_dist(gen));
    }
    for (auto _ : state) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        if (bv.size() == n) {
            state.PauseTiming();
            for (int i = 0; i < 128; ++i) {
                bv.insert(index_dist(gen), bit_dist(gen));
            }
            state.ResumeTiming();
        }
        bv.remove(index_dist(gen));
    }
}
BENCHMARK(BM_Remove<ads::NaiveDynamicBitVector>)->RangeMultiplier(2)->Range(8, 1<<15);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 32>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 64>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Remove<ads::DynamicBitVector<uint64_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<24);



template <class BV>
static void BM_RemoveMid(benchmark::State& state) {
    size_t n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BV bv;
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        bv.insert(index_dist(gen), bit_dist(gen));
    }
    for (auto _ : state) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());
        if (bv.size() == n) {
            state.PauseTiming();
            for (int i = 0; i < 128; ++i) {
                bv.insert(index_dist(gen), bit_dist(gen));
            }
            state.ResumeTiming();
        }
        bv.remove(bv.size() / 2);
    }
}
BENCHMARK(BM_RemoveMid<ads::NaiveDynamicBitVector>)->RangeMultiplier(2)->Range(8, 1<<15);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 32>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 64>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_RemoveMid<ads::DynamicBitVector<uint64_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<24);

template <class BV>
static void BM_Rank(benchmark::State& state) {
    size_t n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BV bv;
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());

        bv.insert(index_dist(gen), bit_dist(gen));
    }

    std::uniform_int_distribution<size_t> index_dist(0, bv.size());

    for (auto _ : state) {
        benchmark::DoNotOptimize(bv.rank(index_dist(gen), bit_dist(gen)));
    }
}

BENCHMARK(BM_Rank<ads::NaiveDynamicBitVector>)->RangeMultiplier(2)->Range(8, 1<<15);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 32>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 64>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Rank<ads::DynamicBitVector<uint64_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<24);


template <class BV>
static void BM_Select(benchmark::State& state) {
    size_t n = state.range(0);
    std::minstd_rand gen;
    std::bernoulli_distribution bit_dist;

    BV bv;
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> index_dist(0, bv.size());

        bv.insert(index_dist(gen), bit_dist(gen));
    }

    auto num_ones = bv.rank(bv.size(), true);
    std::uniform_int_distribution<size_t> index_dist(1, num_ones);

    for (auto _ : state) {
        benchmark::DoNotOptimize(bv.select(index_dist(gen), true));
    }
}

BENCHMARK(BM_Select<ads::NaiveDynamicBitVector>)->RangeMultiplier(2)->Range(8, 1<<15);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint8_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint8_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 2>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 4>>)->RangeMultiplier(2)->Range(8, 1<<22);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 8>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 16>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 32>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 64>>)->RangeMultiplier(2)->Range(8, 1<<24);
BENCHMARK(BM_Select<ads::DynamicBitVector<uint64_t, 128>>)->RangeMultiplier(2)->Range(8, 1<<24);

BENCHMARK_MAIN();