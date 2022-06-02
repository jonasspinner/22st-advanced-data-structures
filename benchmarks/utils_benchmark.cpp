#include <benchmark/benchmark.h>
#include <random>

#include "utils.h"

static void BM_Select(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(ads::select(14514284786278117030ULL, 30));
    }
}
BENCHMARK(BM_Select);

static void BM_Select64_V1(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(ads::select64_v1(14514284786278117030ULL, 30));
    }
}
BENCHMARK(BM_Select64_V1);

static void BM_Select64_V2(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(ads::select64_v2(14514284786278117030ULL, 30));
    }
}
BENCHMARK(BM_Select64_V2);


static void BM_Select64_V3(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(ads::select64_v3(14514284786278117030ULL, 30));
    }
}
BENCHMARK(BM_Select64_V3);

static void BM_Random_Select(benchmark::State& state) {
    std::minstd_rand gen;
    std::uniform_int_distribution<uint64_t> dist;
    for (auto _ : state) {
        auto word = dist(gen);
        std::uniform_int_distribution<int> index_dist(1, std::popcount(word));
        benchmark::DoNotOptimize(ads::select(word, index_dist(gen)));
    }
}
BENCHMARK(BM_Random_Select);

static void BM_Random_Select64_V1(benchmark::State& state) {
    std::minstd_rand gen;
    std::uniform_int_distribution<uint64_t> dist;
    for (auto _ : state) {
        auto word = dist(gen);
        std::uniform_int_distribution<int> index_dist(1, std::popcount(word));
        benchmark::DoNotOptimize(ads::select64_v1(word, index_dist(gen)));
    }
}
BENCHMARK(BM_Random_Select64_V1);

static void BM_Random_Select64_V2(benchmark::State& state) {
    std::minstd_rand gen;
    std::uniform_int_distribution<uint64_t> dist;
    for (auto _ : state) {
        auto word = dist(gen);
        std::uniform_int_distribution<int> index_dist(1, std::popcount(word));
        benchmark::DoNotOptimize(ads::select64_v2(word, index_dist(gen)));
    }
}
BENCHMARK(BM_Random_Select64_V2);


static void BM_Random_Select64_V3(benchmark::State& state) {
    std::minstd_rand gen;
    std::uniform_int_distribution<uint64_t> dist;
    for (auto _ : state) {
        auto word = dist(gen);
        std::uniform_int_distribution<int> index_dist(1, std::popcount(word));
        benchmark::DoNotOptimize(ads::select64_v3(word, index_dist(gen)));
    }
}
BENCHMARK(BM_Random_Select64_V3);

static void BM_Random_Baseline(benchmark::State& state) {
    std::minstd_rand gen;
    std::uniform_int_distribution<uint64_t> dist;
    for (auto _ : state) {
        auto word = dist(gen);
        std::uniform_int_distribution<int> index_dist(1, std::popcount(word));
        benchmark::DoNotOptimize(word + index_dist(gen));
    }
}
BENCHMARK(BM_Random_Baseline);

static void BM_Popcount(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(std::popcount(14514284786278117030ULL));
    }
}
BENCHMARK(BM_Popcount);

BENCHMARK_MAIN();
