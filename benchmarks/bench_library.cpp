#include <benchmark/benchmark.h>
#include <fock/library.h>

static void BM_Gosper(benchmark::State& state) {
    int n = state.range(0);
    int k = state.range(1);

    for (auto _ : state) {
        auto result = generate_basis_states(n, k);
        benchmark::DoNotOptimize(result);
    }
}

// Register benchmarks with different parameters
BENCHMARK(BM_Gosper)
    ->Args({10, 3})
    ->Args({15, 5})
    ->Args({20, 8})
    ->Args({25, 10});

BENCHMARK_MAIN();
