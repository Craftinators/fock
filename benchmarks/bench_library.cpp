#include <benchmark/benchmark.h>
#include <fock/library.hpp>

static void BM_GenerateBasisStates(benchmark::State& state) {
    const int64_t n = state.range(0);
    const int64_t k = state.range(1);

    // ReSharper disable once CppDFAUnreadVariable
    for (auto _ : state) {
        auto result = generate_basis_states(n, k);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

// Helper function to generate parameter ranges
static void AddParameterRanges(benchmark::internal::Benchmark* b) {
    for (int n = 16; n < 21; ++n)
    {
        for (int k = n / 2 - 2; k < n / 2 + 3; ++k)
        {
            b->Args({n, k});
        }
    }
}

BENCHMARK(BM_GenerateBasisStates)->Apply(AddParameterRanges);

BENCHMARK_MAIN();