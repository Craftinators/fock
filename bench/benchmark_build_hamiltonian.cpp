#include <benchmark/benchmark.h>
#include <Eigen/Sparse>
#include "fock/library.h"

// ReSharper disable once CppInconsistentNaming
static void BM_BuildHamiltonian(benchmark::State &state) {
    const auto lattice_length = static_cast<unsigned int>(state.range(0));
    const auto num_fermions = static_cast<unsigned int>(state.range(1));

    // ReSharper disable once CppDFAUnreadVariable
    for (auto _: state) {
        auto mat = build_hamiltonian(lattice_length, num_fermions);
        benchmark::DoNotOptimize(mat);
    }
}

#define REGISTER_PAIR(L,N) ->Args({L,N})

#define REGISTER_PAIRS() REGISTER_PAIR(24, 12)

BENCHMARK(BM_BuildHamiltonian)
REGISTER_PAIRS()
->
Unit (benchmark::kMillisecond);

BENCHMARK_MAIN();