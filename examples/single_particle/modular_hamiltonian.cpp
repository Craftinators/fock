#include "fock/single_particle.h"

#include <matplot/matplot.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

int main()
{
    constexpr unsigned lattice_length = 24;
    constexpr unsigned fermion_count = 12;
    constexpr unsigned subregion_length = 12;

#if FOCK_DEBUG_BUILD
    spdlog::set_level(spdlog::level::debug);
#else
    spdlog::set_level(spdlog::level::info);
#endif
    spdlog::stopwatch stopwatch;

    stopwatch.reset();
    const auto subregion_correlation_matrix = fock::build_subregion_correlation_matrix(
        lattice_length,
        fermion_count,
        subregion_length);
    auto time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stopwatch.elapsed()).count();
    spdlog::debug("Built subregion_correlation_matrix in {0:d} us", time_elapsed);

    stopwatch.reset();
    auto modular_hamiltonian = fock::build_modular_hamiltonian(subregion_correlation_matrix);
    spdlog::debug("Built modular_hamiltonian in {0:d} us", time_elapsed);

    // matplotplusplus does not support Eigen matrices out of the box, so sad :(
    std::vector matrix_elements(subregion_length, std::vector<double>(subregion_length));
    for (int i = 0; i < subregion_length; ++i)
    {
        for (int j = 0; j < subregion_length; ++j)
        {
            matrix_elements[i][j] = modular_hamiltonian(i, j);
        }
    }

    matplot::figure();
    matplot::imagesc(matrix_elements);
    matplot::colorbar();
    matplot::title("Modular Hamiltonian h_{ij} for L=" +
        std::to_string(lattice_length) + ", N=" +
        std::to_string(fermion_count) + ", A=" +
        std::to_string(subregion_length) + " system");
    matplot::show();
}