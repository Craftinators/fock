#include "fock/single_particle.h"

#include <matplot/matplot.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

int main()
{
    constexpr unsigned lattice_length = 100;
    constexpr unsigned fermion_count = 50;
    constexpr unsigned maximum_subregion_size = lattice_length;

    std::vector<double> x,y;

#if FOCK_DEBUG_BUILD
    spdlog::set_level(spdlog::level::debug);
#else
    spdlog::set_level(spdlog::level::info);
#endif
    spdlog::stopwatch stopwatch;

    for (int subregion_length = 0; subregion_length < maximum_subregion_size + 1; ++subregion_length)
    {
        stopwatch.reset();

        auto subregion_correlation_matrix = fock::build_subregion_correlation_matrix(
            lattice_length,
            fermion_count,
            subregion_length);
        double entropy = fock::entanglement_entropy(subregion_correlation_matrix);

        x.push_back(subregion_length);
        y.push_back(entropy);

        auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stopwatch.elapsed()).count();
        spdlog::debug("subregion_length = {0:03d}, time_elapsed = {1:03d} ms", subregion_length, time_elapsed);
    }

    matplot::figure();
    matplot::plot(x, y);
    matplot::xlabel("Subregion length");
    matplot::ylabel("Entanglement Entropy");
    matplot::title("Entanglement Entropy for L=" +
        std::to_string(lattice_length) + ", N=" +
        std::to_string(fermion_count) + " system");
    matplot::show();
}
