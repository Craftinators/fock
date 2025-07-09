#include "fock/single_particle.h"

#include <matplot/matplot.h>
#include <spdlog/spdlog.h>

int main()
{
    constexpr unsigned maxA = 100;
    spdlog::info("Building entropy curve...");
    std::vector<double> A, S;
    for (unsigned a = 1; a <= maxA; ++a) {
        auto C = fock::build_subregion_correlation_matrix(100, 50, a);
        double s = fock::entanglement_entropy(C);
        A.push_back(a);
        S.push_back(s);
        spdlog::info("A = {}, S = {}", a, s);
    }
    matplot::figure();
    matplot::plot(A, S);
    matplot::xlabel("Subregion size A");
    matplot::ylabel("Entanglement entropy S");
    matplot::title("S(A) for L=100, N=50");
    matplot::show();
}
