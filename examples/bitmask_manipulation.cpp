#include <iostream>

#include "fock/library.hpp"

int main()
{
    constexpr std::uint32_t num_sites = 4;
    constexpr std::uint32_t num_filled_sites = 2;
    Eigen::MatrixXd hamiltonian;

    fock::generate_hamiltonian(num_sites, num_filled_sites, hamiltonian);
    std::cout << hamiltonian;
}
