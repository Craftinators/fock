#include <iostream>

#include "fock/library.hpp"

int main()
{
    constexpr std::uint32_t num_sites = 4;
    constexpr std::uint32_t num_filled_sites = 2;

    Eigen::MatrixXd hamiltonian;
    fock::generate_hamiltonian(num_sites, num_filled_sites, hamiltonian);
    std::cout << hamiltonian << "\n\n";

    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hamiltonian);
    const Eigen::VectorXd ground_state = es.eigenvectors().col(0);
    std::cout << ground_state << "\n\n";

    Eigen::MatrixXd rdm;
    fock::generate_rdm(num_sites, num_filled_sites, ground_state, {0, 1}, rdm);
    std::cout << rdm;
}
