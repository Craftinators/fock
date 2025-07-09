#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <Eigen/Dense>

namespace fock
{
    Eigen::MatrixXd build_subregion_correlation_matrix(unsigned lattice_length, unsigned fermion_count, unsigned subregion_length);

    double entanglement_entropy(const Eigen::MatrixXd& subregion_correlation_matrix);

    Eigen::MatrixXd build_modular_hamiltonian(const Eigen::MatrixXd& subregion_correlation_matrix);
}

#endif //FOCK_LIBRARY_H
