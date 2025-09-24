#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <Eigen/Sparse>

Eigen::SparseMatrix<double> build_hamiltonian(unsigned int lattice_length, unsigned int num_fermions);

#endif // FOCK_LIBRARY_H