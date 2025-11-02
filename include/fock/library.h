#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <vector>

std::vector<unsigned long long> generate_full_hilbert_space(unsigned int num_sites);
std::vector<unsigned long long> generate_hilbert_subspace(unsigned int num_sites, unsigned int num_filled_sites, unsigned long long &num_states);

// TODO separate this from code above

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

Eigen::Matrix<boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic>
build_hamiltonian(unsigned int num_sites, unsigned int num_filled_sites);

#endif // FOCK_LIBRARY_H