#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <vector>

std::vector<unsigned long long> generate_full_hilbert_space(unsigned int num_sites);
std::vector<unsigned long long> generate_hilbert_subspace(unsigned int num_sites, unsigned int num_filled_sites);

#endif // FOCK_LIBRARY_H