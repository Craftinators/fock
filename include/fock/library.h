#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <cstdint>
#include <vector>

std::vector<std::uint64_t> generate_basis_states(unsigned int lattice_size, unsigned int fermion_count);

#endif //FOCK_LIBRARY_H