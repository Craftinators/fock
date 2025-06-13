#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <cstdint>
#include <vector>

std::vector<std::uint64_t> generate_basis_states(unsigned int lattice_size, unsigned int fermion_count);

std::vector<uint64_t> get_adjacent_states(uint64_t state, unsigned lattice_size, bool use_periodic_boundaries);

#endif //FOCK_LIBRARY_H