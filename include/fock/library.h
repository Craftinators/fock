#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <vector>

std::vector<std::uint64_t> generate_full_hilbert_space(std::uint32_t num_sites);

std::vector<std::uint64_t> generate_hilbert_subspace(std::uint32_t num_sites, std::uint32_t num_filled_sites,
                                                     std::uint64_t &num_states);

bool is_valid_region_composition(std::uint32_t num_filled_sites,
                                 std::uint64_t subregion_state_bitmask, std::uint64_t complement_subregion_state_bitmask);

bool are_neighbors(std::uint64_t initial_state_bitmask, std::uint64_t final_state_bitmask);

std::uint64_t get_composite_state_bitmask(const std::vector<std::uint32_t> &subregion_indices,
                                          std::uint32_t num_sites, std::uint64_t subregion_state_bitmask, std::uint64_t complement_subregion_state_bitmask);

// TODO separate this from code above

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

Eigen::Matrix<boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic>
build_hamiltonian(unsigned int num_sites, unsigned int num_filled_sites);

#endif // FOCK_LIBRARY_H