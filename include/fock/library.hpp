#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <numeric>
#include <stdexcept>

#include <Eigen/Dense>

namespace fock
{
    inline void generate_basis(const std::uint32_t num_sites, std::vector<uint64_t>& basis_states)
    {
        const std::size_t num_states = 1 << num_sites;
        basis_states.resize(num_states);
        std::iota(basis_states.begin(), basis_states.end(), 0ULL);
    }

    inline void generate_basis(const std::uint32_t num_sites, const std::uint32_t num_filled_sites, std::vector<uint64_t>& basis_states)
    {
        std::uint64_t current_state_bitmask = (1ULL << num_filled_sites) - 1ULL;
        while (current_state_bitmask < 1ULL << num_sites)
        {
            basis_states.push_back(current_state_bitmask);
            // https://programmingforinsomniacs.blogspot.com/2018/03/gospers-hack-explained.html
            const std::uint64_t c = current_state_bitmask & -current_state_bitmask;
            const std::uint64_t r = current_state_bitmask + c;
            current_state_bitmask = ((r ^ current_state_bitmask) >> 2) / c | r;
        }
    }

    // See https://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class for Eigen::Ref
    inline void generate_hamiltonian(const std::uint32_t num_sites, const std::uint32_t num_filled_sites, Eigen::Ref<Eigen::MatrixXd> hamiltonian)
    {
        std::vector<uint64_t> basis_states;
        generate_basis(num_sites, num_filled_sites, basis_states);
        const std::size_t num_states = basis_states.size();

        hamiltonian = Eigen::MatrixXd::Zero(num_states, num_states); // NOLINT(*-narrowing-conversions)
        for (int initial_state_index = 0; initial_state_index < num_states; ++initial_state_index)
        {
            for (int final_state_index = 0; final_state_index < num_states; ++final_state_index)
            {
                // If states differ by a singular hop e.g., 0011 -> 0101
                const std::uint64_t differing_bitmask = basis_states[initial_state_index] ^ basis_states[final_state_index];
                if (!(__builtin_popcountll(differing_bitmask) == 2 && differing_bitmask & differing_bitmask >> 1ULL)) continue;
                hamiltonian(final_state_index, initial_state_index) = -1;
            }
        }
    }

    // See https://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class for Eigen::Ref
    inline void generate_rdm(std::uint32_t num_sites,
                         std::uint32_t num_filled_sites,
                         const Eigen::Ref<const Eigen::VectorXd>&  eigenstate,
                         const std::vector<std::uint32_t>& subregion_indices,
                         Eigen::Ref<Eigen::MatrixXd> rdm)
    {
        // TODO: implementation
    }
}

#endif // FOCK_LIBRARY_H