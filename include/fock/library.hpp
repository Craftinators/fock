#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <numeric>

#include <Eigen/Dense>

namespace fock
{
    inline void generate_basis(const std::uint32_t num_sites, std::vector<uint64_t>& basis_states)
    {
        basis_states.resize(1 << num_sites);
        std::iota(basis_states.begin(), basis_states.end(), 0ULL);
    }

    inline void generate_basis(const std::uint32_t num_sites, const std::uint32_t num_filled_sites, std::vector<uint64_t>& basis_states)
    {
        const std::uint32_t minimum_difference = std::min(num_filled_sites, num_sites - num_filled_sites);
        std::uint64_t accumulation = 1;
        for (std::uint32_t i = 1; i <= minimum_difference; ++i)
        {
            accumulation = accumulation * (num_sites - i + 1) / i;
        }
        basis_states.reserve(accumulation);
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

    inline void generate_hamiltonian(const std::uint32_t num_sites, const std::uint32_t num_filled_sites, Eigen::MatrixXd& hamiltonian)
    {
        std::vector<uint64_t> basis_states;
        generate_basis(num_sites, num_filled_sites, basis_states);
        const std::size_t num_states = basis_states.size();

        hamiltonian = Eigen::MatrixXd::Zero(num_states, num_states); // NOLINT(*-narrowing-conversions)
        for (std::uint64_t initial_state_index = 0; initial_state_index < num_states; ++initial_state_index)
        {
            for (std::uint64_t final_state_index = 0; final_state_index < num_states; ++final_state_index)
            {
                // If states differ by a singular hop e.g., 0011 -> 0101
                const std::uint64_t differing_bitmask = basis_states[initial_state_index] ^ basis_states[final_state_index];
                if (!(std::popcount(differing_bitmask) == 2 && differing_bitmask & differing_bitmask >> 1ULL)) continue;
                hamiltonian(final_state_index, initial_state_index) = -1; // NOLINT(*-narrowing-conversions)
            }
        }
    }

    inline std::uint64_t composite_state_bitmask(const std::uint64_t subregion_state_bitmask, const std::uint64_t complement_subregion_state_bitmask, const std::uint32_t num_sites, const std::vector<std::uint32_t>& subregion_indices)
    {
        const std::size_t subregion_size = subregion_indices.size();

        std::uint64_t specified_mask = 0;
        for (std::size_t i = 0; i < subregion_size; ++i) {
            const std::uint32_t target_pos_from_right = num_sites - 1u - subregion_indices[i];
            specified_mask |= 1ULL << target_pos_from_right;
        }

        std::uint64_t merged_result = 0;
        for (std::size_t i = 0; i < subregion_size; ++i) {
            const std::uint32_t target_pos_from_right = num_sites - 1u - subregion_indices[i];
            if ((subregion_state_bitmask >> i & 1ULL) == 0) continue;
            merged_result |= 1ULL << target_pos_from_right;
        }

        std::uint64_t current_complement_bits = complement_subregion_state_bitmask;
        std::uint64_t gaps = ~specified_mask;
        while (gaps) {
            const std::uint64_t lsb = gaps & -gaps;
            if ((current_complement_bits & 1ULL) != 0) merged_result |= lsb;
            current_complement_bits >>= 1;
            gaps ^= lsb;
        }
        return merged_result;
    }

    // See https://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class for Eigen::Ref
    inline void generate_rdm(const std::uint32_t num_sites,
                         const std::uint32_t num_filled_sites,
                         const Eigen::Ref<const Eigen::VectorXd>&  eigenstate,
                         const std::vector<std::uint32_t>& subregion_indices,
                         Eigen::MatrixXd& rdm)
    {
        Eigen::MatrixXd density_matrix = eigenstate * eigenstate.transpose();

        const std::size_t num_subregion_sites = subregion_indices.size();
        const std::size_t num_complement_subregion_sites = num_sites - num_subregion_sites;

        std::vector<std::uint64_t> subregion_states;
        std::vector<std::uint64_t> complement_subregion_states;
        generate_basis(num_subregion_sites, subregion_states);
        generate_basis(num_complement_subregion_sites, complement_subregion_states);

        std::vector<std::uint64_t> composite_states;
        generate_basis(num_sites, num_filled_sites, composite_states);

        rdm = Eigen::MatrixXd::Zero(1 << num_complement_subregion_sites, 1 << num_complement_subregion_sites);
        for (std::size_t final_state_index = 0; final_state_index < 1 << num_subregion_sites; ++final_state_index)
        {
            for (std::size_t initial_state_index = 0; initial_state_index < 1 << num_subregion_sites; ++initial_state_index)
            {
                const auto x = subregion_states[final_state_index];
                const auto y= subregion_states[initial_state_index];
                if (__builtin_popcountll(x) != __builtin_popcountll(y)) continue; // if x and y dont have same count, then theres no way they form valid state
                for (std::size_t j = 0; j < 1 << num_complement_subregion_sites; ++j)
                {
                    const auto z = complement_subregion_states[j];
                    if (__builtin_popcountll(x) + __builtin_popcountll(z) != num_filled_sites || __builtin_popcountll(y) + __builtin_popcountll(z) != num_filled_sites) continue;
                    const auto a = composite_state_bitmask(x, z,
                            num_sites, subregion_indices);
                    const auto b = composite_state_bitmask(y, z,
                            num_sites, subregion_indices);
                    const auto a_idx = std::distance(composite_states.begin(), std::ranges::find(composite_states, a));
                    const auto b_idx = std::distance(composite_states.begin(), std::ranges::find(composite_states, b));
                    rdm(final_state_index, initial_state_index) += density_matrix(a_idx, b_idx);
                }
            }
        }
    }
}

#endif // FOCK_LIBRARY_H