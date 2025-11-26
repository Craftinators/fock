#include <unordered_map>
#import <Eigen/Sparse>

Eigen::SparseMatrix<double> build_hamiltonian(const std::uint32_t num_sites, const std::uint32_t num_filled_sites)
{
    std::uint64_t num_states;
    const std::vector<std::uint64_t> states = generate_hilbert_subspace(num_sites, num_filled_sites, num_states);

    std::unordered_map<std::uint64_t, std::uint64_t> state_to_index;
    state_to_index.reserve(num_states);
    for (std::uint64_t state_index = 0; state_index < num_states; ++state_index)
    {
        state_to_index[states[state_index]] = state_index;
    }

    std::vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(2 * num_states * num_filled_sites * (num_sites - num_filled_sites) / num_sites);

    for (std::uint64_t initial_state_index = 0; initial_state_index < num_states; ++initial_state_index)
    {
        const std::uint64_t initial_state_bitmask = states[initial_state_index];
        for (std::uint64_t current_site = 0; current_site < num_sites - 1; ++current_site)
        {
            const bool is_current_site_occupied = initial_state_bitmask >> current_site & 1;
            if (const bool is_next_site_occupied = initial_state_bitmask >> (current_site + 1) & 1;
                is_current_site_occupied && !is_next_site_occupied)
            {
                const std::uint64_t final_state_bitmask = initial_state_bitmask & ~(1ULL << current_site) | 1ULL << (current_site + 1);
                const std::uint64_t final_state_index = state_to_index.at(final_state_bitmask);

                triplets.emplace_back(initial_state_index, final_state_index, -1.0);
                triplets.emplace_back(final_state_index, initial_state_index, -1.0);
            }
        }
    }

    // TODO Revisit this warning
    Eigen::SparseMatrix<double> hamiltonian(num_states, num_states); // NOLINT(*-narrowing-conversions)
    hamiltonian.setFromTriplets(triplets.begin(), triplets.end());
    return hamiltonian;
}
