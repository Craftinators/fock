#include <vector>
#include <numeric>
#include <stdexcept>
#include <cstdint>

inline std::uint64_t get_next_bit_permutation(const std::uint64_t current_bitmask)
{
    const std::uint64_t least_significant_bit = current_bitmask & -current_bitmask;
    const std::uint64_t next_rippled_bitmask = current_bitmask + least_significant_bit;
    return ((next_rippled_bitmask ^ current_bitmask) >> 2) / least_significant_bit | next_rippled_bitmask;
}

inline std::uint64_t get_binomial_coefficient(const std::uint32_t n, std::uint32_t k)
{
    if (k > n)
    {
        return 0;
    }

    if (k > n - k)
    {
        k = n - k;
    }

    if (k == 0)
    {
        return 1ULL;
    }

    std::uint64_t result = 1ULL;
    for (std::uint32_t i = 0; i < k; ++i)
    {
        result = result * (n - i) / (i + 1ULL);
    }
    return result;
}

inline std::uint64_t pdep_implementation(std::uint64_t source, std::uint64_t mask)
{
    std::uint64_t destination = 0;
    while (source && mask)
    {
        const unsigned num_trailing_zeros = __builtin_ctzll(mask);
        if (source & 1ULL)
        {
            destination |= 1ULL << num_trailing_zeros;
        }
        source >>= 1ULL;
        mask &= mask - 1ULL;
    }
    return destination;
}

std::vector<std::uint64_t> generate_full_hilbert_space(const std::uint32_t num_sites)
{
    if (num_sites >= 64)
    {
        throw std::runtime_error("Number of sites must be less than 64 to fit in std::uint64_t and std::vector.");
    }

    const std::size_t num_states = 1ULL << num_sites;
    std::vector<std::uint64_t> states(num_states);
    std::iota(states.begin(), states.end(), 0ULL);
    return states;
}

std::vector<std::uint64_t> generate_hilbert_subspace(const std::uint32_t num_sites, const std::uint32_t num_filled_sites,
                                                     std::uint64_t &num_states)
{
    if (num_filled_sites > num_sites)
    {
        throw std::logic_error("Number of filled sites cannot be greater than the number of sites.");
    }

    num_states = get_binomial_coefficient(num_sites, num_filled_sites);
    std::vector<std::uint64_t> states(num_states);

    if (num_filled_sites == 0)
    {
        if (num_states == 1ULL) states[0] = 0ULL;
        return states;
    }

    std::uint64_t current_state_bitmask = (1ULL << num_filled_sites) - 1ULL;
    for (std::size_t state_index = 0; state_index < num_states; ++state_index)
    {
        states[state_index] = current_state_bitmask;

        if (state_index < num_states - 1ULL)
        {
            current_state_bitmask = get_next_bit_permutation(current_state_bitmask);
        }
    }
    return states;
}

bool is_valid_region_composition(const std::uint32_t num_filled_sites,
                                 const std::uint64_t subregion_state_bitmask, const std::uint64_t complement_subregion_state_bitmask)
{
    return __builtin_popcountll(subregion_state_bitmask) + __builtin_popcountll(complement_subregion_state_bitmask) == num_filled_sites;
}

bool are_neighbors(const std::uint64_t initial_state_bitmask, const std::uint64_t final_state_bitmask)
{
    const std::uint64_t differing_bitmask = initial_state_bitmask ^ final_state_bitmask;
    // __builtin_popcountll(differing_bitmask) == 2 ensures that the states differ by exactly 2 bits
    // and differing_bitmask & differing_bitmask >> 1ULL ensures that these differing bits are adjacent.
    return __builtin_popcountll(differing_bitmask) == 2 && differing_bitmask & differing_bitmask >> 1ULL;
}

std::uint64_t get_composite_state_bitmask(const std::vector<std::uint32_t> &subregion_indices, const std::uint32_t num_sites,
                                          const std::uint64_t subregion_state_bitmask, const std::uint64_t complement_subregion_state_bitmask)
{
    std::uint64_t subregion_mask = 0;
    for (const std::uint32_t x: subregion_indices)
    {
        subregion_mask |= 1ULL << (num_sites - x - 1ULL);
    }
    const std::uint64_t complement_subregion_mask = (1ULL << num_sites) - 1ULL ^ subregion_mask;
#if defined(__BMI2__) && (defined(__x86_64__) || defined(_M_X64))
    return _pdep_u64(subregion_state_bitmask, subregion_mask) |
           _pdep_u64(complement_subregion_state_bitmask, complement_subregion_mask);
#else
    return pdep_implementation(subregion_state_bitmask, subregion_mask) |
           pdep_implementation(complement_subregion_state_bitmask, complement_subregion_mask);
#endif
}

// TODO --- Separate this from the code above! ---

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

Eigen::SparseMatrix<double> build_reduced_density_matrix(const Eigen::VectorXd state,
                                                         const std::vector<std::uint32_t> &subregion_indices,
                                                         const std::uint32_t num_sites,
                                                         const std::uint32_t num_filled_sites)
{
    // --- High level:
    // We compute rho_A = Tr_B |Psi><Psi| where |Psi> is given by `state` expressed in the
    // particle-number-conserving basis produced by generate_hilbert_subspace(num_sites,num_filled_sites).
    //
    // Indexing:
    //  - For the full system the basis states are enumerated by bitmasks (uint64_t) with exactly
    //    num_filled_sites bits set. `generate_hilbert_subspace` returns that list in `states`.
    //  - The subregion is defined by `subregion_indices` (sites belonging to A). We enumerate
    //    all local subregion bitmasks (0..2^m-1) and complement bitmasks (0..2^(n-m)-1).
    //  - We only consider pairs (subregion_mask, complement_mask) whose popcounts add to num_filled_sites.
    //
    // Algorithm:
    //  1. Recreate the full-subspace basis `states` and a map state->index (so we can look up amplitudes).
    //  2. Determine which local subregion bitmasks are *possible* (i.e. there exists a complement filling
    //     making the total particle number correct). These form the basis of rho_A (dimension = #allowed masks).
    //  3. For each complement bitmask b:
    //       - compute the vector v_b whose components are psi_{a,b} for all allowed a
    //       - add v_b * v_b^T to rho_A
    //  4. Return rho_A as an Eigen::SparseMatrix<double>.
    //
    // Complexity: loops over all complement configurations and for each builds a small vector over
    // compatible subregion configurations. Works well when subregion is small; if subregion is large
    // the memory/time grows combinatorially.

    // Basic checks
    if (num_sites >= 64)
    {
        throw std::runtime_error("Number of sites must be less than 64.");
    }

    // Reconstruct the full Hilbert subspace basis and map states -> linear indices (same as build_hamiltonian).
    std::uint64_t num_states_u64;
    const std::vector<std::uint64_t> full_states = generate_hilbert_subspace(num_sites, num_filled_sites, num_states_u64);
    const Eigen::Index num_states = static_cast<Eigen::Index>(num_states_u64);

    if (static_cast<Eigen::Index>(state.size()) != num_states)
    {
        throw std::invalid_argument("Size of `state` does not match number of states in the specified subspace.");
    }

    std::unordered_map<std::uint64_t, std::uint64_t> state_to_index;
    state_to_index.reserve(num_states_u64);
    for (std::uint64_t i = 0; i < num_states_u64; ++i)
    {
        state_to_index[full_states[i]] = i;
    }

    // Subregion size and complement size
    const std::uint32_t m = static_cast<std::uint32_t>(subregion_indices.size());
    const std::uint32_t complement_size = num_sites - m;

    // Determine allowed subregion local bitmasks:
    // A subregion mask is allowed if there exists a complement filling giving total num_filled_sites.
    // That is, kA = popcount(subregion_mask) must satisfy 0 <= num_filled_sites - kA <= complement_size.
    std::vector<std::uint64_t> allowed_sub_masks;
    allowed_sub_masks.reserve(1u << std::min<std::uint32_t>(m, 20u)); // heuristic reserve
    for (std::uint64_t sub_mask = 0ULL; sub_mask < (1ULL << m); ++sub_mask)
    {
        const std::uint32_t kA = __builtin_popcountll(sub_mask);
        if (kA <= num_filled_sites && (num_filled_sites - kA) <= complement_size)
        {
            allowed_sub_masks.push_back(sub_mask);
        }
    }

    // Map subregion local mask -> index in reduced matrix
    std::unordered_map<std::uint64_t, std::uint32_t> submask_to_index;
    submask_to_index.reserve(allowed_sub_masks.size());
    for (std::uint32_t i = 0; i < allowed_sub_masks.size(); ++i)
    {
        submask_to_index[allowed_sub_masks[i]] = i;
    }

    const Eigen::Index dimA = static_cast<Eigen::Index>(allowed_sub_masks.size());
    // Dense accumulation for simplicity and clarity, then convert to sparse at the end.
    Eigen::MatrixXd rhoA = Eigen::MatrixXd::Zero(dimA, dimA);

    // Iterate over all complement local bitmasks. For each complement mask b compute kB = popcount(b),
    // then kA = num_filled_sites - kB. If kA is in range, gather amplitudes psi_{a,b} for all subregion
    // masks a with popcount(a) == kA and accumulate outer product.
    if (complement_size == 0)
    {
        // Trivial complement: the reduced density matrix equals outer product of the state expressed in
        // the subregion basis. There is exactly one complement configuration (empty).
        // Build vector v_a = psi_{a,empty} for all allowed a and do v v^T.
        Eigen::VectorXd v = Eigen::VectorXd::Zero(dimA);
        for (std::uint32_t ai = 0; ai < allowed_sub_masks.size(); ++ai)
        {
            const std::uint64_t a_mask = allowed_sub_masks[ai];
            const std::uint64_t full_mask = get_composite_state_bitmask(subregion_indices, num_sites, a_mask, 0ULL);
            const std::uint64_t full_index = state_to_index.at(full_mask);
            v(static_cast<Eigen::Index>(ai)) = state(static_cast<Eigen::Index>(full_index));
        }
        rhoA = v * v.transpose();
    } else
    {
        const std::uint64_t complement_iter_limit = (1ULL << complement_size);
        for (std::uint64_t comp_mask = 0ULL; comp_mask < complement_iter_limit; ++comp_mask)
        {
            const std::uint32_t kB = __builtin_popcountll(comp_mask);
            // Required number of particles in A for this complement:
            if (kB > num_filled_sites) continue;
            const std::uint32_t kA_required = num_filled_sites - kB;
            if (kA_required > m) continue; // impossible

            // Gather amplitudes for all allowed subregion masks a that have popcount == kA_required
            // We'll build a small vector v where v[i] corresponds to allowed_sub_masks[i] (if its popcount matches).
            std::vector<std::pair<Eigen::Index, double> > comp_amplitudes; // (reduced_index, amplitude)
            comp_amplitudes.reserve(get_binomial_coefficient(m, kA_required));
            for (std::uint32_t ai = 0; ai < allowed_sub_masks.size(); ++ai)
            {
                const std::uint64_t a_mask = allowed_sub_masks[ai];
                if (__builtin_popcountll(a_mask) != kA_required) continue;

                // Compose full-system bitmask: place bits of `a_mask` in subregion positions and
                // bits of `comp_mask` in complement positions.
                const std::uint64_t full_mask = get_composite_state_bitmask(subregion_indices, num_sites, a_mask, comp_mask);

                // Look up the linear index in the full subspace and fetch amplitude from `state`.
                // We expect the full_mask to be present in the subspace (since popcounts sum to num_filled_sites).
                auto it = state_to_index.find(full_mask);
                if (it == state_to_index.end())
                {
                    // This should not happen for a correctly-constructed mapping; skip defensively.
                    continue;
                }
                const std::uint64_t full_index = it->second;
                const double amp = state(static_cast<Eigen::Index>(full_index));
                comp_amplitudes.emplace_back(static_cast<Eigen::Index>(ai), amp);
            }

            // Accumulate outer product v * v^T into rhoA, where v contains amplitudes for this complement.
            const std::size_t L = comp_amplitudes.size();
            for (std::size_t i = 0; i < L; ++i)
            {
                const Eigen::Index ri = comp_amplitudes[i].first;
                const double ai_val = comp_amplitudes[i].second;
                for (std::size_t j = 0; j < L; ++j)
                {
                    const Eigen::Index rj = comp_amplitudes[j].first;
                    const double aj_val = comp_amplitudes[j].second;
                    rhoA(ri, rj) += ai_val * aj_val; // conjugation if complex
                }
            }
        }
    }

    // Convert dense reduced density matrix to sparse and return.
    Eigen::SparseMatrix<double> rhoA_sparse(static_cast<Eigen::Index>(dimA), static_cast<Eigen::Index>(dimA));
    rhoA_sparse = rhoA.sparseView(); // efficient conversion
    return rhoA_sparse;
}
