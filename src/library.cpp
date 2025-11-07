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
