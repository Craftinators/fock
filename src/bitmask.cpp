#include <vector>
#include <numeric>
#include <stdexcept>

#if defined(__BMI2__) && (defined(__x86_64__) || defined(_M_X64))
#include <immintrin.h>
#endif

/// @brief Computes the next lexicographical bit permutation with the same Hamming weight.
/// @details This algorithm is commonly known as "Gosper's Hack". It finds the next integer
///          that has the exact same number of set bits (1s) as the current bitmask.
/// @param current_bitmask The current state or bit pattern.
/// @return The next integer in the sequence with the same number of set bits.
inline std::uint64_t next_bit_permutation(const std::uint64_t current_bitmask)
{
    const std::uint64_t least_significant_bit = current_bitmask & -current_bitmask;
    const std::uint64_t next_rippled_bitmask = current_bitmask + least_significant_bit;
    return ((next_rippled_bitmask ^ current_bitmask) >> 2) / least_significant_bit | next_rippled_bitmask;
}

/// @brief Calculates the binomial coefficient "n choose k".
/// @details Computes the number of ways to choose k elements from a set of n elements.
///          Includes optimizations for symmetry and overflow protection during calculation.
/// @param n The total number of items.
/// @param k The number of items to choose.
/// @return The number of combinations, or 0 if k > n.
inline std::uint64_t binomial_coefficient(const std::uint32_t n, std::uint32_t k)
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

/// @brief Software fallback implementation of the Parallel Bit Deposit (PDEP) instruction.
/// @details Deposits bits from the low bits of the source operand into the destination operand
///          at locations specified by the mask. This is used when hardware BMI2 support is unavailable.
/// @param source The contiguous bits to deposit.
/// @param mask The bitmask specifying where to place the bits.
/// @return The resulting bitmask with bits distributed according to the mask.
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

/// @brief Generates the full Hilbert space for a given number of sites.
/// @details Creates a vector containing all integers from 0 to 2^N-1.
/// @param num_sites The number of quantum sites (or bits).
/// @return A vector containing all possible states.
/// @throws std::runtime_error if num_sites >= 64 (cannot fit in uint64).
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

/// @brief Generates a subspace of the Hilbert space with a fixed particle number.
/// @details Generates all bitmasks of length `num_sites` that strictly contain
///          `num_filled_sites` set bits (1's).
/// @param num_sites The total number of sites.
/// @param num_filled_sites The number of particles (set bits).
/// @param[out] num_states Reference to store the total count of generated states.
/// @return A vector containing the valid states in the subspace.
/// @throws std::logic_error if num_filled_sites > num_sites.
std::vector<std::uint64_t> generate_hilbert_subspace(const std::uint32_t num_sites, const std::uint32_t num_filled_sites,
                                                     std::uint64_t &num_states)
{
    if (num_filled_sites > num_sites)
    {
        throw std::logic_error("Number of filled sites cannot be greater than the number of sites.");
    }

    num_states = binomial_coefficient(num_sites, num_filled_sites);
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
            current_state_bitmask = next_bit_permutation(current_state_bitmask);
        }
    }
    return states;
}

/// @brief Validates if two subregions combine to form the correct total particle count.
/// @details Checks if the sum of the population counts of the two subregions
///          equals the expected total filled sites.
/// @param num_filled_sites The expected total number of particles.
/// @param subregion_state_bitmask Bitmask for the first subregion.
/// @param complement_subregion_state_bitmask Bitmask for the second subregion.
/// @return True if the particle counts sum up correctly, false otherwise.
bool is_valid_region_composition(const std::uint32_t num_filled_sites,
                                 const std::uint64_t subregion_state_bitmask, const std::uint64_t complement_subregion_state_bitmask)
{
    return __builtin_popcountll(subregion_state_bitmask) + __builtin_popcountll(complement_subregion_state_bitmask) == num_filled_sites;
}

/// @brief Determines if two states are adjacent neighbors.
/// @details Checks if two states differ by exactly one particle hop between *adjacent* sites.
///          Specifically, the Hamming distance must be 2 (one bit flipped 0->1, another 1->0),
///          and the differing bits must be physically adjacent in the bitmask representation.
/// @param initial_state_bitmask The starting state.
/// @param final_state_bitmask The target state.
/// @return True if the states represent a nearest-neighbor hop.
bool are_neighbors(const std::uint64_t initial_state_bitmask, const std::uint64_t final_state_bitmask)
{
    const std::uint64_t differing_bitmask = initial_state_bitmask ^ final_state_bitmask;
    // __builtin_popcountll(differing_bitmask) == 2 ensures that the states differ by exactly 2 bits
    // and differing_bitmask & differing_bitmask >> 1ULL ensures that these differing bits are adjacent.
    return __builtin_popcountll(differing_bitmask) == 2 && differing_bitmask & differing_bitmask >> 1ULL;
}

/// @brief Reconstructs a full system state from subregion configurations.
/// @details Maps bits from a subregion and its complement into their correct positions
///          in the full system Hilbert space using Parallel Bit Deposit (PDEP).
/// @param subregion_indices A vector of indices belonging to the primary subregion.
/// @param num_sites The total size of the system.
/// @param subregion_state_bitmask The compressed state of the primary subregion.
/// @param complement_subregion_state_bitmask The compressed state of the remaining sites.
/// @return The full 64-bit state mask combining both regions.
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
    // Use hardware acceleration if available
    return _pdep_u64(subregion_state_bitmask, subregion_mask) |
           _pdep_u64(complement_subregion_state_bitmask, complement_subregion_mask);
#else
    return pdep_implementation(subregion_state_bitmask, subregion_mask) |
           pdep_implementation(complement_subregion_state_bitmask, complement_subregion_mask);
#endif
}

/// @brief Counts how many particles were created (bits flipped 0 -> 1).
/// @param initial_state_bitmask The starting state.
/// @param final_state_bitmask The ending state.
/// @return The number of sites that were empty in initial but full in final.
std::uint32_t num_sites_created(const std::uint64_t initial_state_bitmask, const std::uint64_t final_state_bitmask)
{
    return __builtin_popcountll(~initial_state_bitmask & final_state_bitmask);
}

/// @brief Counts how many particles were annihilated (bits flipped 1 -> 0).
/// @param initial_state_bitmask The starting state.
/// @param final_state_bitmask The ending state.
/// @return The number of sites that were full in initial but empty in final.
std::uint32_t num_sites_annihilated(const std::uint64_t initial_state_bitmask, const std::uint64_t final_state_bitmask)
{
    return __builtin_popcountll(initial_state_bitmask & ~final_state_bitmask);
}