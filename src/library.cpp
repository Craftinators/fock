#include <cstdint>
#include <vector>
#include <fock/library.h>
#include <Eigen/Sparse>

constexpr uint64_t binomial_coefficient(const unsigned n, unsigned k) noexcept
{
    if (k > n - k) k = n - k;

    // TODO
    // For n = 64, k = 32, result will exceed 2^64 causing overflow issues.
    // Handle arbitrary precision later (Look at mini-gmp?)
    uint64_t result = 1;
    for (unsigned i = 0; i < k; ++i)
    {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

std::vector<uint64_t> generate_basis_states(const unsigned lattice_size, const unsigned fermion_count)
{
    // Validation
    if (fermion_count > lattice_size) return {}; // Can't have more fermions than sites
    if (fermion_count == 0) return { UINT64_C(0) }; // Just the |000...000> state
    if (fermion_count == lattice_size) return { (UINT64_C(1) << lattice_size) - 1 }; // Just the |111...111> state

    // Pre-calculate size using binomial coefficient
    const size_t size = binomial_coefficient(lattice_size, fermion_count);
    std::vector<uint64_t> basis_states;
    basis_states.reserve(size);

    uint64_t mask = (UINT64_C(1) << fermion_count) - 1;
    const uint64_t limit = UINT64_C(1) << lattice_size;

    while (mask < limit)
    {
        basis_states.push_back(mask);

        // Gosper's Hack https://iamkate.com/code/hakmem-item-175/
        const uint64_t c = mask & -mask;
        const uint64_t r = mask + c;
        mask = ((r ^ mask) >> 2) / c | r;
    }

    return basis_states;
}

std::vector<uint64_t> get_adjacent_states(const uint64_t state, const unsigned lattice_size, const bool use_periodic_boundaries)
{
    std::vector<uint64_t> adjacent_states;
    if (lattice_size == 0) return adjacent_states;

    for (unsigned i = 0; i < lattice_size; ++i)
    {
        if (state >> i & 1) // Check if there is a fermion at the current site 'i'
        {
            if (use_periodic_boundaries || i < lattice_size - 1)
            {
                // With PBC, (i+1) % lattice_size wraps around.
                // With OBC, the check `i < lattice_size - 1` ensures this is just `i + 1`.

                // If the right neighbor site is empty...
                if (const unsigned right_neighbor_idx = (i + 1) % lattice_size; !(state >> right_neighbor_idx & 1))
                {
                    // Move the fermion from 'i' to the right neighbor.
                    uint64_t new_state = state;
                    new_state &= ~(UINT64_C(1) << i);                  // Clear the bit at position 'i'
                    new_state |= UINT64_C(1) << right_neighbor_idx;  // Set the bit at the neighbor's position
                    adjacent_states.push_back(new_state);
                }
            }

            if (use_periodic_boundaries || i > 0)
            {
                // With PBC, this handles the wrap-around from site 0 to the end.
                // With OBC, the check `i > 0` prevents this case, making it just `i - 1`.

                // If the left neighbor site is empty...
                if (const unsigned left_neighbor_idx = i == 0 ? lattice_size - 1 : i - 1; !(state >> left_neighbor_idx & 1))
                {
                     // Move the fermion from 'i' to the left neighbor.
                    uint64_t new_state = state;
                    new_state &= ~(UINT64_C(1) << i);              // Clear the bit at position 'i'
                    new_state |= UINT64_C(1) << left_neighbor_idx; // Set the bit at the neighbor's position
                    adjacent_states.push_back(new_state);
                }
            }
        }
    }

    return adjacent_states;
}