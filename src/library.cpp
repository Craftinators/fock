#include <cstdint>
#include <vector>
#include <fock/library.h>

constexpr std::uint64_t binomial_coefficient(const unsigned n, unsigned k) noexcept
{
    if (k > n - k) k = n - k;

    __uint128_t result = 1;
    for (unsigned i = 0; i < k; ++i)
    {
        result = result * (n - i) / (i + 1);
    }
    return static_cast<uint64_t>(result);
}

std::vector<std::uint64_t> generate_basis_states(const unsigned lattice_size, const unsigned fermion_count)
{
    // Validation
    if (fermion_count > lattice_size) return {}; // Can't have more fermions than sites
    if (fermion_count == 0) return { UINT64_C(0) }; // Just the |000...000> state
    if (fermion_count == lattice_size) return { (UINT64_C(1) << lattice_size) - 1 }; // Just the |111...111> state

    // Pre-calculate size using binomial coefficient
    const size_t size = binomial_coefficient(lattice_size, fermion_count);
    std::vector<std::uint64_t> result;
    result.reserve(size);

    std::uint64_t mask = (UINT64_C(1) << fermion_count) - 1;
    const std::uint64_t limit = UINT64_C(1) << lattice_size;

    while (mask < limit)
    {
        result.push_back(mask);

        // Gosper's Hack https://iamkate.com/code/hakmem-item-175/
        const std::uint64_t c = mask & -mask;
        const std::uint64_t r = mask + c;
        mask = ((r ^ mask) >> 2) / c | r;
    }

    return result;
}