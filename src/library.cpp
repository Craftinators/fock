#include <cstdint>
#include <vector>
#include <fock/library.h>

constexpr std::uint64_t binomial_coefficient(const unsigned int n, unsigned int k) noexcept
{
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;

    if (k > n - k) k = n - k;

    std::uint64_t result = 1;
    for (int i = 0; i < k; ++i)
    {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

std::vector<std::uint64_t> generate_basis_states(const unsigned int lattice_size, const unsigned int fermion_count)
{
    if (fermion_count <= 0 || fermion_count > lattice_size) return {};

    // Pre-calculate size using binomial coefficient
    const size_t size = binomial_coefficient(lattice_size, fermion_count);
    std::vector<std::uint64_t> result;
    result.reserve(size);

    std::uint64_t mask = (static_cast<std::uint64_t>(1) << fermion_count) - 1;
    const std::uint64_t limit = static_cast<std::uint64_t>(1) << lattice_size;

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