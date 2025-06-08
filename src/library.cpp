#include <cstdint>
#include <vector>
#include <fock/library.h>

constexpr size_t binomial_coefficient(int n, int k)
{
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;

    if (k > n - k) k = n - k;

    size_t result = 1;
    for (int i = 0; i < k; ++i)
    {
        result = result * (n - i ) / (i + 1);
    }
    return result;
}

std::vector<std::uint64_t> generate_basis_states(int n, int k)
{
    if (k <= 0 || k > n) return {};

    // Pre-calculate size using binomial coefficient
    const size_t size = binomial_coefficient(n, k);
    std::vector<std::uint64_t> result;
    result.reserve(size);

    std::uint64_t mask = (static_cast<std::uint64_t>(1) << k) - 1;
    const std::uint64_t limit = static_cast<std::uint64_t>(1) << n;

    while (mask < limit)
    {
        result.push_back(mask);

        const std::uint64_t c = mask & -mask;
        const std::uint64_t r = mask + c;
        mask = ((r ^ mask) >> 2) / c | r;
    }

    return result;
}