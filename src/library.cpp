#include <cstdint>
#include <vector>
#include <fock/library.h>

std::vector<std::uint64_t> generate_basis_states(int n, int k)
{
    std::vector<std::uint64_t> result;
    if (k <= 0 || k > n) return result;

    // initial mask: k the lowest bits set
    std::uint64_t mask = (static_cast<std::uint64_t>(1) << k) - 1;
    const std::uint64_t limit = static_cast<std::uint64_t>(1) << n;

    while (mask < limit)
    {
        result.push_back(mask);

        // Gosper’s hack to get next mask with same popcount
        const std::uint64_t c = mask & -mask;
        const std::uint64_t r = mask + c;
        mask = ((r ^ mask) >> 2) / c | r;
    }

    return result;
}