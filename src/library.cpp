#include <vector>
#include <numeric>
#include <stdexcept>

inline unsigned long long get_next_bit_permutation(const unsigned long long current_bitmask)
{
    const unsigned long long least_significant_bit = current_bitmask & -current_bitmask;
    const unsigned long long next_rippled_bitmask = current_bitmask + least_significant_bit;
    return ((next_rippled_bitmask ^ current_bitmask) >> 2) / least_significant_bit | next_rippled_bitmask;
}

unsigned long long get_binomial_coefficient(const unsigned int n, unsigned int k)
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
        return 1;
    }

    unsigned long long result = 1;
    for (unsigned int i = 0; i < k; ++i)
    {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

std::vector<unsigned long long> generate_full_hilbert_space(const unsigned int num_sites)
{
    if (num_sites >= 64)
    {
        throw std::runtime_error("Number of sites must be less than 64 to fit in unsigned long long and std::vector.");
    }

    const std::size_t num_states = 1ULL << num_sites;
    std::vector<unsigned long long> states(num_states);
    std::iota(states.begin(), states.end(), 0ULL);
    return states;
}

std::vector<unsigned long long> generate_hilbert_subspace(const unsigned int num_sites, const unsigned int num_filled_sites)
{
    if (num_filled_sites > num_sites)
    {
        throw std::logic_error("Number of filled sites cannot be greater than the number of sites.");
    }

    const std::size_t num_states = get_binomial_coefficient(num_sites, num_filled_sites);
    std::vector<unsigned long long> states(num_states);

    if (num_filled_sites == 0)
    {
        if (num_states == 1) states[0] = 0ULL;
        return states;
    }

    unsigned long long current_state_bitmask = (1ULL << num_filled_sites) - 1ULL;
    for (std::size_t state_index = 0; state_index < num_states; ++state_index)
    {
        states[state_index] = current_state_bitmask;

        if (state_index < num_states - 1)
        {
            current_state_bitmask = get_next_bit_permutation(current_state_bitmask);
        }
    }
    return states;
}