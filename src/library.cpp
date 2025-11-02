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

std::vector<unsigned long long> generate_hilbert_subspace(const unsigned int num_sites, const unsigned int num_filled_sites,
                                                          unsigned long long &num_states)
{
    if (num_filled_sites > num_sites)
    {
        throw std::logic_error("Number of filled sites cannot be greater than the number of sites.");
    }

    num_states = get_binomial_coefficient(num_sites, num_filled_sites);
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

bool is_valid_region_composition(const unsigned int num_filled_sites,
                                 const unsigned long long subregion_state_bitmask, const unsigned long long complement_subregion_state_bitmask)
{
    return __builtin_popcountll(subregion_state_bitmask) + __builtin_popcountll(complement_subregion_state_bitmask) == num_filled_sites;
}

bool are_neighbors(const unsigned long long initial_state_bitmask, const unsigned long long final_state_bitmask)
{
    const unsigned long long differing_bitmask = initial_state_bitmask ^ final_state_bitmask;
    // __builtin_popcountll(differing_bitmask) == 2 ensures that the states differ by exactly 2 bits
    // and differing_bitmask & differing_bitmask >> 1ULL ensures that these differing bits are adjacent.
    return __builtin_popcountll(differing_bitmask) == 2 && differing_bitmask & differing_bitmask >> 1ULL;
}

// TODO separate this from code above

#define EIGEN_DONT_VECTORIZE
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

const boost::multiprecision::cpp_dec_float_50 zero("0");
const boost::multiprecision::cpp_dec_float_50 minus_one("-1");

Eigen::Matrix<boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic>
build_hamiltonian(const unsigned int num_sites, const unsigned int num_filled_sites)
{
    unsigned long long num_states;
    const auto states = generate_hilbert_subspace(num_sites, num_filled_sites, num_states);

    Eigen::Matrix<boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic>
            hamiltonian(num_states, num_states);
    hamiltonian.setZero();

    for (Eigen::Index i = 0; i < num_states; ++i)
    {
        for (Eigen::Index j = i + 1; j < num_states; ++j)
        {
            if (are_neighbors(states[i], states[j]))
            {
                hamiltonian(i, j) = hamiltonian(j, i) = minus_one;
            }
        }
    }

    return hamiltonian;
}
