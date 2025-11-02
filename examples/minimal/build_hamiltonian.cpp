#include <iostream>

#include "fock/library.h"

int main()
{
    std::cout << "Full Hilbert Space:\n";
    for (const auto state_bitmask: generate_full_hilbert_space(4))
    {
        std::cout << std::bitset<4>(state_bitmask);
        if (__builtin_popcount(state_bitmask) == 2)
            std::cout << " (IN SUBSPACE)";
        std::cout << "\n";
    }
    std::cout << "Hilbert Subspace:\n";
    for (const auto state_bitmask: generate_hilbert_subspace(4, 2))
    {
        std::cout << std::bitset<4>(state_bitmask) << "\n";
    }
}