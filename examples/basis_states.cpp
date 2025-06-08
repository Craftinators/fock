#include <bitset>
#include <iostream>
#include <fock/library.h>

int main()
{
    constexpr int n = 7, k = 4;

    for (const auto states = generate_basis_states(n, k);
        const auto state: states)
    {
        std::cout << std::bitset<n>(state) << std::endl;
    }
}