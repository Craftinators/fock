#include <bitset>
#include <iostream>
#include <fock/library.h>

int main()
{
    const int n = 7;
    int k = 4;

    auto states = generate_basis_states(n, k);

    for (auto state: states)
    {
        std::cout << std::bitset<n>(state) << std::endl;
    }
}