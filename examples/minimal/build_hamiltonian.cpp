#include <iostream>

#include "fock/library.h"

int main()
{
    const auto hamiltonian = build_hamiltonian(6, 3);
    std::cout << hamiltonian << std::endl;
}
