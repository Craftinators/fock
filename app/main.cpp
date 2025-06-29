#include <iomanip>
#include <iostream>

#include "fock/library.h"

int main()
{
    const auto srcm =
        fock::build_subregion_correlation_matrix<double>(6, 3, 3);
    std::cout << std::setprecision(16) << fock::build_modular_hamiltonian<double>(srcm) << std::endl;
}
