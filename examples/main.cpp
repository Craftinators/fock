#include <iomanip>
#include <iostream>

#include "fock/single_particle.h"

int main()
{
    std::cout << std::setprecision(std::numeric_limits<double>::digits10) << std::fixed;
    std::cout << fock::entanglement_entropy(fock::build_subregion_correlation_matrix(6, 3, 3)) << std::endl;
}
