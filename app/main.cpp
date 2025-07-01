#include <iomanip>
#include <iostream>

#include "fock/library.h"

using precision_type = double;

int main()
{
    std::cout << std::setprecision(std::numeric_limits<precision_type>::digits10) << std::fixed;
    const auto subregion_correlation_matrix = fock::build_subregion_correlation_matrix<precision_type>(6, 3, 3);
    std::cout << fock::entanglement_entropy<precision_type>(subregion_correlation_matrix) << std::endl;
}
