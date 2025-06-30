  #include <iomanip>
#include <iostream>

#include "fock/library.h"

int main()
{
    using type = double;

    std::cout << std::setprecision(std::numeric_limits<type>::digits10) << std::fixed;
    const auto subregion_correlation_matrix = fock::build_subregion_correlation_matrix<type>(6, 3, 3);
    std::cout << fock::entanglement_entropy<type>(subregion_correlation_matrix) << std::endl;
}
