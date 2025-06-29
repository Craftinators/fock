#include <iomanip>
#include <iostream>

#include "fock/library.h"

int main()
{
    const auto srcm =
        fock::build_subregion_correlation_matrix<boost::multiprecision::cpp_dec_float_100>(6, 3, 3);
    std::cout << std::setprecision(100) << fock::entanglement_entropy<boost::multiprecision::cpp_dec_float_100>(srcm) << std::endl;
}
