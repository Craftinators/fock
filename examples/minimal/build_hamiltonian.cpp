#include <iostream>

#include "fock/library.h"

int main()
{
    const auto hamiltonian = build_hamiltonian(4, 2);

    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic> > solver(hamiltonian);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Eigen decomposition failed\n";
        return 1;
    }

    const Eigen::Matrix<boost::multiprecision::cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic> &eigenstates = solver.eigenvectors();

    std::cout << std::setprecision(std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits10)
            << eigenstates << std::endl;
}
