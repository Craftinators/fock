#include "fock/single_particle.h"

namespace fock
{
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &subregion_correlation_matrix)
    {
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solver(subregion_correlation_matrix);
        auto eigenvalues = solver.eigenvalues();
        auto eigenvectors = solver.eigenvectors();

        // TODO Do we need to be setting elements to zero, or is ... entanglement_energies(...) enough?
        Eigen::Matrix<T, Eigen::Dynamic, 1> entanglement_spectrum = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(subregion_correlation_matrix.cols());

        for (std::size_t i = 0; i < subregion_correlation_matrix.cols(); ++i)
        {
            auto eigenvalue = eigenvalues(i);
            entanglement_spectrum(i) = log((T(1) - eigenvalue) / eigenvalue);
        }

        return eigenvectors * entanglement_spectrum.asDiagonal() * eigenvectors.transpose();
    }

    template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian<double>(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&);
    template Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian<boost::multiprecision::cpp_dec_float_100>(const Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic>&);
}