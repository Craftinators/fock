#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <unordered_map>

#include <boost/multiprecision/cpp_dec_float.hpp> // TODO This is just for the app, !! temporary !!

namespace fock
{
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix(const unsigned lattice_size, const unsigned fermion_count, const unsigned subregion_size)
    {
        const unsigned mode_count = 2 * lattice_size + 2;

        const T alpha_coefficient = boost::math::constants::pi<T>() / T(mode_count);
        const T beta_coefficient = alpha_coefficient * T(2 * fermion_count + 1);

        // Precompute sine functions into lookup tables
        // TODO Apparently, a vector is faster than an unordered map since the keys are consecutive integers
        std::unordered_map<std::size_t, T> alpha_sine_lookup, beta_sine_lookup;
        for (std::size_t i = 1; i < mode_count - 1; ++i)
        {
            alpha_sine_lookup[i] = sin(alpha_coefficient * T(i));
            beta_sine_lookup[i] = sin(beta_coefficient * T(i));
        }

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> correlation_matrix =
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(lattice_size, lattice_size);

        const T normalization_factor = T(1) / T(mode_count); // TODO Check if 1 needs to be wrapped in T(...)

        for (std::size_t i = 0; i < lattice_size; ++i)
        {
            for (std::size_t j = i; j < lattice_size; ++j) // The correlation matrix is symmetric
            {
                const unsigned difference = j - i;
                const unsigned sum  = i + j + 2;

                T matrix_element = T(0);
                if (difference == 0) // Handle singularity when i == j
                {
                    matrix_element = beta_coefficient / alpha_coefficient - beta_sine_lookup[sum] / alpha_sine_lookup[sum];
                } else
                {
                    matrix_element = beta_sine_lookup[difference] / alpha_sine_lookup[difference] - beta_sine_lookup[sum] / alpha_sine_lookup[sum];
                }
                correlation_matrix(i, j) = correlation_matrix(j, i) = normalization_factor * matrix_element;
            }
        }

        return correlation_matrix.topLeftCorner(subregion_size, subregion_size);
    }

    template<typename T>
    T entanglement_entropy(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subregion_correlation_matrix)
    {
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solver(subregion_correlation_matrix);
        auto eigenvalues = solver.eigenvalues();

        T entropy = T(0);

        for (std::size_t i = 0; i < subregion_correlation_matrix.cols(); ++i)
        {
            const T eigenvalue = eigenvalues(i);
            const T epsilon = std::numeric_limits<T>::epsilon();

            if (eigenvalue < epsilon || eigenvalue > T(1) - epsilon) continue;
            entropy -= eigenvalue * log(eigenvalue) + (T(1) - eigenvalue) * log(T(1) - eigenvalue);
        }

        return entropy;
    }

    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subregion_correlation_matrix)
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
}

#endif //FOCK_LIBRARY_H
