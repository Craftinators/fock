#include "fock/single_particle.h"
#include <numbers>

namespace fock
{
    Eigen::MatrixXd build_subregion_correlation_matrix(
        const unsigned lattice_length,
        const unsigned fermion_count,
        const unsigned subregion_size)
    {
        // Don't bother if subregion is 0, that's just a 0 by 0 matrix
        if (subregion_size == 0)
        {
            return Eigen::MatrixXd::Zero(0, 0);
        }

        Eigen::MatrixXd subregion_correlation_matrix = Eigen::MatrixXd::Zero(subregion_size, subregion_size);

        const double normalization_factor = 1.0 / (2 * (lattice_length + 1));
        const double lower_sine_coefficient = std::numbers::pi * normalization_factor;
        const double upper_sine_coefficient = lower_sine_coefficient * (2 * fermion_count + 1);

        // Precompute Dirichlet kernels
        std::unordered_map<unsigned, double> difference_kernel_values;
        // When difference == 0, it results in singularity, thus we do not need to precompute this value.
        // The largest value occurs when j == subregion_size - 1 (maximum) and i == 0 (minimum)
        for (unsigned i = 1; i < subregion_size; ++i)
        {
            difference_kernel_values[i] = sin(upper_sine_coefficient * i) / sin(lower_sine_coefficient * i);
        }

        std::unordered_map<unsigned, double> sum_kernel_values;
        // Minimum value occurs when i == j == 0, which results in sum == 2, whereas maximum occurs
        // when i == j == subregion_size - 1, when results in sum == 2 * subregion_size
        for (unsigned i = 2; i < 2 * subregion_size + 1; ++i)
        {
            sum_kernel_values[i] = sin(upper_sine_coefficient * i) / sin(lower_sine_coefficient * i);
        }

        for (unsigned i = 0; i < subregion_size; ++i)
        {
            for (unsigned j = i; j < subregion_size; ++j)
            {
                unsigned difference = j - i;
                unsigned sum = i + j + 2; // "+2" is an artifact of 0-indexing

                double matrix_element = 0;

                if (difference == 0) // Handle singularity
                {
                    matrix_element = upper_sine_coefficient / lower_sine_coefficient - sum_kernel_values[sum];
                }
                else
                {
                    matrix_element = difference_kernel_values[difference] - sum_kernel_values[sum];
                }

                subregion_correlation_matrix(i, j) = subregion_correlation_matrix(j, i) = normalization_factor * matrix_element;
            }
        }

        return subregion_correlation_matrix;
    }

    double entanglement_entropy(const Eigen::MatrixXd &subregion_correlation_matrix)
    {
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(subregion_correlation_matrix);
        auto eigenvalues = solver.eigenvalues();

        double entropy = 0;

        for (int i = 0; i < subregion_correlation_matrix.cols(); ++i)
        {
            const double eigenvalue = eigenvalues(i);

            if (constexpr double epsilon = std::numeric_limits<double>::epsilon();
                eigenvalue < epsilon || eigenvalue > 1 - epsilon) continue;

            entropy -= eigenvalue * log(eigenvalue) + (1 - eigenvalue) * log(1 - eigenvalue);
        }

        return entropy;
    }

    Eigen::MatrixXd build_modular_hamiltonian(const Eigen::MatrixXd &subregion_correlation_matrix)
    {
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(subregion_correlation_matrix);
        auto eigenvalues = solver.eigenvalues();
        auto eigenvectors = solver.eigenvectors();

        Eigen::VectorXd entanglement_spectrum(subregion_correlation_matrix.cols());

        for (int i = 0; i < subregion_correlation_matrix.cols(); ++i)
        {
            const auto eigenvalue = eigenvalues(i);
            entanglement_spectrum(i) = log((1 - eigenvalue) / eigenvalue);
        }

        return eigenvectors * entanglement_spectrum.asDiagonal() * eigenvectors.transpose();
    }
}