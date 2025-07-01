#include "fock/library.h"

namespace fock
{
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    fock::build_subregion_correlation_matrix(const unsigned lattice_size,
        const unsigned fermion_count,
        const unsigned subregion_size)
    {
        // TODO Replace with exceptions, this just temporary
        assert(lattice_size != 0);
        assert(fermion_count <= lattice_size);
        assert(subregion_size <= lattice_size);

        // Don't bother if subregion is 0, that's just a 0 by 0 matrix
        // TODO Should I just error here? I mean realistically what are you doing with a 0 by 0 matrix lol
        if (subregion_size == 0)
        {
            return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(0, 0);
        }

        // TODO Replace with perf code later, I need to take a closer look at it
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> subregion_correlation_matrix =
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(lattice_size, lattice_size);

        T normalization_factor = T(1) / T(2 * (lattice_size + 1));
        T lower_sine_coefficient = boost::math::constants::pi<T>() * normalization_factor;
        T upper_sine_coefficient = lower_sine_coefficient * T(2 * fermion_count + 1);

        // Precompute Dirichlet kernels
        std::unordered_map<unsigned, T> difference_kernel_values;
        // When difference == 0, it results in singularity, thus we do not need to precompute this value.
        // The largest value occurs when j == subregion_size - 1 (maximum) and i == 0 (minimum)
        for (unsigned i = 1; i < subregion_size; ++i)
        {
            difference_kernel_values[i] = sin(upper_sine_coefficient * i) / sin(lower_sine_coefficient * i);
        }

        std::unordered_map<unsigned, T> sum_kernel_values;
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

                T matrix_element = T(0);

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

    template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix<double>(unsigned, unsigned, unsigned);
    template Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix<boost::multiprecision::cpp_dec_float_100>(unsigned, unsigned, unsigned);
}