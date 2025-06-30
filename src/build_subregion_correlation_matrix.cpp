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
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> correlation_matrix =
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(lattice_size, lattice_size);

        T pi = boost::math::constants::pi<T>();
        T normalization_factor = T(2) / T(lattice_size + 1);
        T sine_value = pi / T(lattice_size + 1);
        for (int i = 0; i < lattice_size; ++i)
        {
            for (int j = 0; j < lattice_size; ++j)
            {
                T matrix_element = T(0);
                for (int k = 1; k < fermion_count + 1; ++k)
                {
                    matrix_element += sin(sine_value * T(k) * T(i + 1)) * sin(sine_value * T(k) * T(j + 1));
                }
                correlation_matrix(i, j) = normalization_factor * matrix_element;
            }
        }

        return correlation_matrix.topLeftCorner(subregion_size, subregion_size);
    }

    template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix<double>(unsigned, unsigned, unsigned);
    template Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix<boost::multiprecision::cpp_dec_float_100>(unsigned, unsigned, unsigned);
}