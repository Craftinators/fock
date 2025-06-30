#include "fock/library.h"

namespace fock
{
    template<typename T>
    T entanglement_entropy(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &subregion_correlation_matrix)
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

    template double entanglement_entropy<double>(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&);
    template boost::multiprecision::cpp_dec_float_100 entanglement_entropy<boost::multiprecision::cpp_dec_float_100>(const Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic>&);
}