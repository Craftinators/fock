#ifndef FOCK_LIBRARY_H
#define FOCK_LIBRARY_H

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace fock
{
    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix(unsigned lattice_size, unsigned fermion_count, unsigned subregion_size);

    extern template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix<double>(unsigned, unsigned, unsigned);
    extern template Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic> build_subregion_correlation_matrix<boost::multiprecision::cpp_dec_float_100>(unsigned, unsigned, unsigned);

    template<typename T>
    T entanglement_entropy(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& subregion_correlation_matrix);

    extern template double entanglement_entropy(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& subregion_correlation_matrix);
    extern template boost::multiprecision::cpp_dec_float_100 entanglement_entropy(const Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic>& subregion_correlation_matrix);

    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& subregion_correlation_matrix);

    extern template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& subregion_correlation_matrix);
    extern template Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic> build_modular_hamiltonian(const Eigen::Matrix<boost::multiprecision::cpp_dec_float_100, Eigen::Dynamic, Eigen::Dynamic>& subregion_correlation_matrix);
}

#endif //FOCK_LIBRARY_H
