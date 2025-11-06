#include <iostream>

#include "fock/library.h"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

int main()
{
    const auto hamiltonian = build_sparse_hamiltonian_sparse(12, 6);

    Spectra::SparseGenMatProd<double> operation(hamiltonian);
    Spectra::GenEigsSolver eigen_solver(operation, 15, 20);

    eigen_solver.init();
    eigen_solver.compute(Spectra::SortRule::SmallestReal);

    Eigen::VectorXcd eigenvalues;
    Eigen::MatrixXcd eigenvectors;
    if(eigen_solver.info() == Spectra::CompInfo::Successful)
    {
        eigenvalues = eigen_solver.eigenvalues();
        eigenvectors = eigen_solver.eigenvectors();
    }

    std::cout << "\n13th eigenvalue: " << eigenvalues(12) << "\n13th eigenvector:\n" << eigenvectors.col(12) <<
        "\n\n14th eigenvalue: " << eigenvalues(13) << "\n14th eigenvector:\n" << eigenvectors.col(13) << std::endl;

    return 0;
}
