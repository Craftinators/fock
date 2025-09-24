#include <fock/library.h>
#include <matplot/matplot.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cstring>

// Convert any Eigen dense matrix to vector<vector<T>> efficiently (row-major copy)
template<typename Derived>
std::vector<std::vector<typename Derived::Scalar> >
eigenDenseToVectorRowMajor(const Eigen::MatrixBase<Derived> &matrix) {
    using T = Derived::Scalar;
    const Eigen::Index rows = matrix.rows();
    const Eigen::Index cols = matrix.cols();

    // Make a row-major copy so each row is contiguous
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_matrix = matrix;

    std::vector<std::vector<T> > out;
    out.resize(static_cast<size_t>(rows));
    for (Eigen::Index i = 0; i < rows; ++i) {
        out[i].resize(static_cast<size_t>(cols));
        // fast copy from contiguous row storage
        std::memcpy(out[i].data(),
                    row_matrix.data() + i * cols,
                    static_cast<size_t>(cols) * sizeof(T));
    }
    return out;
}

int main() {
    // build_hamiltonian returns Eigen::SparseMatrix<double>
    const Eigen::SparseMatrix<double> hamiltonian = build_hamiltonian(6, 3);

    const auto dense_hamiltonian = Eigen::MatrixXd(hamiltonian);
    const auto image_data = eigenDenseToVectorRowMajor(dense_hamiltonian);

    using namespace matplot;
    image(image_data, true);
    colorbar();
    show();
    return 0;
}