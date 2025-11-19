#include <iostream>

#include "fock/library.h"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

// --- Include your header that contains all functions here ---
// #include "your_header.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>

#include "../../cmake-build-debug/_deps/matplotplusplus-src/source/matplot/matplot.h"

// Sweep eps, compute K = -log(rhoA), dump diagnostics to CSV.
// - psi_n, psi_m: Eigen::VectorXd (same ordering as generate_hilbert_subspace / build_reduced_density_matrix expects).
// - subregion: vector<uint32_t> as before.
// - num_sites, num_filled_sites: ints
// - out_filename: where CSV is written
// - eps_start, eps_end, eps_step: sweep parameters (eps values are positive magnitudes)
void sweep_epsilon_and_dump(
    const Eigen::VectorXd &psi_n,
    const Eigen::VectorXd &psi_m,
    const std::vector<uint32_t> &subregion,
    const uint32_t num_sites,
    const uint32_t num_filled_sites,
    const std::string &out_filename,
    const double eps_start = 1e-6,
    const double eps_end = 1e-1,
    const double eps_step = 1e-2)
{
    // Basic checks
    if (psi_n.size() != psi_m.size())
        throw std::invalid_argument("psi_n and psi_m must have same size.");
    if (eps_end <= eps_start) throw std::invalid_argument("eps_end must be > eps_start.");
    if (eps_step <= 0.0) throw std::invalid_argument("eps_step must be > 0.");

    // Build K0 (epsilon = 0)
    {
        // ensure psi_n normalized? Eigen solver usually returns normalized eigenvectors, but we normalize anyway.
    }
    Eigen::VectorXd psi0 = psi_n;
    psi0 /= psi0.norm();

    // Reduced density at epsilon = 0
    Eigen::SparseMatrix<double> rhoA0_sparse = build_reduced_density_matrix(psi0, subregion, num_sites, num_filled_sites);
    Eigen::MatrixXd rhoA0 = Eigen::MatrixXd(rhoA0_sparse);
    // Symmetrize to counter numeric asymmetry
    rhoA0 = 0.5 * (rhoA0 + rhoA0.transpose());

    std::cout << rhoA0_sparse.nonZeros() << std::endl;
    std::cout << rhoA0_sparse.cols() * rhoA0_sparse.rows() << std::endl;

    // Diagonalize rhoA0, regularize small eigenvalues, form K0 = -U log(lam) U^T
    const double eps_reg = 1e-12; // regularization floor for eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> r0_solver(rhoA0);
    if (r0_solver.info() != Eigen::Success) throw std::runtime_error("Diagonalization of rhoA0 failed");

    Eigen::VectorXd lam0 = r0_solver.eigenvalues();
    Eigen::MatrixXd U0 = r0_solver.eigenvectors();
    for (int i = 0; i < lam0.size(); ++i) if (lam0(i) < eps_reg) lam0(i) = eps_reg;
    Eigen::VectorXd loglam0 = lam0.array().log();
    Eigen::MatrixXd K0 = -U0 * loglam0.asDiagonal() * U0.transpose();

    // We'll write a CSV with a header
    std::ofstream csv(out_filename);
    if (!csv) throw std::runtime_error("Failed to open output CSV: " + out_filename);
    csv << std::setprecision(12) << std::scientific;
    csv << "epsilon,fro_norm,ratio_lin,ratio_quad,op_norm,K00,K01,traceK,trace_rhoA\n";

    const double tiny = 1e-300; // guard against division by zero when eps very small

    // Sweep
    for (double eps = eps_start; eps <= eps_end + 1e-12; eps += eps_step)
    {
        // Build normalized superposition (assume real eps for now; change to complex if needed)
        const double norm_factor = std::sqrt(1.0 + eps * eps);
        Eigen::VectorXd psi_eps = (psi_n + eps * psi_m) / norm_factor;

        // Build reduced density matrix
        Eigen::SparseMatrix<double> rhoA_sparse = build_reduced_density_matrix(psi_eps, subregion, num_sites, num_filled_sites);
        Eigen::MatrixXd rhoA = Eigen::MatrixXd(rhoA_sparse);
        rhoA = 0.5 * (rhoA + rhoA.transpose()); // symmetrize

        // Diagonalize rhoA
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> rsolver(rhoA);
        if (rsolver.info() != Eigen::Success)
        {
            std::cerr << "Warning: diagonalization of rhoA failed at eps=" << eps << " — skipping\n";
            continue;
        }
        Eigen::VectorXd lam = rsolver.eigenvalues();
        Eigen::MatrixXd U = rsolver.eigenvectors();

        // Regularize small eigenvalues before log
        for (int i = 0; i < lam.size(); ++i) if (lam(i) < eps_reg) lam(i) = eps_reg;
        Eigen::VectorXd loglam = lam.array().log();
        Eigen::MatrixXd K = -U * loglam.asDiagonal() * U.transpose();

        // Delta K = K - K0 (additive constants cancel)
        // (Note: K0 and K are in the same reduced basis ordering returned by build_reduced_density_matrix)
        // If K and K0 have different sizes, something is wrong; assert same size:
        if (K.rows() != K0.rows() || K.cols() != K0.cols())
            throw std::runtime_error("Dimension mismatch between K and K0.");

        Eigen::MatrixXd dK = K - K0;

        // Frobenius norm
        double fro_norm = std::sqrt((dK.array() * dK.array()).sum());

        // operator norm (for symmetric dK, largest absolute eigenvalue)
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> dKsolver(dK);
        double op_norm = 0.0;
        if (dKsolver.info() == Eigen::Success)
        {
            op_norm = dKsolver.eigenvalues().cwiseAbs().maxCoeff();
        } else
        {
            // fallback: use matrix norm via singular values (more expensive); here we just set NaN
            op_norm = std::numeric_limits<double>::quiet_NaN();
        }

        // selected elements (check size)
        double K00 = 0.0;
        double K01 = 0.0;
        if (K.rows() > 0 && K.cols() > 0) K00 = K(0, 0);
        if (K.rows() > 1 && K.cols() > 1) K01 = K(0, 1);

        double traceK = K.trace();
        double traceRhoA = rhoA.trace(); // should be 1.0

        // Guard against dividing by zero when eps extremely small
        double ratio_lin = (eps > tiny) ? (fro_norm / eps) : std::numeric_limits<double>::infinity();
        double ratio_quad = (eps > tiny) ? (fro_norm / (eps * eps)) : std::numeric_limits<double>::infinity();

        // write CSV line
        csv << eps << "," << fro_norm << "," << ratio_lin << "," << ratio_quad << ","
                << op_norm << "," << K00 << "," << K01 << "," << traceK << "," << traceRhoA << "\n";
    }

    csv.close();
    std::cout << "Wrote sweep results to " << out_filename << "\n";
}


int main()
{
    const uint32_t num_sites = 12;
    const uint32_t num_filled_sites = 6;

    Eigen::SparseMatrix<double> hamiltonian = build_hamiltonian(num_sites, num_filled_sites);
    auto dense_hamiltonian = Eigen::MatrixXd(hamiltonian);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(dense_hamiltonian);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Diagonalization failed!\n";
        return 1;
    }

    // pick two excited states you want (indices must be valid)
    Eigen::VectorXd excited_state_a = solver.eigenvectors().col(2);
    Eigen::VectorXd excited_state_b = solver.eigenvectors().col(3);
    std::vector<uint32_t> subregion = {0, 1};

    sweep_epsilon_and_dump(excited_state_a, excited_state_b, subregion,
                           num_sites, num_filled_sites,
                           "eps_sweep_results.csv",
                           1e-10, 1e-5, 1e-10); // sweep from 1e-6 to 1e-1 step 1e-3

    // --- Read CSV for plotting ---
    std::ifstream file("eps_sweep_results.csv");
    std::string line;
    std::vector<double> eps_vals, fro_vals;

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;

        while (std::getline(ss, field, ','))
        {
            fields.push_back(field);
        }

        if (fields.size() < 2) continue;

        eps_vals.push_back(std::stod(fields[0]));
        fro_vals.push_back(std::stod(fields[1]));
    }

    // --- Plot using Matplot++ ---
    using namespace matplot;

    auto fig = figure(true);
    plot(eps_vals, fro_vals);
    xlabel("epsilon");
    ylabel("fro_norm");
    title("Frobenius Norm vs Epsilon");
    grid(true);

    show();

    return 0;
}

