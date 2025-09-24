#include <fock/library.h> // Assuming this contains gospers_hack
#include <cinttypes>
#include <limits>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Helper functions (unchanged)
uint64_t gospers_hack(const uint64_t v) {
    const uint64_t c = v & -v;
    const uint64_t r = v + c;
    return (((r ^ v) >> 2) / c) | r;
}

uint64_t binomial_coefficient(const uint32_t n, uint32_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (n - k < k) k = n - k;
    __int128 res = 1;
    for (uint32_t i = 1; i <= k; ++i) {
        res = res * static_cast<__int128>(n - k + i) / static_cast<__int128>(i);
        if (res > static_cast<__int128>(std::numeric_limits<uint64_t>::max())) {
            throw std::overflow_error("binomial result doesn't fit in uint64_t");
        }
    }
    return static_cast<uint64_t>(res);
}

/**
 * @brief Builds the sparse Hamiltonian matrix for a 1D tight-binding model.
 *
 * This function uses an efficient generation method to construct the Hamiltonian,
 * avoiding the O(N^2) complexity of checking all state pairs. It models
 * nearest-neighbor hopping with open boundary conditions.
 *
 * @param lattice_length The number of sites on the lattice.
 * @param num_fermions The number of fermions (and electrons) in the system.
 * @return A sparse matrix representing the Hamiltonian.
 */
Eigen::SparseMatrix<double> build_hamiltonian(const unsigned int lattice_length, const unsigned int num_fermions) {
    // 1. Build the basis and a map for O(1) state-to-index lookups
    const uint64_t hilbert_size = binomial_coefficient(lattice_length, num_fermions);
    std::vector<uint64_t> basis_states;
    std::unordered_map<uint64_t, uint32_t> state_to_index;

    basis_states.reserve(hilbert_size);
    state_to_index.reserve(hilbert_size);

    uint64_t current_basis_state = (1ULL << num_fermions) - 1ULL;
    for (uint32_t i = 0; i < hilbert_size; ++i) {
        basis_states.push_back(current_basis_state);
        state_to_index[current_basis_state] = i;
        if (i < hilbert_size - 1) {
            current_basis_state = gospers_hack(current_basis_state);
        }
    }

    // 2. Build the Hamiltonian by generating connections
    // We create a list of non-zero entries (triplets) to build the sparse matrix
    std::vector<Eigen::Triplet<double> > triplets;

    for (uint32_t i = 0; i < hilbert_size; ++i) {
        const uint64_t initial_state = basis_states[i];

        // Iterate through lattice sites to find possible hops (p -> p+1)
        // This loop implies open boundary conditions, matching the original code.
        for (unsigned int p = 0; p < lattice_length - 1; ++p) {
            const uint64_t occupied_site_mask = 1ULL << p;

            // A hop is possible if a fermion at site 'p' can move to an empty site 'p+1'
            // We check this by XORing the state with the masks. If the result has fewer
            // bits set, it means one was occupied and one was empty.
            if (const uint64_t adjacent_site_mask = 1ULL << (p + 1);
                __builtin_popcountll(initial_state & (occupied_site_mask | adjacent_site_mask)) == 1) {
                const uint64_t final_state = initial_state ^ (occupied_site_mask | adjacent_site_mask);

                // Find the index 'j' of the resulting state using our map
                const uint32_t j = state_to_index.at(final_state);

                // Add non-zero elements for H(i, j) and H(j, i)
                // We only need to check for i < j to fill the matrix symmetrically
                // without redundant work.
                if (i < j) {
                    triplets.emplace_back(i, j, -1.0);
                    triplets.emplace_back(j, i, -1.0);
                }
            }
        }
    }

    // 3. Construct the sparse matrix from the triplets
    // Explicitly cast hilbert_size to Eigen::Index to resolve the warning
    const auto matrix_size = static_cast<Eigen::Index>(hilbert_size);
    Eigen::SparseMatrix<double> hamiltonian(matrix_size, matrix_size);
    hamiltonian.setFromTriplets(triplets.begin(), triplets.end());

    return hamiltonian;
}