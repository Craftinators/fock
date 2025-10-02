#include <fock/library.h>
#include <cinttypes>
#include <limits>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

inline uint64_t gospers_hack(const uint64_t v) {
    const uint64_t c = v & -v;
    const uint64_t r = v + c;
    return ((r ^ v) >> 2) / c | r;
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

Eigen::SparseMatrix<double> build_hamiltonian(const unsigned int lattice_length, const unsigned int num_fermions) {
    const uint64_t hilbert_size_u64 = binomial_coefficient(lattice_length, num_fermions);

    const auto matrix_size = static_cast<std::size_t>(hilbert_size_u64);
    const auto eigen_matrix_size = static_cast<Eigen::Index>(matrix_size);

    std::vector<uint64_t> basis_states;
    basis_states.reserve(matrix_size);

    std::unordered_map<uint64_t, uint64_t> state_to_index;
    state_to_index.reserve(matrix_size);
    state_to_index.max_load_factor(0.7f);

    uint64_t current = num_fermions == 0 ? 0ULL : (1ULL << num_fermions) - 1ULL;
    for (uint64_t idx = 0; idx < static_cast<uint64_t>(hilbert_size_u64); ++idx) {
        basis_states.push_back(current);
        state_to_index.emplace(current, idx);
        if (idx + 1 < hilbert_size_u64) current = gospers_hack(current);
    }

    const uint64_t est_edges_per_state = std::min<uint64_t>(num_fermions, lattice_length > 0 ? lattice_length - 1 : 0);
    const uint64_t estimated_pairs = hilbert_size_u64 * est_edges_per_state;
    const std::size_t reserve_sz = estimated_pairs > std::numeric_limits<std::size_t>::max() / 2
                                       ? std::numeric_limits<std::size_t>::max()
                                       : static_cast<std::size_t>(estimated_pairs * 2);
    std::vector<Eigen::Triplet<double> > triplets;
    triplets.reserve(std::min<std::size_t>(reserve_sz, 1ull << 30));

    for (uint64_t i = 0; i < static_cast<uint64_t>(matrix_size); ++i) {
        const uint64_t state = basis_states[static_cast<std::size_t>(i)];
        uint64_t x = state;
        while (x) {
            const unsigned int p = __builtin_ctzll(x);
            x &= x - 1;
            if (p + 1 >= lattice_length) continue;
            const bool bit_p = (state >> p & 1ULL) != 0ULL;
            if (const bool bit_p1 = (state >> (p + 1) & 1ULL) != 0ULL; !bit_p || bit_p1) continue;
            const uint64_t new_state = state ^ (1ULL << p | 1ULL << (p + 1));
            auto it = state_to_index.find(new_state);
            if (it == state_to_index.end()) continue;
            if (const uint64_t j = it->second; i < j) {
                triplets.emplace_back(static_cast<int>(i), static_cast<int>(j), -1.0);
                triplets.emplace_back(static_cast<int>(j), static_cast<int>(i), -1.0);
            }
        }
    }

    Eigen::SparseMatrix<double> hamiltonian(eigen_matrix_size, eigen_matrix_size);
    hamiltonian.setFromTriplets(triplets.begin(), triplets.end());
    return hamiltonian;
}