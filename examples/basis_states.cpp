#include <bitset>
#include <iostream>
#include <fock/library.hpp>

#include <fock/lattice.hpp>

int main()
{
    constexpr int lattice_size = 4, fermion_count = 2;

    for (const auto states = generate_basis_states(lattice_size, fermion_count);
        const auto state: states)
    {
        std::cout << std::bitset<lattice_size>(state) << std::endl;

        unsigned i = 0;
        for (const auto adjacent_states = get_adjacent_states(state, lattice_size, true);
            const auto adjacent_state: adjacent_states)
        {
            std::cout << "Adjacent State " << ++i << ": " << std::bitset<lattice_size>(adjacent_state) << std::endl;
        }
    }

    const Lattice lattice({ lattice_size, lattice_size },
        { BoundaryCondition::OPEN, BoundaryCondition::OPEN });
    std::cout << lattice.site_count_ << std::endl;
}
