#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <vector>

enum class BoundaryCondition
{
    OPEN,
    PERIODIC,
};

class Lattice
{
public:
    Lattice(const std::vector<unsigned>& dimensions, const std::vector<BoundaryCondition>& boundary_conditions);
public:
    std::vector<unsigned> dimensions_;
    std::vector<BoundaryCondition> boundary_conditions_;
    unsigned site_count_;
};

#endif //LATTICE_HPP
