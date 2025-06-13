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
    Lattice(std::initializer_list<std::size_t> dimension_sizes,
        std::initializer_list<BoundaryCondition> boundary_conditions);

    // Copy, equality, and destruction
    Lattice(const Lattice&)            = default;
    Lattice(Lattice&&)                 = default;
    Lattice& operator=(const Lattice&) = default;
    Lattice& operator=(Lattice&&)      = default;
    ~Lattice()                         = default;

    // Accessors
    [[nodiscard]]
    std::size_t dimension(const std::size_t axis) const noexcept { return dimension_sizes_[axis]; }
    [[nodiscard]]
    BoundaryCondition boundary_condition(const std::size_t axis) const noexcept { return boundary_conditions_[axis]; }
    [[nodiscard]]
    std::size_t rank() const noexcept { return dimension_sizes_.size(); }
    [[nodiscard]]
    std::size_t site_count() const noexcept { return site_count_; }

private:
    std::vector<std::size_t>          dimension_sizes_;
    std::vector<BoundaryCondition>    boundary_conditions_;
    std::size_t                       site_count_;

    static std::size_t safe_prod(const std::vector<std::size_t>& dimension_sizes);
};

#endif //LATTICE_HPP
