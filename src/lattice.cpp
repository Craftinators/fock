#include <limits>
#include <sstream>
#include <stdexcept>
#include <fock/lattice.hpp>

Lattice::Lattice(const std::initializer_list<std::size_t> dimension_sizes, const std::initializer_list<BoundaryCondition> boundary_conditions)
    : dimension_sizes_{ dimension_sizes }, boundary_conditions_{ boundary_conditions }, site_count_{ safe_prod(dimension_sizes_) }
{
    if (dimension_sizes_.size() != boundary_conditions_.size())
    {
        std::ostringstream oss;
        oss << "Dimension count (" << dimension_sizes_.size()
            << ") does not match BC count ("
            << boundary_conditions_.size() << ").";
        throw std::invalid_argument{oss.str()};
    }
}

std::size_t Lattice::safe_prod(const std::vector<std::size_t> &dimension_sizes)
{
    if (dimension_sizes.empty())
    {
        throw std::invalid_argument{"Lattice dimensions must be non-empty."};
    }

    std::size_t accumulation = 1;
    for (std::size_t i = 0; i < dimension_sizes.size(); ++i)
    {
        std::size_t dimension_size = dimension_sizes[i];
        if (dimension_size == 0)
        {
            std::ostringstream oss;
            oss << "Dimension at axis " << i << " must be > 0; got 0.";
            throw std::invalid_argument{oss.str()};
        }

        if (dimension_size > std::numeric_limits<std::size_t>::max() / accumulation)
        {
            throw std::overflow_error{"Product of dimensions overflows size_t."};
        }

        accumulation *= dimension_size;
    }

    return accumulation;
}
