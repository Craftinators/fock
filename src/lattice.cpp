#include <sstream>
#include <stdexcept>
#include <fock/lattice.hpp>

Lattice::Lattice(const std::vector<unsigned> &dimensions, const std::vector<BoundaryCondition> &boundary_conditions)
: dimensions_(dimensions), boundary_conditions_(boundary_conditions), site_count_(0)
{
    if (dimensions_.empty())
    {
        throw std::invalid_argument("dimensions must be non-empty");
    }

    if (dimensions_.size() != boundary_conditions_.size())
    {
        std::ostringstream oss;
        oss << "Mismatch between number of dimensions (" << dimensions_.size()
                << ") and number of boundary conditions specified (" << boundary_conditions_.size() << ").";
        throw std::invalid_argument(oss.str());
    }

    site_count_ = 1;
    for (int i = 0; i < dimensions_.size(); ++i)
    {
        if (dimensions_[i] <= 0)
        {
            std::ostringstream oss;
            oss << "Dimension size at axis " << i << " must be positive, got " << dimensions_[i] << ".";
            throw std::invalid_argument(oss.str());
        }
        site_count_ *= dimensions_[i];
    }

    // Should never run
    if (site_count_ == 0 && !dimensions_.empty()) {
        throw std::runtime_error("Calculated zero sites for a non-empty dimension set; overflow or invalid input.");
    }
}

