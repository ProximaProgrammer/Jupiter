#include "Grid.hpp"

#include <stdexcept>

namespace Jupiter
{

Grid::Grid(
    int nr,
    int ntheta,
    int nphi)
    :
    nr_(nr),
    ntheta_(ntheta),
    nphi_(nphi),
    cells_(nr * ntheta * nphi)
{
}

std::size_t Grid::Flatten(
    int ir,
    int itheta,
    int iphi) const
{
    return static_cast<std::size_t>(
        (ir * ntheta_ + itheta) * nphi_ + iphi
    );
}

Cell& Grid::At(
    int ir,
    int itheta,
    int iphi)
{
    return cells_.at(
        Flatten(ir, itheta, iphi)
    );
}

const Cell& Grid::At(
    int ir,
    int itheta,
    int iphi) const
{
    return cells_.at(
        Flatten(ir, itheta, iphi)
    );
}

int Grid::NR() const
{
    return nr_;
}

int Grid::NTheta() const
{
    return ntheta_;
}

int Grid::NPhi() const
{
    return nphi_;
}

}