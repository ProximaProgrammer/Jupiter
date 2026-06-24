#pragma once

#include <vector>

#include "Cell.hpp"

namespace Jupiter
{

class Grid
{
public:

    Grid(
        int nr,
        int ntheta,
        int nphi);

    Cell& At(
        int ir,
        int itheta,
        int iphi);

    const Cell& At(
        int ir,
        int itheta,
        int iphi) const;

    int NR() const;
    int NTheta() const;
    int NPhi() const;

private:

    int nr_;
    int ntheta_;
    int nphi_;

    std::vector<Cell> cells_;

    std::size_t Flatten(
        int ir,
        int itheta,
        int iphi) const;
};

}