#pragma once

#include <cstddef>
#include <vector>
#include "Cell.hpp"

namespace Jupiter
{

class Grid
{
public:
    Grid(int nr, int ntheta, int nphi, double r_inner, double r_outer);

    Cell& At(int ir, int itheta, int iphi);
    const Cell& At(int ir, int itheta, int iphi) const;

    bool InBounds(int ir, int itheta, int iphi) const;
    std::size_t Flatten(int ir, int itheta, int iphi) const;

    int NR() const { return nr_; }
    int NTheta() const { return ntheta_; }
    int NPhi() const { return nphi_; }
    std::size_t Size() const { return cells_.size(); }

    std::vector<Cell>& Cells() { return cells_; }
    const std::vector<Cell>& Cells() const { return cells_; }

private:
    int nr_ = 0;
    int ntheta_ = 0;
    int nphi_ = 0;
    std::vector<Cell> cells_;
};

} // namespace Jupiter
