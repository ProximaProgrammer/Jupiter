#include "Grid.hpp"
#include "../common/Constants.hpp"

#include <cmath>
#include <stdexcept>

namespace Jupiter
{

Grid::Grid(int nr, int ntheta, int nphi, double r_inner, double r_outer)
    : nr_(nr), ntheta_(ntheta), nphi_(nphi), cells_(static_cast<std::size_t>(nr) * ntheta * nphi)
{
    if (nr <= 0 || ntheta <= 0 || nphi <= 0) throw std::invalid_argument("Grid dimensions must be positive");
    if (!(r_outer > r_inner && r_inner >= 0.0)) throw std::invalid_argument("Invalid radial bounds");

    const double dtheta = Constants::pi / static_cast<double>(ntheta);
    const double dphi = 2.0 * Constants::pi / static_cast<double>(nphi);

    // Log radial spacing gives useful resolution in the observable exterior while still spanning deep layers.
    const bool log_r = r_inner > 0.0;
    const double log_a = log_r ? std::log(r_inner) : 0.0;
    const double log_b = log_r ? std::log(r_outer) : 0.0;

    for (int ir = 0; ir < nr_; ++ir)
    {
        double r0 = 0.0, r1 = 0.0;
        if (log_r)
        {
            r0 = std::exp(log_a + (log_b - log_a) * static_cast<double>(ir) / nr_);
            r1 = std::exp(log_a + (log_b - log_a) * static_cast<double>(ir + 1) / nr_);
        }
        else
        {
            r0 = r_inner + (r_outer - r_inner) * static_cast<double>(ir) / nr_;
            r1 = r_inner + (r_outer - r_inner) * static_cast<double>(ir + 1) / nr_;
        }

        for (int it = 0; it < ntheta_; ++it)
        {
            const double t0 = dtheta * it;
            const double t1 = dtheta * (it + 1);
            const double cos0 = std::cos(t0);
            const double cos1 = std::cos(t1);
            const double sin0 = std::sin(t0);
            const double sin1 = std::sin(t1);

            for (int ip = 0; ip < nphi_; ++ip)
            {
                const double p0 = dphi * ip;
                const double p1 = dphi * (ip + 1);
                Cell& c = At(ir, it, ip);
                c.index = {ir, it, ip};
                auto& g = c.geometry;
                g.r_inner = r0;
                g.r_outer = r1;
                g.theta_lower = t0;
                g.theta_upper = t1;
                g.phi_lower = p0;
                g.phi_upper = p1;
                g.r_center = 0.5 * (r0 + r1);
                g.theta_center = 0.5 * (t0 + t1);
                g.phi_center = 0.5 * (p0 + p1);
                g.dr = r1 - r0;
                g.dtheta = dtheta;
                g.dphi = dphi;
                g.volume = (p1 - p0) * (cos0 - cos1) * (r1*r1*r1 - r0*r0*r0) / 3.0;
                g.area_r_minus = r0*r0 * (p1 - p0) * (cos0 - cos1);
                g.area_r_plus = r1*r1 * (p1 - p0) * (cos0 - cos1);
                g.area_theta_minus = 0.5 * (r1*r1 - r0*r0) * std::max(0.0, sin0) * (p1 - p0);
                g.area_theta_plus = 0.5 * (r1*r1 - r0*r0) * std::max(0.0, sin1) * (p1 - p0);
                g.area_phi_minus = 0.5 * (r1*r1 - r0*r0) * (t1 - t0);
                g.area_phi_plus = g.area_phi_minus;
                g.length_r = g.dr;
                g.length_theta = std::max(1.0, g.r_center * dtheta);
                g.length_phi = std::max(1.0, g.r_center * std::sin(g.theta_center) * dphi);
            }
        }
    }
}

bool Grid::InBounds(int ir, int itheta, int iphi) const
{
    return ir >= 0 && ir < nr_ && itheta >= 0 && itheta < ntheta_ && iphi >= 0 && iphi < nphi_;
}

std::size_t Grid::Flatten(int ir, int itheta, int iphi) const
{
    if (!InBounds(ir, itheta, iphi)) throw std::out_of_range("Grid index out of bounds");
    return static_cast<std::size_t>((ir * ntheta_ + itheta) * nphi_ + iphi);
}

Cell& Grid::At(int ir, int itheta, int iphi) { return cells_.at(Flatten(ir, itheta, iphi)); }
const Cell& Grid::At(int ir, int itheta, int iphi) const { return cells_.at(Flatten(ir, itheta, iphi)); }

} // namespace Jupiter
