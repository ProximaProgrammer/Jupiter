#include "Dynamics.hpp"
#include "EOS.hpp"
#include "../common/Constants.hpp"
#include "../common/Vec3.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Jupiter
{

void Dynamics::UpdatePrimitives(Grid& grid)
{
    for (Cell& c : grid.Cells()) c.primitive = EOS::ConservedToPrimitive(c.conserved, c.geometry).primitive;
}

static double pressure_at(const Grid& grid, int ir, int it, int ip)
{
    return grid.At(ir, it, ip).primitive.pressure;
}

static int wrap_phi(int ip, int nphi)
{
    if (ip < 0) return ip + nphi;
    if (ip >= nphi) return ip - nphi;
    return ip;
}

void Dynamics::ApplySourceTerms(Grid& grid, double dt, const DynamicsConfig& cfg)
{
    UpdatePrimitives(grid);
    std::vector<ConservedState> delta(grid.Size());

    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const Cell& c = grid.At(ir, it, ip);
        const auto& g = c.geometry;
        const auto& p = c.primitive;
        if (c.conserved.mass <= 0.0) continue;

        double ar = 0.0, at = 0.0, ap = 0.0;

        if (cfg.include_gravity)
        {
            ar += -Constants::G * Constants::M_J / std::max(g.r_center*g.r_center, 1.0);
        }

        if (cfg.include_pressure_gradient)
        {
            const int irm = std::max(0, ir - 1);
            const int irp = std::min(grid.NR() - 1, ir + 1);
            const int itm = std::max(0, it - 1);
            const int itp = std::min(grid.NTheta() - 1, it + 1);
            const int ipm = wrap_phi(ip - 1, grid.NPhi());
            const int ipp = wrap_phi(ip + 1, grid.NPhi());

            const double dpdr = (pressure_at(grid, irp, it, ip) - pressure_at(grid, irm, it, ip)) /
                std::max(grid.At(irp, it, ip).geometry.r_center - grid.At(irm, it, ip).geometry.r_center, 1.0);
            const double dpdt = (pressure_at(grid, ir, itp, ip) - pressure_at(grid, ir, itm, ip)) /
                std::max(grid.At(ir, itp, ip).geometry.theta_center - grid.At(ir, itm, ip).geometry.theta_center, 1.0e-12);
            // Periodic phi coordinate: use the local angular step instead of wrapped center difference.
            const double dpdp = (pressure_at(grid, ir, it, ipp) - pressure_at(grid, ir, it, ipm)) /
                std::max(2.0 * g.dphi, 1.0e-12);

            ar += -dpdr / std::max(p.density, 1.0e-30);
            at += -(1.0 / std::max(g.r_center, 1.0)) * dpdt / std::max(p.density, 1.0e-30);
            ap += -(1.0 / std::max(g.r_center * std::sin(g.theta_center), 1.0)) * dpdp / std::max(p.density, 1.0e-30);
        }

        if (cfg.rotating_frame_sources)
        {
            // Centrifugal acceleration: Omega^2 r sin(theta) [sin(theta)e_r + cos(theta)e_theta].
            const double st = std::sin(g.theta_center), ct = std::cos(g.theta_center);
            ar += Constants::Omega_J * Constants::Omega_J * g.r_center * st * st;
            at += Constants::Omega_J * Constants::Omega_J * g.r_center * st * ct;

            // Coriolis acceleration in local basis from -2 Omega x u.
            const Vec3 u_cart = SphericalVelocityToCartesian(g.theta_center, g.phi_center,
                                                            p.velocity_r, p.velocity_theta, p.velocity_phi);
            const Vec3 omega{0.0, 0.0, Constants::Omega_J};
            const Vec3 a_cor = -2.0 * Cross(omega, u_cart);
            ar += Dot(a_cor, SphericalBasisR(g.theta_center, g.phi_center));
            at += Dot(a_cor, SphericalBasisTheta(g.theta_center, g.phi_center));
            ap += Dot(a_cor, SphericalBasisPhi(g.phi_center));
        }

        ConservedState d;
        d.momentum_r = c.conserved.mass * ar * dt;
        d.momentum_theta = c.conserved.mass * at * dt;
        d.momentum_phi = c.conserved.mass * ap * dt;
        const double power = c.conserved.mass * (p.velocity_r * ar + p.velocity_theta * at + p.velocity_phi * ap);
        d.total_energy = power * dt;
        delta[grid.Flatten(ir, it, ip)] += d;
    }

    for (std::size_t i = 0; i < grid.Size(); ++i)
    {
        grid.Cells()[i].conserved += delta[i];
        ClipConserved(grid.Cells()[i].conserved);
    }
    UpdatePrimitives(grid);
}

static double transfer_fraction(double v, double dt, double length, double eps_max)
{
    if (v <= 0.0) return 0.0;
    return std::clamp(v * dt / std::max(length, 1.0), 0.0, eps_max);
}

static void transfer(Grid& grid, std::vector<ConservedState>& delta,
                     int ir0, int it0, int ip0, int ir1, int it1, int ip1, double eps)
{
    if (eps <= 0.0 || !grid.InBounds(ir1, it1, ip1)) return;
    const std::size_t a = grid.Flatten(ir0, it0, ip0);
    const std::size_t b = grid.Flatten(ir1, it1, ip1);
    const ConservedState moved = grid.Cells()[a].conserved * eps;
    delta[a] -= moved;
    delta[b] += moved;
}

void Dynamics::ShiftEpsilon(Grid& grid, double dt, const DynamicsConfig& cfg)
{
    UpdatePrimitives(grid);
    std::vector<ConservedState> delta(grid.Size());

    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const Cell& c = grid.At(ir, it, ip);
        const auto& p = c.primitive;
        const auto& g = c.geometry;
        if (c.conserved.mass <= 0.0) continue;

        // Face-neighbor conservative shifts only. Corner-adjacent terms are intentionally ignored O(epsilon^2).
        const double erp = transfer_fraction( p.velocity_r, dt, g.length_r, cfg.epsilon_max);
        const double erm = transfer_fraction(-p.velocity_r, dt, g.length_r, cfg.epsilon_max);
        transfer(grid, delta, ir, it, ip, ir + 1, it, ip, erp);
        transfer(grid, delta, ir, it, ip, ir - 1, it, ip, erm);

        const double etp = transfer_fraction( p.velocity_theta, dt, g.length_theta, cfg.epsilon_max);
        const double etm = transfer_fraction(-p.velocity_theta, dt, g.length_theta, cfg.epsilon_max);
        transfer(grid, delta, ir, it, ip, ir, it + 1, ip, etp);
        transfer(grid, delta, ir, it, ip, ir, it - 1, ip, etm);

        const int ipp = wrap_phi(ip + 1, grid.NPhi());
        const int ipm = wrap_phi(ip - 1, grid.NPhi());
        const double epp = transfer_fraction( p.velocity_phi, dt, g.length_phi, cfg.epsilon_max);
        const double epm = transfer_fraction(-p.velocity_phi, dt, g.length_phi, cfg.epsilon_max);
        transfer(grid, delta, ir, it, ip, ir, it, ipp, epp);
        transfer(grid, delta, ir, it, ip, ir, it, ipm, epm);
    }

    for (std::size_t i = 0; i < grid.Size(); ++i)
    {
        grid.Cells()[i].conserved += delta[i];
        ClipConserved(grid.Cells()[i].conserved);
    }
    UpdatePrimitives(grid);
}

} // namespace Jupiter
