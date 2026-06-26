#include "physics/Dynamics.h"
#include "physics/EOS.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace jrt {
namespace {

double sqr(double x) { return x * x; }

double speed(const Vec3& u) { return std::sqrt(u.r*u.r + u.th*u.th + u.ph*u.ph); }

double minCellLength(const Grid& g, std::size_t i, std::size_t j) {
    const double r = g.radius(i);
    const double s = std::max(0.08, std::sin(g.theta(j)));
    return std::min({g.dr, r * g.dtheta, r * s * g.dphi});
}

std::size_t km(const Grid& g, std::size_t k) { return (k == 0) ? (g.nphTot - 1) : (k - 1); }
std::size_t kp(const Grid& g, std::size_t k) { return (k + 1) % g.nphTot; }

double laplacianScalar(const Grid& g, const std::vector<Cell>& old, std::size_t i, std::size_t j, std::size_t k, char which) {
    const auto value = [&](std::size_t a, std::size_t b, std::size_t c) -> double {
        const Cell& cell = old[g.idx(a, b, c)];
        if (which == 'T') return cell.T;
        if (which == 'r') return cell.u.r;
        if (which == 't') return cell.u.th;
        return cell.u.ph;
    };
    const double r = g.radius(i);
    const double th = g.theta(j);
    const double sinth = std::max(0.08, std::sin(th));
    const double radial = (value(i+1,j,k) - 2.0*value(i,j,k) + value(i-1,j,k)) / (g.dr*g.dr);
    const double merid = (value(i,j+1,k) - 2.0*value(i,j,k) + value(i,j-1,k)) / (r*r*g.dtheta*g.dtheta);
    const double zonal = (value(i,j,kp(g,k)) - 2.0*value(i,j,k) + value(i,j,km(g,k))) / (r*r*sinth*sinth*g.dphi*g.dphi);
    return radial + merid + zonal;
}

} // namespace

StepSummary advanceDynamics(Grid& g, double requestedDt) {
    g.applyBoundaryConditions();
    g.updateDerivedOpacities();

    StepSummary summary;
    summary.dt = requestedDt;
    for (std::size_t i = g.iBegin(); i < g.iEnd(); ++i) {
        for (std::size_t j = g.jBegin(); j < g.jEnd(); ++j) {
            for (std::size_t k = g.kBegin(); k < g.kEnd(); ++k) {
                const Cell& c = g.at(i, j, k);
                const double cs = soundSpeedIdealGas(g.cfg.gamma, c.p, c.rho);
                const double sig = cs + speed(c.u);
                summary.maxSignalSpeed = std::max(summary.maxSignalSpeed, sig);
                summary.maxCorotatingSpeed = std::max(summary.maxCorotatingSpeed, speed(c.u));
                summary.dt = std::min(summary.dt, g.cfg.cfl * minCellLength(g, i, j) / std::max(1.0, sig));
                const auto vi = g.inertialVelocityCartesian(i, j, k);
                summary.maxInertialSpeed = std::max(summary.maxInertialSpeed, std::sqrt(sqr(vi[0]) + sqr(vi[1]) + sqr(vi[2])));
            }
        }
    }
    summary.dt = std::clamp(summary.dt, g.cfg.min_dt_s, g.cfg.max_dt_s);

    const std::vector<Cell> old = g.cells;
    const double dt = summary.dt;
    const double Omega = g.cfg.omega_rad_s;

    for (std::size_t i = g.iBegin(); i < g.iEnd(); ++i) {
        const double r = g.radius(i);
        const double invr = 1.0 / r;
        for (std::size_t j = g.jBegin(); j < g.jEnd(); ++j) {
            const double th = g.theta(j);
            const double sinth = std::max(0.08, std::sin(th));
            const double costh = std::cos(th);
            for (std::size_t k = g.kBegin(); k < g.kEnd(); ++k) {
                const Cell& c = old[g.idx(i, j, k)];
                Cell& n = g.at(i, j, k);

                const double dPdr = (old[g.idx(i+1,j,k)].p - old[g.idx(i-1,j,k)].p) / (2.0 * g.dr);
                const double dPdth = (old[g.idx(i,j+1,k)].p - old[g.idx(i,j-1,k)].p) / (2.0 * g.dtheta);
                const double dPdph = (old[g.idx(i,j,kp(g,k))].p - old[g.idx(i,j,km(g,k))].p) / (2.0 * g.dphi);

                const double drhodr = (old[g.idx(i+1,j,k)].rho - old[g.idx(i-1,j,k)].rho) / (2.0 * g.dr);
                const double drhodth = (old[g.idx(i,j+1,k)].rho - old[g.idx(i,j-1,k)].rho) / (2.0 * g.dtheta);
                const double drhodph = (old[g.idx(i,j,kp(g,k))].rho - old[g.idx(i,j,km(g,k))].rho) / (2.0 * g.dphi);

                const double dTdr = (old[g.idx(i+1,j,k)].T - old[g.idx(i-1,j,k)].T) / (2.0 * g.dr);
                const double dTdth = (old[g.idx(i,j+1,k)].T - old[g.idx(i,j-1,k)].T) / (2.0 * g.dtheta);
                const double dTdph = (old[g.idx(i,j,kp(g,k))].T - old[g.idx(i,j,km(g,k))].T) / (2.0 * g.dphi);

                const double divu =
                    (sqr(r + 0.5*g.dr) * old[g.idx(i+1,j,k)].u.r - sqr(r - 0.5*g.dr) * old[g.idx(i-1,j,k)].u.r) / (2.0 * g.dr * r * r)
                    + ((std::sin(th + 0.5*g.dtheta) * old[g.idx(i,j+1,k)].u.th - std::sin(th - 0.5*g.dtheta) * old[g.idx(i,j-1,k)].u.th) / (2.0 * g.dtheta * r * sinth))
                    + (old[g.idx(i,j,kp(g,k))].u.ph - old[g.idx(i,j,km(g,k))].u.ph) / (2.0 * g.dphi * r * sinth);

                const double rho = std::max(1.0e-10, c.rho);
                const double a_cent_th = Omega*Omega * r * sinth * costh;

                // Background-state split: the analytic radial pressure profile
                // carries the enormous near-hydrostatic gravity/centrifugal balance.
                // The explicit velocity update responds to deviations from that
                // background instead of launching a global artificial collapse.
                const double pDeep = 3.0e6;
                const double dP0dr = -c.p * std::log(pDeep / g.cfg.top_pressure_pa) / (g.rOuter - g.rInner);

                // Coriolis acceleration in spherical components for a frame rotating
                // about +z: -2 Omega x u.
                const double a_cor_r  =  2.0 * Omega * sinth * c.u.ph;
                const double a_cor_th =  2.0 * Omega * costh * c.u.ph;
                const double a_cor_ph = -2.0 * Omega * (costh * c.u.th + sinth * c.u.r);

                Vec3 a;
                a.r  = -(dPdr - dP0dr) / rho + a_cor_r;
                a.th = -(dPdth * invr) / rho + a_cent_th + a_cor_th;
                a.ph = -(dPdph * invr / sinth) / rho + a_cor_ph;

                // Explicit turbulent viscosity and Rayleigh drag keep the first usable
                // physics patch stable while still letting jets/shear evolve.
                a.r  += g.cfg.kinematic_viscosity_m2_s * laplacianScalar(g, old, i, j, k, 'r') - c.u.r / g.cfg.velocity_drag_time_s;
                a.th += g.cfg.kinematic_viscosity_m2_s * laplacianScalar(g, old, i, j, k, 't') - c.u.th / g.cfg.velocity_drag_time_s;
                a.ph += g.cfg.kinematic_viscosity_m2_s * laplacianScalar(g, old, i, j, k, 'p') - c.u.ph / g.cfg.velocity_drag_time_s;

                n.u.r  = c.u.r  + dt * a.r;
                n.u.th = c.u.th + dt * a.th;
                n.u.ph = c.u.ph + dt * a.ph;

                // Upwind-lite primitive advection plus compressive heating/cooling.
                const double advRho = c.u.r * drhodr + c.u.th * invr * drhodth + c.u.ph * invr / sinth * drhodph;
                n.rho = std::max(1.0e-10, c.rho + dt * (-advRho - c.rho * divu));

                const double advT = c.u.r * dTdr + c.u.th * invr * dTdth + c.u.ph * invr / sinth * dTdph;
                const double lapT = laplacianScalar(g, old, i, j, k, 'T');
                const double eta = std::clamp((r - g.rInner) / (g.rOuter - g.rInner), 0.0, 1.0);
                const double tauRad = g.cfg.radiative_time_top_s * eta + g.cfg.radiative_time_deep_s * (1.0 - eta);
                const double Teq = equilibriumTemperatureAtRadius(r, g.rInner, g.rOuter, g.cfg.top_temperature_k, g.cfg.deep_temperature_k);
                const double dTdt = -advT - (g.cfg.gamma - 1.0) * c.T * divu
                                  + g.cfg.thermal_diffusivity_m2_s * lapT
                                  - (c.T - Teq) / std::max(1.0, tauRad);
                n.T = std::clamp(c.T + dt * dTdt, 45.0, 6000.0);
                n.p = pressureIdealGas(n.rho, n.T, g.cfg.mean_molecular_weight_kg);

                // Bound velocities to keep first-order explicit steps from exploding if
                // the user experiments with very coarse grids or very large dt.
                const double cs = soundSpeedIdealGas(g.cfg.gamma, n.p, n.rho);
                const double vmax = 0.80 * cs;
                const double sp = speed(n.u);
                if (sp > vmax && sp > 0.0) {
                    const double q = vmax / sp;
                    n.u.r *= q; n.u.th *= q; n.u.ph *= q;
                }
            }
        }
    }

    g.applyBoundaryConditions();
    g.updateDerivedOpacities();
    return summary;
}

} // namespace jrt
