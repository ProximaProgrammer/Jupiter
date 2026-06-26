#include "core/Grid.h"
#include "physics/EOS.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace jrt {

namespace {

double clampPositive(double x, double floor) { return std::max(x, floor); }

} // namespace

Grid::Grid(const SimConfig& c) : cfg(c) {
    nrTot = cfg.nr + 2 * cfg.ng;
    nthTot = cfg.ntheta + 2 * cfg.ng;
    nphTot = cfg.nphi;
    rOuter = cfg.r_outer_m;
    rInner = cfg.r_inner_fraction * cfg.r_outer_m;
    dr = (rOuter - rInner) / static_cast<double>(cfg.nr);
    dtheta = pi / static_cast<double>(cfg.ntheta);
    dphi = 2.0 * pi / static_cast<double>(cfg.nphi);
    cells.assign(nrTot * nthTot * nphTot, Cell{});
}

std::size_t Grid::idx(std::size_t i, std::size_t j, std::size_t k) const {
    k %= nphTot;
    return (i * nthTot + j) * nphTot + k;
}

Cell& Grid::at(std::size_t i, std::size_t j, std::size_t k) { return cells[idx(i, j, k)]; }
const Cell& Grid::at(std::size_t i, std::size_t j, std::size_t k) const { return cells[idx(i, j, k)]; }

bool Grid::isInterior(std::size_t i, std::size_t j, std::size_t k) const {
    return i >= iBegin() && i < iEnd() && j >= jBegin() && j < jEnd() && k < kEnd();
}

double Grid::radius(std::size_t i) const {
    const double ii = static_cast<double>(i) - static_cast<double>(cfg.ng) + 0.5;
    return rInner + ii * dr;
}

double Grid::theta(std::size_t j) const {
    const double jj = static_cast<double>(j) - static_cast<double>(cfg.ng) + 0.5;
    return std::clamp(jj * dtheta, 0.5 * dtheta, pi - 0.5 * dtheta);
}

double Grid::phi(std::size_t k) const {
    return (static_cast<double>(k) + 0.5) * dphi;
}

double Grid::latitude(std::size_t j) const { return 0.5 * pi - theta(j); }

double Grid::cellVolume(std::size_t i, std::size_t j) const {
    const double r0 = radius(i) - 0.5 * dr;
    const double r1 = radius(i) + 0.5 * dr;
    const double th0 = theta(j) - 0.5 * dtheta;
    const double th1 = theta(j) + 0.5 * dtheta;
    return (r1*r1*r1 - r0*r0*r0) / 3.0 * (std::cos(th0) - std::cos(th1)) * dphi;
}

void Grid::initializeJovianEnvelope() {
    // Temperature/composition include zonal belts and waves so the second outer shell
    // is visually meaningful even before long-time dynamics is added.
    for (std::size_t i = 0; i < nrTot; ++i) {
        const double r = radius(i);
        const double eta = std::clamp((r - rInner) / (rOuter - rInner), 0.0, 1.0);
        const double Tbase = cfg.top_temperature_k + (cfg.deep_temperature_k - cfg.top_temperature_k) * std::pow(1.0 - eta, 1.15);
        for (std::size_t j = 0; j < nthTot; ++j) {
            const double lat = latitude(j);
            const double band = std::sin(10.0 * lat) + 0.35 * std::sin(22.0 * lat + 0.3);
            const double jet = cfg.jet_speed_m_s * std::sin(13.0 * lat) * std::exp(-0.8 * (1.0 - eta));
            for (std::size_t k = 0; k < nphTot; ++k) {
                const double ph = phi(k);
                Cell& c = at(i, j, k);
                const double wave = 0.012 * std::sin(4.0 * ph + 3.0 * lat) + 0.006 * std::cos(7.0 * ph - 2.0 * lat);
                c.T = clampPositive(Tbase * (1.0 + 0.018 * band + wave), 45.0);
                c.x_h2 = 0.862;
                c.x_he = 0.136;
                c.x_ch4 = std::max(2.0e-4, 2.1e-3 * (1.0 + 0.18 * band + 0.08 * std::sin(3.0 * ph)));
                c.u.r = cfg.perturbation_speed_m_s * 0.20 * std::sin(2.0 * ph) * std::cos(lat) * std::exp(-2.0 * eta);
                c.u.th = cfg.perturbation_speed_m_s * 0.15 * std::cos(3.0 * ph + lat) * std::sin(2.0 * lat);
                c.u.ph = jet + cfg.perturbation_speed_m_s * std::sin(5.0 * ph + 2.0 * lat) * std::cos(lat);
            }
        }
    }

    // Smooth pressure profile for the outer Jovian envelope. A full hydrostatic
    // table will replace this later, but this avoids an unstable runaway from
    // explicit one-cell hydrostatic integration over hundreds of km per cell.
    const double pDeep = 3.0e6; // Pa at r_inner for this outer-shell patch
    const double logRatio = std::log(pDeep / cfg.top_pressure_pa);
    for (std::size_t i = iBegin(); i < iEnd(); ++i) {
        const double eta = std::clamp((radius(i) - rInner) / (rOuter - rInner), 0.0, 1.0);
        const double pBase = cfg.top_pressure_pa * std::exp((1.0 - eta) * logRatio);
        for (std::size_t j = jBegin(); j < jEnd(); ++j) {
            const double lat = latitude(j);
            const double bandP = 1.0 + 0.006 * std::sin(10.0 * lat);
            for (std::size_t k = kBegin(); k < kEnd(); ++k) {
                Cell& c = at(i, j, k);
                c.p = std::max(10.0, pBase * bandP);
                c.rho = densityIdealGas(c.p, c.T, cfg.mean_molecular_weight_kg);
            }
        }
    }

    applyBoundaryConditions();
    updateDerivedOpacities();
}

void Grid::applyBoundaryConditions() {
    // Radial ghosts: inner is reflecting for u_r; outer is open/radiative and is
    // intentionally the layer not visualized by the shell exporter.
    for (std::size_t j = jBegin(); j < jEnd(); ++j) {
        for (std::size_t k = kBegin(); k < kEnd(); ++k) {
            for (std::size_t g = 0; g < cfg.ng; ++g) {
                Cell inner = at(cfg.ng, j, k);
                inner.u.r = std::max(0.0, -inner.u.r);
                at(g, j, k) = inner;

                Cell outer = at(outerInteriorI(), j, k);
                outer.u.r = std::min(0.0, outer.u.r);
                outer.T = 0.985 * outer.T + 0.015 * cfg.top_temperature_k;
                outer.p = std::max(0.25 * cfg.top_pressure_pa, outer.p);
                outer.rho = densityIdealGas(outer.p, outer.T, cfg.mean_molecular_weight_kg);
                at(outerGhostI() + g, j, k) = outer;
            }
        }
    }

    // Polar ghosts: reflect meridional velocity. This avoids tiny singular cells
    // while preserving a spherical coordinate layout.
    for (std::size_t i = 0; i < nrTot; ++i) {
        for (std::size_t k = 0; k < nphTot; ++k) {
            Cell north = at(i, cfg.ng, k);
            north.u.th *= -1.0;
            at(i, 0, k) = north;
            Cell south = at(i, cfg.ng + cfg.ntheta - 1, k);
            south.u.th *= -1.0;
            at(i, cfg.ng + cfg.ntheta, k) = south;
        }
    }
}

void Grid::updateDerivedOpacities() {
    for (std::size_t i = iBegin(); i < iEnd(); ++i) {
        const double eta = std::clamp((radius(i) - rInner) / (rOuter - rInner), 0.0, 1.0);
        for (std::size_t j = jBegin(); j < jEnd(); ++j) {
            const double lat = latitude(j);
            for (std::size_t k = kBegin(); k < kEnd(); ++k) {
                Cell& c = at(i, j, k);
                c.rho = densityIdealGas(c.p, c.T, cfg.mean_molecular_weight_kg);
                const double cloud = 1.0 + 0.35 * std::abs(std::sin(10.0 * lat));
                c.grey_opacity = 0.006 * cloud + 0.020 * c.x_ch4 + 0.0025 * std::sqrt(std::max(0.0, c.rho));
                c.h3p_emissivity = std::exp(-std::pow((eta - 0.98) / 0.07, 2.0)) * (1.0 + 0.5 * std::abs(std::sin(lat)));
            }
        }
    }
}

std::array<double, 3> Grid::cartesianPosition(std::size_t i, std::size_t j, std::size_t k) const {
    const double r = radius(i);
    const double th = theta(j);
    const double ph = phi(k);
    return {r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph), r * std::cos(th)};
}

std::array<double, 3> Grid::inertialVelocityCartesian(std::size_t i, std::size_t j, std::size_t k) const {
    const Cell& c = at(i, j, k);
    const double th = theta(j);
    const double ph = phi(k);
    const std::array<double, 3> er = {std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th)};
    const std::array<double, 3> eth = {std::cos(th) * std::cos(ph), std::cos(th) * std::sin(ph), -std::sin(th)};
    const std::array<double, 3> eph = {-std::sin(ph), std::cos(ph), 0.0};

    std::array<double, 3> v = {
        c.u.r * er[0] + c.u.th * eth[0] + c.u.ph * eph[0],
        c.u.r * er[1] + c.u.th * eth[1] + c.u.ph * eph[1],
        c.u.r * er[2] + c.u.th * eth[2] + c.u.ph * eph[2]
    };
    // Add v_rot = Omega x r. This is the Doppler-relevant inertial velocity.
    const auto x = cartesianPosition(i, j, k);
    v[0] += -cfg.omega_rad_s * x[1];
    v[1] +=  cfg.omega_rad_s * x[0];
    return v;
}

} // namespace jrt
