#pragma once

#include "config/SimConfig.h"

#include <array>
#include <cstddef>
#include <vector>

namespace jrt {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double kBoltzmann = 1.380649e-23;
constexpr double speedOfLight = 299792458.0;

struct Vec3 {
    double r = 0.0;
    double th = 0.0;
    double ph = 0.0;
};

struct Cell {
    // Primitive state in the rotating Eulerian frame.
    double rho = 0.0;      // kg m^-3
    double p = 0.0;        // Pa
    double T = 0.0;        // K
    Vec3 u;                // spherical components in co-rotating frame [m s^-1]

    // Lightweight composition tracers by mass fraction.
    double x_h2 = 0.862;
    double x_he = 0.136;
    double x_ch4 = 2.0e-3;

    // Local opacity proxies used by the RT and shell exporter.
    double grey_opacity = 0.0;  // m^2 kg^-1
    double h3p_emissivity = 0.0;
};

struct Grid {
    SimConfig cfg;
    std::size_t nrTot = 0;
    std::size_t nthTot = 0;
    std::size_t nphTot = 0;
    double rInner = 0.0;
    double rOuter = 0.0;
    double dr = 0.0;
    double dtheta = 0.0;
    double dphi = 0.0;
    std::vector<Cell> cells;

    explicit Grid(const SimConfig& c);

    std::size_t idx(std::size_t i, std::size_t j, std::size_t k) const;
    Cell& at(std::size_t i, std::size_t j, std::size_t k);
    const Cell& at(std::size_t i, std::size_t j, std::size_t k) const;

    bool isInterior(std::size_t i, std::size_t j, std::size_t k) const;
    std::size_t iBegin() const { return cfg.ng; }
    std::size_t iEnd() const { return cfg.ng + cfg.nr; }
    std::size_t jBegin() const { return cfg.ng; }
    std::size_t jEnd() const { return cfg.ng + cfg.ntheta; }
    std::size_t kBegin() const { return 0; }
    std::size_t kEnd() const { return cfg.nphi; }
    std::size_t outerInteriorI() const { return cfg.ng + cfg.nr - 1; }
    std::size_t outerGhostI() const { return cfg.ng + cfg.nr; }

    double radius(std::size_t i) const;
    double theta(std::size_t j) const;
    double phi(std::size_t k) const;
    double latitude(std::size_t j) const;
    double cellVolume(std::size_t i, std::size_t j) const;

    void initializeJovianEnvelope();
    void applyBoundaryConditions();
    void updateDerivedOpacities();

    std::array<double, 3> cartesianPosition(std::size_t i, std::size_t j, std::size_t k) const;
    std::array<double, 3> inertialVelocityCartesian(std::size_t i, std::size_t j, std::size_t k) const;
};

} // namespace jrt
