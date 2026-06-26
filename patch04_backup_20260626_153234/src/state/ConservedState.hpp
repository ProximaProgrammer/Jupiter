#pragma once

#include <array>
#include <algorithm>
#include "Species.hpp"

namespace Jupiter
{

struct ConservedState
{
    double mass = 0.0;
    double momentum_r = 0.0;
    double momentum_theta = 0.0;
    double momentum_phi = 0.0;
    double total_energy = 0.0;
    double magnetic_r = 0.0;
    double magnetic_theta = 0.0;
    double magnetic_phi = 0.0;
    std::array<double, NumSpecies> species_mass{};

    ConservedState& operator+=(const ConservedState& b)
    {
        mass += b.mass;
        momentum_r += b.momentum_r;
        momentum_theta += b.momentum_theta;
        momentum_phi += b.momentum_phi;
        total_energy += b.total_energy;
        magnetic_r += b.magnetic_r;
        magnetic_theta += b.magnetic_theta;
        magnetic_phi += b.magnetic_phi;
        for (std::size_t i = 0; i < NumSpecies; ++i) species_mass[i] += b.species_mass[i];
        return *this;
    }

    ConservedState& operator-=(const ConservedState& b)
    {
        mass -= b.mass;
        momentum_r -= b.momentum_r;
        momentum_theta -= b.momentum_theta;
        momentum_phi -= b.momentum_phi;
        total_energy -= b.total_energy;
        magnetic_r -= b.magnetic_r;
        magnetic_theta -= b.magnetic_theta;
        magnetic_phi -= b.magnetic_phi;
        for (std::size_t i = 0; i < NumSpecies; ++i) species_mass[i] -= b.species_mass[i];
        return *this;
    }
};

inline ConservedState operator+(ConservedState a, const ConservedState& b) { a += b; return a; }
inline ConservedState operator-(ConservedState a, const ConservedState& b) { a -= b; return a; }

inline ConservedState operator*(const ConservedState& a, double f)
{
    ConservedState r;
    r.mass = a.mass * f;
    r.momentum_r = a.momentum_r * f;
    r.momentum_theta = a.momentum_theta * f;
    r.momentum_phi = a.momentum_phi * f;
    r.total_energy = a.total_energy * f;
    r.magnetic_r = a.magnetic_r * f;
    r.magnetic_theta = a.magnetic_theta * f;
    r.magnetic_phi = a.magnetic_phi * f;
    for (std::size_t i = 0; i < NumSpecies; ++i) r.species_mass[i] = a.species_mass[i] * f;
    return r;
}

inline ConservedState operator*(double f, const ConservedState& a) { return a * f; }

inline void ClipConserved(ConservedState& u)
{
    if (u.mass < 0.0) u.mass = 0.0;
    if (u.total_energy < 0.0) u.total_energy = 0.0;
    for (double& m : u.species_mass) if (m < 0.0) m = 0.0;
}

} // namespace Jupiter
