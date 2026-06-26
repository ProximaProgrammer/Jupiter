#pragma once

#include <array>
#include "Species.hpp"

namespace Jupiter
{

struct PrimitiveState
{
    double density = 0.0;
    double pressure = 0.0;
    double temperature = 0.0;
    double velocity_r = 0.0;
    double velocity_theta = 0.0;
    double velocity_phi = 0.0;
    double sound_speed = 0.0;
    double gamma_eff = 1.4;
    double mean_molar_mass = 2.3e-3;
    double cold_pressure = 0.0;
    double electron_degeneracy_pressure = 0.0;
    double magnetic_r = 0.0;
    double magnetic_theta = 0.0;
    double magnetic_phi = 0.0;
    std::array<double, NumSpecies> mass_fraction{};
};

} // namespace Jupiter
