#pragma once

#include <vector>

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

    std::vector<double> species_mass;
};

}