#pragma once

#include "../core/Grid.hpp"

namespace Jupiter
{

struct InitialConditionConfig
{
    double surface_temperature_K = 165.0;
    double deep_temperature_K = 2.5e4;
    double surface_pressure_Pa = 1.0e5;
    double deep_pressure_Pa = 5.0e12;
    double surface_density_kg_m3 = 0.16;
    double deep_density_kg_m3 = 1400.0;
};

class InitialConditions
{
public:
    static void FillJupiterLike(Grid& grid, const InitialConditionConfig& cfg = {});
};

} // namespace Jupiter
