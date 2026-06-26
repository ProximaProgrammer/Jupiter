#pragma once

#include "../core/Grid.hpp"

namespace Jupiter
{

struct DynamicsConfig
{
    double epsilon_max = 0.01;
    bool include_gravity = true;
    bool include_pressure_gradient = true;
    bool rotating_frame_sources = true;
};

class Dynamics
{
public:
    static void UpdatePrimitives(Grid& grid);
    static void ApplySourceTerms(Grid& grid, double dt, const DynamicsConfig& cfg = {});
    static void ShiftEpsilon(Grid& grid, double dt, const DynamicsConfig& cfg = {});
};

} // namespace Jupiter
