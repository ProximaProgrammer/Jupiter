#pragma once

#include "../core/Grid.hpp"

namespace Jupiter
{

struct DynamicsConfig
{
    double epsilon_max = 0.025;
    bool use_finite_volume_advection = true;

    bool include_gravity = true;
    bool include_pressure_gradient = true;
    bool rotating_frame_sources = true;
    bool include_viscosity = true;
    bool include_thermal_diffusion = true;
    bool include_convective_buoyancy = true;
    bool include_radiative_relaxation = true;
    bool include_chemistry = true;
    bool enforce_boundaries = true;

    double kinematic_viscosity_m2_s = 2.5e7;
    double thermal_diffusivity_m2_s = 2.5e7;
    double convective_buoyancy_strength = 0.08;
    double velocity_drag_timescale_s = 2.5e7;
    double top_radiative_timescale_s = 2.0e5;
    double deep_radiative_timescale_s = 2.0e9;
    double chemistry_timescale_s = 5.0e4;
};

class Dynamics
{
public:
    static void UpdatePrimitives(Grid& grid);
    static void ApplySourceTerms(Grid& grid, double dt, const DynamicsConfig& cfg = {});
    static void AdvectFiniteVolume(Grid& grid, double dt, const DynamicsConfig& cfg = {});
    static void ShiftEpsilon(Grid& grid, double dt, const DynamicsConfig& cfg = {});
    static void ApplyBoundaryConditions(Grid& grid, const DynamicsConfig& cfg = {});
};

} // namespace Jupiter
