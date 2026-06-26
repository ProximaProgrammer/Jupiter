#pragma once

#include <string>
#include "../common/Constants.hpp"
#include "../physics/Dynamics.hpp"
#include "../physics/RadiativeTransfer.hpp"
#include "../pipeline/InitialConditions.hpp"

namespace Jupiter
{

struct SimConfig
{
    int nr = 96;
    int ntheta = 64;
    int nphi = 128;

    double r_inner_fraction = 0.20;
    double r_outer_fraction = 1.00;

    double requested_dt_seconds = 0.10;
    int steps = 120;
    bool use_cfl_dt = true;
    double cfl_safety = 0.35;

    // Conservative transport. epsilon_max is an extra safety limiter on the
    // fraction of a donor cell that can cross one face in a single explicit step.
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

    // These are deliberately eddy/effective parameters for a coarse global grid,
    // not molecular values. They are used to make the simulation visibly evolve
    // while remaining stable on a laptop-sized grid.
    double kinematic_viscosity_m2_s = 2.5e7;
    double thermal_diffusivity_m2_s = 2.5e7;
    double convective_buoyancy_strength = 0.08;
    double velocity_drag_timescale_s = 2.5e7;
    double top_radiative_timescale_s = 2.0e5;
    double deep_radiative_timescale_s = 2.0e9;
    double chemistry_timescale_s = 5.0e4;

    InitialConditionConfig initial;

    std::string raw_exomol_dir = "data/raw/exomol";
    std::string raw_cia_dir = "data/raw/cia";
    std::string runtime_asset_dir = "data/runtime";
    bool use_runtime_tables = true;

    bool include_doppler = true;
    bool include_rotation_velocity = true;
    bool use_cuda_rt = false;
    double opacity_scale = 1.0;
    double h3plus_emission_scale = 1.0;
    int ray_theta_index = -1;
    int ray_phi_index = 0;

    std::string output_dir = "data/output";
    bool clean_output_on_start = true;
    bool write_spectrum = true;
    bool write_diagnostics = true;
    bool visualizer_enabled = true;
    int visualizer_every = 4;
    int visualizer_phi_index = 0;
    std::string visualizer_format = "ppm";
    // Comma-separated: temperature,velocity,density,pressure,entropy,opacity_proxy.
    std::string visualizer_modes = "temperature,velocity,density,pressure";
    bool visualizer_dashboard = true;
    int visualizer_panel_width = 320;
    int visualizer_panel_height = 480;
    int visualizer_margin = 24;
    bool visualizer_write_individual = false;
    bool visualizer_autoscale_each_frame = true;

    double RInnerMeters() const { return r_inner_fraction * Constants::R_J; }
    double ROuterMeters() const { return r_outer_fraction * Constants::R_J; }

    DynamicsConfig Dynamics() const;
    RTConfig RT() const;
};

SimConfig LoadSimConfig(const std::string& path);
void PrintSimConfig(const SimConfig& cfg);

} // namespace Jupiter
