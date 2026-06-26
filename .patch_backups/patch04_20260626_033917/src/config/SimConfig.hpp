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

    double requested_dt_seconds = 0.05;
    int steps = 3;
    bool use_cfl_dt = true;
    double cfl_safety = 0.35;

    double epsilon_max = 0.01;
    bool include_gravity = true;
    bool include_pressure_gradient = true;
    bool rotating_frame_sources = true;

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
    bool write_spectrum = true;
    bool write_diagnostics = true;
    bool visualizer_enabled = true;
    int visualizer_every = 1;
    int visualizer_phi_index = 0;
    std::string visualizer_format = "ppm";
    // Comma-separated: temperature,velocity,density,pressure,opacity_proxy.
    std::string visualizer_modes = "temperature,velocity";

    double RInnerMeters() const { return r_inner_fraction * Constants::R_J; }
    double ROuterMeters() const { return r_outer_fraction * Constants::R_J; }

    DynamicsConfig Dynamics() const;
    RTConfig RT() const;
};

SimConfig LoadSimConfig(const std::string& path);
void PrintSimConfig(const SimConfig& cfg);

} // namespace Jupiter
