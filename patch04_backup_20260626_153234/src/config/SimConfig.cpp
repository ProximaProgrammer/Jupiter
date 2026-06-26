#include "SimConfig.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace Jupiter
{

namespace
{

std::string Trim(std::string s)
{
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

std::string Lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

bool ParseBool(const std::string& raw)
{
    const std::string s = Lower(Trim(raw));
    return s == "1" || s == "true" || s == "yes" || s == "on";
}

} // namespace

DynamicsConfig SimConfig::Dynamics() const
{
    DynamicsConfig cfg;
    cfg.epsilon_max = epsilon_max;
    cfg.use_finite_volume_advection = use_finite_volume_advection;
    cfg.include_gravity = include_gravity;
    cfg.include_pressure_gradient = include_pressure_gradient;
    cfg.rotating_frame_sources = rotating_frame_sources;
    cfg.include_viscosity = include_viscosity;
    cfg.include_thermal_diffusion = include_thermal_diffusion;
    cfg.include_convective_buoyancy = include_convective_buoyancy;
    cfg.include_radiative_relaxation = include_radiative_relaxation;
    cfg.include_chemistry = include_chemistry;
    cfg.enforce_boundaries = enforce_boundaries;
    cfg.kinematic_viscosity_m2_s = kinematic_viscosity_m2_s;
    cfg.thermal_diffusivity_m2_s = thermal_diffusivity_m2_s;
    cfg.convective_buoyancy_strength = convective_buoyancy_strength;
    cfg.velocity_drag_timescale_s = velocity_drag_timescale_s;
    cfg.top_radiative_timescale_s = top_radiative_timescale_s;
    cfg.deep_radiative_timescale_s = deep_radiative_timescale_s;
    cfg.chemistry_timescale_s = chemistry_timescale_s;
    return cfg;
}

RTConfig SimConfig::RT() const
{
    RTConfig cfg;
    cfg.include_doppler = include_doppler;
    cfg.include_rotation_velocity = include_rotation_velocity;
    cfg.use_runtime_tables = use_runtime_tables;
    cfg.use_cuda_rt_if_available = use_cuda_rt;
    cfg.opacity_scale = opacity_scale;
    cfg.h3plus_emission_scale = h3plus_emission_scale;
    return cfg;
}

SimConfig LoadSimConfig(const std::string& path)
{
    SimConfig cfg;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open config file: " + path);

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line))
    {
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = Trim(line);
        if (line.empty()) continue;

        const std::size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = Trim(line.substr(0, eq));
        const std::string val = Trim(line.substr(eq + 1));
        if (!key.empty()) kv[key] = val;
    }

    auto getI = [&](const char* k, int& x) { if (kv.count(k)) x = std::stoi(kv.at(k)); };
    auto getD = [&](const char* k, double& x) { if (kv.count(k)) x = std::stod(kv.at(k)); };
    auto getB = [&](const char* k, bool& x) { if (kv.count(k)) x = ParseBool(kv.at(k)); };
    auto getS = [&](const char* k, std::string& x) { if (kv.count(k)) x = kv.at(k); };

    getI("nr", cfg.nr);
    getI("ntheta", cfg.ntheta);
    getI("nphi", cfg.nphi);
    getD("r_inner_fraction", cfg.r_inner_fraction);
    getD("r_outer_fraction", cfg.r_outer_fraction);

    getD("requested_dt_seconds", cfg.requested_dt_seconds);
    getI("steps", cfg.steps);
    getB("use_cfl_dt", cfg.use_cfl_dt);
    getD("cfl_safety", cfg.cfl_safety);

    getD("epsilon_max", cfg.epsilon_max);
    getB("use_finite_volume_advection", cfg.use_finite_volume_advection);
    getB("include_gravity", cfg.include_gravity);
    getB("include_pressure_gradient", cfg.include_pressure_gradient);
    getB("rotating_frame_sources", cfg.rotating_frame_sources);
    getB("include_viscosity", cfg.include_viscosity);
    getB("include_thermal_diffusion", cfg.include_thermal_diffusion);
    getB("include_convective_buoyancy", cfg.include_convective_buoyancy);
    getB("include_radiative_relaxation", cfg.include_radiative_relaxation);
    getB("include_chemistry", cfg.include_chemistry);
    getB("enforce_boundaries", cfg.enforce_boundaries);
    getD("kinematic_viscosity_m2_s", cfg.kinematic_viscosity_m2_s);
    getD("thermal_diffusivity_m2_s", cfg.thermal_diffusivity_m2_s);
    getD("convective_buoyancy_strength", cfg.convective_buoyancy_strength);
    getD("velocity_drag_timescale_s", cfg.velocity_drag_timescale_s);
    getD("top_radiative_timescale_s", cfg.top_radiative_timescale_s);
    getD("deep_radiative_timescale_s", cfg.deep_radiative_timescale_s);
    getD("chemistry_timescale_s", cfg.chemistry_timescale_s);

    getD("surface_temperature_K", cfg.initial.surface_temperature_K);
    getD("deep_temperature_K", cfg.initial.deep_temperature_K);
    getD("surface_pressure_Pa", cfg.initial.surface_pressure_Pa);
    getD("deep_pressure_Pa", cfg.initial.deep_pressure_Pa);
    getD("surface_density_kg_m3", cfg.initial.surface_density_kg_m3);
    getD("deep_density_kg_m3", cfg.initial.deep_density_kg_m3);
    getD("initial_thermal_perturbation", cfg.initial.thermal_perturbation_amplitude);
    getD("initial_velocity_perturbation_m_s", cfg.initial.velocity_perturbation_m_s);
    getI("initial_seed_modes", cfg.initial.seed_modes);

    getS("raw_exomol_dir", cfg.raw_exomol_dir);
    getS("raw_cia_dir", cfg.raw_cia_dir);
    getS("runtime_asset_dir", cfg.runtime_asset_dir);
    getB("use_runtime_tables", cfg.use_runtime_tables);

    getB("include_doppler", cfg.include_doppler);
    getB("include_rotation_velocity", cfg.include_rotation_velocity);
    getB("use_cuda_rt", cfg.use_cuda_rt);
    getD("opacity_scale", cfg.opacity_scale);
    getD("h3plus_emission_scale", cfg.h3plus_emission_scale);
    getI("ray_theta_index", cfg.ray_theta_index);
    getI("ray_phi_index", cfg.ray_phi_index);

    getS("output_dir", cfg.output_dir);
    getB("clean_output_on_start", cfg.clean_output_on_start);
    getB("write_spectrum", cfg.write_spectrum);
    getB("write_diagnostics", cfg.write_diagnostics);
    getB("visualizer_enabled", cfg.visualizer_enabled);
    getI("visualizer_every", cfg.visualizer_every);
    getI("visualizer_phi_index", cfg.visualizer_phi_index);
    getS("visualizer_format", cfg.visualizer_format);
    getS("visualizer_modes", cfg.visualizer_modes);
    getB("visualizer_dashboard", cfg.visualizer_dashboard);
    getI("visualizer_panel_width", cfg.visualizer_panel_width);
    getI("visualizer_panel_height", cfg.visualizer_panel_height);
    getI("visualizer_margin", cfg.visualizer_margin);
    getB("visualizer_write_individual", cfg.visualizer_write_individual);
    getB("visualizer_autoscale_each_frame", cfg.visualizer_autoscale_each_frame);

    if (cfg.nr < 4 || cfg.ntheta < 2 || cfg.nphi < 1) throw std::runtime_error("Config grid resolution is too small.");
    if (cfg.steps < 0) throw std::runtime_error("Config steps must be nonnegative.");
    if (!(cfg.r_outer_fraction > cfg.r_inner_fraction && cfg.r_inner_fraction >= 0.0)) throw std::runtime_error("Invalid radial fractions.");
    if (cfg.visualizer_every <= 0) cfg.visualizer_every = 1;
    if (cfg.visualizer_panel_width < 64) cfg.visualizer_panel_width = 64;
    if (cfg.visualizer_panel_height < 64) cfg.visualizer_panel_height = 64;
    if (cfg.visualizer_margin < 0) cfg.visualizer_margin = 0;

    return cfg;
}

void PrintSimConfig(const SimConfig& cfg)
{
    const long long cells = static_cast<long long>(cfg.nr) * cfg.ntheta * cfg.nphi;
    std::cout << "JupiterRT config\n"
              << "  grid: " << cfg.nr << " x " << cfg.ntheta << " x " << cfg.nphi
              << " = " << cells << " cells\n"
              << "  radius: " << cfg.r_inner_fraction << ".." << cfg.r_outer_fraction << " R_J\n"
              << "  requested dt: " << cfg.requested_dt_seconds << " s"
              << (cfg.use_cfl_dt ? " with CFL cap" : " without CFL cap") << "\n"
              << "  steps: " << cfg.steps << "\n"
              << "  physics: FV=" << (cfg.use_finite_volume_advection ? "on" : "off")
              << ", convection=" << (cfg.include_convective_buoyancy ? "on" : "off")
              << ", diffusion=" << (cfg.include_thermal_diffusion ? "on" : "off")
              << ", viscosity=" << (cfg.include_viscosity ? "on" : "off")
              << ", RT relax=" << (cfg.include_radiative_relaxation ? "on" : "off") << "\n"
              << "  runtime assets: " << cfg.runtime_asset_dir << "\n"
              << "  output: " << cfg.output_dir << "\n"
              << "  visualizer modes: " << cfg.visualizer_modes
              << (cfg.visualizer_dashboard ? " plus dashboard" : "") << "\n";
}

} // namespace Jupiter
