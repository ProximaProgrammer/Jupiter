#pragma once

#include <cstddef>
#include <string>

namespace jrt {

struct SimConfig {
    // Grid: one radial ghost layer is reserved at each boundary by default.
    // The shell rendered/exported by Patch 04 is the second-to-outermost radial layer:
    // shell_i = ng + nr - 1. The outermost layer shell_i + 1 is the ghost boundary.
    std::size_t nr = 48;
    std::size_t ntheta = 36;
    std::size_t nphi = 72;
    std::size_t ng = 1;

    // This patch evolves Jupiter's observable outer envelope rather than the full
    // metallic interior. Later patches can extend r_inner_fraction downward.
    double r_inner_fraction = 0.96;
    double r_outer_m = 71492.0e3;       // equatorial radius scale [m]
    double gm_m3_s2 = 1.26686534e17;   // Jupiter GM [m^3 s^-2]
    double omega_rad_s = 1.758518e-4;  // System III-ish angular rate [rad s^-1]

    // H/He atmosphere defaults. EOS.cpp keeps the conversion explicit.
    double gamma = 1.40;
    double mean_molecular_weight_kg = 2.30 * 1.66053906660e-27;

    // Stable default: CFL-limited timestep is computed every step and capped here.
    double cfl = 0.18;
    double max_dt_s = 25.0;
    double min_dt_s = 0.25;
    int steps = 180;
    int output_every = 10;

    // Crude turbulent closures; these are deliberately explicit and visible.
    double kinematic_viscosity_m2_s = 2.0e6;
    double thermal_diffusivity_m2_s = 2.5e6;
    double velocity_drag_time_s = 9.0e5;
    double radiative_time_top_s = 5.0e4;
    double radiative_time_deep_s = 3.0e6;

    // Initial state controls.
    double top_temperature_k = 128.0;
    double deep_temperature_k = 450.0;
    double top_pressure_pa = 7.0e4;
    double jet_speed_m_s = 75.0;
    double perturbation_speed_m_s = 8.0;

    std::string output_dir = "data/output";

    static SimConfig fromArgs(int argc, char** argv);
};

} // namespace jrt
