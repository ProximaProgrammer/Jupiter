#include "Diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace Jupiter
{

void Diagnostics::AppendCSV(const Grid& grid, const std::string& output_dir,
                            int step, double time_seconds, const TimeStepReport& dt_report)
{
    std::filesystem::create_directories(output_dir);
    const std::string path = output_dir + "/diagnostics.csv";
    const bool new_file = !std::filesystem::exists(path);

    std::ofstream out(path, std::ios::app);
    if (!out) throw std::runtime_error("Could not write diagnostics file: " + path);

    if (new_file)
    {
        out << "step,time_seconds,requested_dt,stable_dt,used_dt,cfl_max_speed_m_s,min_cell_length_m,"
               "mass_kg,total_energy_J,radial_momentum_kg_m_s,theta_momentum_kg_m_s,phi_momentum_kg_m_s,"
               "mean_temperature_K,min_temperature_K,max_temperature_K,"
               "mean_density_kg_m3,max_pressure_Pa,mean_speed_m_s,max_speed_m_s,max_mach,"
               "max_abs_vr_m_s,max_abs_vtheta_m_s,max_abs_vphi_m_s,mean_h3plus_mass_fraction\n";
    }

    double mass = 0.0;
    double energy = 0.0;
    double mr = 0.0, mt = 0.0, mp = 0.0;
    double sumT = 0.0;
    double sumRho = 0.0;
    double minT = std::numeric_limits<double>::infinity();
    double maxT = 0.0;
    double maxP = 0.0;
    double sumSpeed = 0.0;
    double maxSpeed = 0.0;
    double maxMach = 0.0;
    double maxVr = 0.0;
    double maxVtheta = 0.0;
    double maxVphi = 0.0;
    double sumH3 = 0.0;

    for (const Cell& c : grid.Cells())
    {
        mass += c.conserved.mass;
        energy += c.conserved.total_energy;
        mr += c.conserved.momentum_r;
        mt += c.conserved.momentum_theta;
        mp += c.conserved.momentum_phi;
        sumT += c.primitive.temperature;
        sumRho += c.primitive.density;
        minT = std::min(minT, c.primitive.temperature);
        maxT = std::max(maxT, c.primitive.temperature);
        maxP = std::max(maxP, c.primitive.pressure);

        const double vr = c.primitive.velocity_r;
        const double vt = c.primitive.velocity_theta;
        const double vp = c.primitive.velocity_phi;
        const double speed = std::sqrt(vr*vr + vt*vt + vp*vp);
        sumSpeed += speed;
        maxSpeed = std::max(maxSpeed, speed);
        maxMach = std::max(maxMach, speed / std::max(c.primitive.sound_speed, 1.0));
        maxVr = std::max(maxVr, std::abs(vr));
        maxVtheta = std::max(maxVtheta, std::abs(vt));
        maxVphi = std::max(maxVphi, std::abs(vp));
        sumH3 += c.primitive.mass_fraction[SIndex(Species::H3plus)];
    }

    const double n = static_cast<double>(grid.Size());
    out << step << ',' << time_seconds << ','
        << dt_report.requested_dt << ',' << dt_report.stable_dt << ',' << dt_report.used_dt << ','
        << dt_report.max_speed_m_s << ',' << dt_report.min_cell_length_m << ','
        << mass << ',' << energy << ',' << mr << ',' << mt << ',' << mp << ','
        << sumT / n << ',' << minT << ',' << maxT << ','
        << sumRho / n << ',' << maxP << ','
        << sumSpeed / n << ',' << maxSpeed << ',' << maxMach << ','
        << maxVr << ',' << maxVtheta << ',' << maxVphi << ',' << sumH3 / n << '\n';
}

} // namespace Jupiter
