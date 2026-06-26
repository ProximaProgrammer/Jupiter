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
               "mass_kg,mean_temperature_K,min_temperature_K,max_temperature_K,"
               "mean_density_kg_m3,max_pressure_Pa,mean_speed_m_s,max_speed_m_s,"
               "max_abs_vr_m_s,max_abs_vtheta_m_s,max_abs_vphi_m_s\n";
    }

    double mass = 0.0;
    double sumT = 0.0;
    double sumRho = 0.0;
    double minT = std::numeric_limits<double>::infinity();
    double maxT = 0.0;
    double maxP = 0.0;
    double sumSpeed = 0.0;
    double maxSpeed = 0.0;
    double maxVr = 0.0;
    double maxVtheta = 0.0;
    double maxVphi = 0.0;

    for (const Cell& c : grid.Cells())
    {
        mass += c.conserved.mass;
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
        maxVr = std::max(maxVr, std::abs(vr));
        maxVtheta = std::max(maxVtheta, std::abs(vt));
        maxVphi = std::max(maxVphi, std::abs(vp));
    }

    const double n = static_cast<double>(grid.Size());
    out << step << ',' << time_seconds << ','
        << dt_report.requested_dt << ',' << dt_report.stable_dt << ',' << dt_report.used_dt << ','
        << dt_report.max_speed_m_s << ',' << dt_report.min_cell_length_m << ','
        << mass << ',' << sumT / n << ',' << minT << ',' << maxT << ','
        << sumRho / n << ',' << maxP << ','
        << sumSpeed / n << ',' << maxSpeed << ','
        << maxVr << ',' << maxVtheta << ',' << maxVphi << '\n';
}

} // namespace Jupiter
