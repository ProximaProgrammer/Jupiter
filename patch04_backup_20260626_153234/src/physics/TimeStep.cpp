#include "TimeStep.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Jupiter
{

TimeStepReport TimeStep::ComputeCFL(const Grid& grid, double requested_dt, double cfl_safety, bool use_cfl_dt)
{
    TimeStepReport report;
    report.requested_dt = requested_dt;
    report.stable_dt = requested_dt;
    report.used_dt = requested_dt;
    report.max_speed_m_s = 0.0;
    report.min_cell_length_m = std::numeric_limits<double>::infinity();

    double candidate = std::numeric_limits<double>::infinity();

    for (const Cell& c : grid.Cells())
    {
        const auto& p = c.primitive;
        const auto& g = c.geometry;
        const double speed = std::max(1.0,
            std::abs(p.velocity_r) + std::abs(p.velocity_theta) + std::abs(p.velocity_phi) + std::max(0.0, p.sound_speed));
        const double length = std::max(1.0, std::min({g.length_r, g.length_theta, g.length_phi}));

        report.max_speed_m_s = std::max(report.max_speed_m_s, speed);
        report.min_cell_length_m = std::min(report.min_cell_length_m, length);
        candidate = std::min(candidate, cfl_safety * length / speed);
    }

    if (!std::isfinite(candidate)) candidate = requested_dt;
    report.stable_dt = candidate;
    report.used_dt = use_cfl_dt ? std::min(requested_dt, candidate) : requested_dt;
    return report;
}

} // namespace Jupiter
