#pragma once

#include "../core/Grid.hpp"

namespace Jupiter
{

struct TimeStepReport
{
    double requested_dt = 0.0;
    double stable_dt = 0.0;
    double used_dt = 0.0;
    double max_speed_m_s = 0.0;
    double min_cell_length_m = 0.0;
};

class TimeStep
{
public:
    static TimeStepReport ComputeCFL(const Grid& grid, double requested_dt, double cfl_safety, bool use_cfl_dt);
};

} // namespace Jupiter
