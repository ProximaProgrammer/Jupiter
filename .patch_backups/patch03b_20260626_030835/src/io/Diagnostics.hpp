#pragma once

#include "../core/Grid.hpp"
#include "../physics/TimeStep.hpp"

#include <string>

namespace Jupiter
{

class Diagnostics
{
public:
    static void AppendCSV(const Grid& grid, const std::string& output_dir,
                          int step, double time_seconds, const TimeStepReport& dt_report);
};

} // namespace Jupiter
