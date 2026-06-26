#pragma once

#include "core/Grid.h"

#include <string>

namespace jrt {

void writeDashboard(const Grid& grid, int step, double time_s, const std::string& outputDir);

} // namespace jrt
