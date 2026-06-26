#pragma once

#include "../core/Grid.hpp"

#include <string>

namespace Jupiter
{

class Visualizer
{
public:
    static void WriteMeridionalTemperaturePPM(const Grid& grid, const std::string& output_dir,
                                              int step, int phi_index);
};

} // namespace Jupiter
