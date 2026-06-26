#pragma once

#include "../core/Grid.hpp"

#include <string>
#include <vector>

namespace Jupiter
{

struct VisualizerConfig
{
    // Comma-separated modes from config/default.sim, for example:
    //   temperature,velocity
    // Supported modes: temperature, velocity, density, pressure, opacity_proxy.
    std::string modes = "temperature,velocity";
    std::string format = "ppm";
    int phi_index = 0;

    // If true, each frame uses its own min/max. Good for debugging subtle evolution.
    // If false, current implementation still uses per-frame autoscale until fixed ranges are added.
    bool autoscale_each_frame = true;
};

class Visualizer
{
public:
    // Backward-compatible helper retained for older call sites.
    static void WriteMeridionalTemperaturePPM(const Grid& grid, const std::string& output_dir,
                                              int step, int phi_index);

    // Primary patch-03 entry point. Writes all configured meridional slice frames.
    static void WriteMeridionalFrames(const Grid& grid, const std::string& output_dir,
                                      int step, const VisualizerConfig& cfg);

    static std::vector<std::string> SupportedModes();
};

} // namespace Jupiter
