#pragma once

#include "../core/Grid.hpp"

#include <string>
#include <vector>

namespace Jupiter
{

struct VisualizerConfig
{
    // Comma-separated modes from config/default.sim, for example:
    //   temperature,velocity,density,pressure
    // Supported modes: temperature, velocity, density, pressure, entropy, opacity_proxy.
    std::string modes = "temperature,velocity,density,pressure";
    std::string format = "ppm";
    int phi_index = 0;

    bool autoscale_each_frame = true;
    bool dashboard = true;
    int panel_width = 320;
    int panel_height = 480;
    int margin = 24;
    bool write_individual = false;
};

class Visualizer
{
public:
    // Backward-compatible helper retained for older call sites.
    static void WriteMeridionalTemperaturePPM(const Grid& grid, const std::string& output_dir,
                                              int step, int phi_index);

    // Writes configured individual mode frames and, by default, a single large
    // dashboard_NNNNNN.ppm containing all panels in one image.
    static void WriteMeridionalFrames(const Grid& grid, const std::string& output_dir,
                                      int step, const VisualizerConfig& cfg);

    static std::vector<std::string> SupportedModes();
};

} // namespace Jupiter
