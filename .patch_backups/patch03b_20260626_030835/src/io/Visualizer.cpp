#include "Visualizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace Jupiter
{

namespace
{

struct RGB { std::uint8_t r, g, b; };

RGB Heat(double q)
{
    q = std::clamp(q, 0.0, 1.0);
    auto b = [](double x) { return static_cast<std::uint8_t>(std::clamp(x, 0.0, 255.0)); };
    return {b(255.0 * q), b(255.0 * (1.0 - std::abs(2.0*q - 1.0))), b(255.0 * (1.0 - q))};
}

std::string FramePath(const std::string& output_dir, int step)
{
    std::ostringstream path;
    path << output_dir << "/frames/temp_meridional_" << std::setw(6) << std::setfill('0') << step << ".ppm";
    return path.str();
}

} // namespace

void Visualizer::WriteMeridionalTemperaturePPM(const Grid& grid, const std::string& output_dir,
                                               int step, int phi_index)
{
    std::filesystem::create_directories(output_dir + "/frames");
    if (grid.NPhi() <= 0) return;

    phi_index = ((phi_index % grid.NPhi()) + grid.NPhi()) % grid.NPhi();

    double tmin = std::numeric_limits<double>::infinity();
    double tmax = 0.0;
    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    {
        const double T = grid.At(ir, it, phi_index).primitive.temperature;
        tmin = std::min(tmin, T);
        tmax = std::max(tmax, T);
    }

    const double lo = std::log(std::max(1.0, tmin));
    const double hi = std::log(std::max(2.0, tmax));
    const double inv = 1.0 / std::max(1.0e-12, hi - lo);

    const int width = grid.NTheta();
    const int height = grid.NR();

    std::ofstream out(FramePath(output_dir, step), std::ios::binary);
    if (!out) throw std::runtime_error("Could not write visualizer PPM frame");
    out << "P6\n" << width << " " << height << "\n255\n";

    // Top row = outer radius, bottom row = deepest shell.
    for (int y = 0; y < height; ++y)
    {
        const int ir = grid.NR() - 1 - y;
        for (int x = 0; x < width; ++x)
        {
            const int it = x;
            const double T = grid.At(ir, it, phi_index).primitive.temperature;
            const double q = (std::log(std::max(1.0, T)) - lo) * inv;
            const RGB rgb = Heat(q);
            out.put(static_cast<char>(rgb.r));
            out.put(static_cast<char>(rgb.g));
            out.put(static_cast<char>(rgb.b));
        }
    }
}

} // namespace Jupiter
