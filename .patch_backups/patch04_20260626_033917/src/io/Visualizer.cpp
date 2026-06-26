#include "Visualizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
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

std::string Trim(std::string s)
{
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

std::string Lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::vector<std::string> SplitModes(const std::string& raw)
{
    std::vector<std::string> modes;
    std::stringstream ss(raw);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        item = Lower(Trim(item));
        if (!item.empty()) modes.push_back(item);
    }
    if (modes.empty()) modes.push_back("temperature");
    return modes;
}

bool IsSupportedMode(const std::string& mode)
{
    return mode == "temperature" || mode == "temp" ||
           mode == "velocity" || mode == "velocity_magnitude" || mode == "speed" ||
           mode == "density" || mode == "rho" ||
           mode == "pressure" ||
           mode == "opacity_proxy" || mode == "opacity";
}

std::string CanonicalMode(std::string mode)
{
    mode = Lower(Trim(mode));
    if (mode == "temp") return "temperature";
    if (mode == "velocity_magnitude" || mode == "speed") return "velocity";
    if (mode == "rho") return "density";
    if (mode == "opacity") return "opacity_proxy";
    return mode;
}

RGB Heat(double q)
{
    q = std::clamp(q, 0.0, 1.0);
    auto b = [](double x) { return static_cast<std::uint8_t>(std::clamp(x, 0.0, 255.0)); };
    return {b(255.0 * q), b(255.0 * (1.0 - std::abs(2.0*q - 1.0))), b(255.0 * (1.0 - q))};
}

RGB SymmetricVelocityColor(double q)
{
    // q=0 cold, q=1 hot; keep same heat map for now so velocity images stay dependency-free.
    return Heat(q);
}

std::string FramePath(const std::string& output_dir, const std::string& mode, int step)
{
    std::ostringstream path;
    path << output_dir << "/frames/" << mode << "_meridional_"
         << std::setw(6) << std::setfill('0') << step << ".ppm";
    return path.str();
}

double VelocityMagnitude(const Cell& c)
{
    const auto& p = c.primitive;
    return std::sqrt(p.velocity_r*p.velocity_r + p.velocity_theta*p.velocity_theta + p.velocity_phi*p.velocity_phi);
}

double FieldValue(const Cell& c, const std::string& mode)
{
    if (mode == "temperature") return c.primitive.temperature;
    if (mode == "velocity") return VelocityMagnitude(c);
    if (mode == "density") return c.primitive.density;
    if (mode == "pressure") return c.primitive.pressure;

    if (mode == "opacity_proxy")
    {
        // Lightweight Rosseland-like proxy for visualization only.
        // Real spectral opacity still lives in RuntimeTables/RadiativeTransfer.
        const double rho = std::max(c.primitive.density, 1.0e-30);
        const double T = std::max(c.primitive.temperature, 1.0);
        return std::pow(rho / 0.16, 0.35) * std::pow(T / 165.0, -0.15);
    }

    return c.primitive.temperature;
}

bool UseLogScale(const std::string& mode)
{
    return mode == "temperature" || mode == "density" || mode == "pressure" || mode == "opacity_proxy";
}

double ScaledValue(double x, double lo, double hi, bool log_scale)
{
    if (log_scale)
    {
        const double llo = std::log(std::max(1.0e-300, lo));
        const double lhi = std::log(std::max(1.0e-299, hi));
        const double lx = std::log(std::max(1.0e-300, x));
        return (lx - llo) / std::max(1.0e-12, lhi - llo);
    }
    return (x - lo) / std::max(1.0e-12, hi - lo);
}

void WriteModePPM(const Grid& grid, const std::string& output_dir, int step, int phi_index, std::string mode)
{
    if (grid.NPhi() <= 0) return;
    mode = CanonicalMode(mode);
    if (!IsSupportedMode(mode)) throw std::runtime_error("Unsupported visualizer mode: " + mode);

    std::filesystem::create_directories(output_dir + "/frames");
    phi_index = ((phi_index % grid.NPhi()) + grid.NPhi()) % grid.NPhi();

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();

    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    {
        const double v = FieldValue(grid.At(ir, it, phi_index), mode);
        if (std::isfinite(v))
        {
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
    }

    if (!std::isfinite(vmin) || !std::isfinite(vmax))
    {
        vmin = 0.0;
        vmax = 1.0;
    }
    if (std::abs(vmax - vmin) < 1.0e-30) vmax = vmin + 1.0;

    const int width = grid.NTheta();
    const int height = grid.NR();

    std::ofstream out(FramePath(output_dir, mode, step), std::ios::binary);
    if (!out) throw std::runtime_error("Could not write visualizer PPM frame");
    out << "P6\n" << width << " " << height << "\n255\n";

    const bool log_scale = UseLogScale(mode);

    // Top row = outer radius, bottom row = deepest shell.
    for (int y = 0; y < height; ++y)
    {
        const int ir = grid.NR() - 1 - y;
        for (int x = 0; x < width; ++x)
        {
            const int it = x;
            const double v = FieldValue(grid.At(ir, it, phi_index), mode);
            const double q = ScaledValue(v, vmin, vmax, log_scale);
            const RGB rgb = (mode == "velocity") ? SymmetricVelocityColor(q) : Heat(q);
            out.put(static_cast<char>(rgb.r));
            out.put(static_cast<char>(rgb.g));
            out.put(static_cast<char>(rgb.b));
        }
    }

    // Sidecar file makes each autoscaled frame interpretable without adding dependencies.
    std::ofstream meta(FramePath(output_dir, mode, step) + ".txt");
    if (meta)
    {
        meta << "mode=" << mode << "\n"
             << "step=" << step << "\n"
             << "phi_index=" << phi_index << "\n"
             << "min=" << vmin << "\n"
             << "max=" << vmax << "\n"
             << "scale=" << (log_scale ? "log" : "linear") << "\n"
             << "width_theta=" << width << "\n"
             << "height_radius=" << height << "\n";
    }
}

} // namespace

std::vector<std::string> Visualizer::SupportedModes()
{
    return {"temperature", "velocity", "density", "pressure", "opacity_proxy"};
}

void Visualizer::WriteMeridionalTemperaturePPM(const Grid& grid, const std::string& output_dir,
                                               int step, int phi_index)
{
    WriteModePPM(grid, output_dir, step, phi_index, "temperature");
}

void Visualizer::WriteMeridionalFrames(const Grid& grid, const std::string& output_dir,
                                       int step, const VisualizerConfig& cfg)
{
    if (Lower(cfg.format) != "ppm")
        throw std::runtime_error("Only visualizer_format=ppm is currently supported.");

    for (const std::string& raw_mode : SplitModes(cfg.modes))
    {
        const std::string mode = CanonicalMode(raw_mode);
        if (!IsSupportedMode(mode))
            throw std::runtime_error("Unsupported visualizer mode: " + raw_mode);
        WriteModePPM(grid, output_dir, step, cfg.phi_index, mode);
    }
}

} // namespace Jupiter
