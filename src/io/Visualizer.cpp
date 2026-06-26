#include "Visualizer.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Jupiter
{

namespace
{

struct RGB
{
    std::uint8_t r = 0;
    std::uint8_t g = 0;
    std::uint8_t b = 0;
};

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

std::vector<std::string> SplitModes(const std::string& csv)
{
    std::vector<std::string> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        item = Trim(item);
        if (!item.empty()) out.push_back(item);
    }
    if (out.empty()) out.push_back("temperature");
    return out;
}

bool IsSupportedMode(const std::string& mode)
{
    return mode == "temperature" || mode == "temp" ||
           mode == "velocity" || mode == "velocity_magnitude" || mode == "speed" ||
           mode == "density" || mode == "rho" ||
           mode == "pressure" ||
           mode == "entropy" || mode == "log_entropy" ||
           mode == "opacity_proxy" || mode == "opacity";
}

std::string CanonicalMode(std::string mode)
{
    mode = Lower(Trim(mode));
    if (mode == "temp") return "temperature";
    if (mode == "velocity_magnitude" || mode == "speed") return "velocity";
    if (mode == "rho") return "density";
    if (mode == "log_entropy") return "entropy";
    if (mode == "opacity") return "opacity_proxy";
    return mode;
}

RGB Heat(double q)
{
    q = std::clamp(q, 0.0, 1.0);
    auto b = [](double x) { return static_cast<std::uint8_t>(std::clamp(x, 0.0, 255.0)); };
    // blue -> cyan/green -> yellow/red.
    return {b(255.0 * q), b(255.0 * (1.0 - std::abs(2.0*q - 1.0))), b(255.0 * (1.0 - q))};
}

RGB VelocityColor(double q)
{
    // Keep a compatible dependency-free palette. Velocity panels use linear scale
    // so changing motion is visible even when absolute temperatures dwarf it.
    return Heat(q);
}

RGB Background()
{
    return {24, 24, 24};
}

RGB Border()
{
    return {86, 86, 86};
}

std::string FramePath(const std::string& output_dir, const std::string& mode, int step)
{
    std::ostringstream path;
    path << output_dir << "/frames/" << mode << "_meridional_"
         << std::setw(6) << std::setfill('0') << step << ".ppm";
    return path.str();
}

std::string DashboardPath(const std::string& output_dir, int step)
{
    std::ostringstream path;
    path << output_dir << "/frames/dashboard_"
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
    if (mode == "entropy")
    {
        const double P = std::max(c.primitive.pressure, 1.0e-300);
        const double rho = std::max(c.primitive.density, 1.0e-300);
        return std::log(P) - c.primitive.gamma_eff * std::log(rho);
    }
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

struct Range
{
    double lo = 0.0;
    double hi = 1.0;
};

Range ComputeRange(const Grid& grid, int phi_index, const std::string& mode)
{
    Range r;
    r.lo = std::numeric_limits<double>::infinity();
    r.hi = -std::numeric_limits<double>::infinity();
    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    {
        const double v = FieldValue(grid.At(ir, it, phi_index), mode);
        if (std::isfinite(v))
        {
            r.lo = std::min(r.lo, v);
            r.hi = std::max(r.hi, v);
        }
    }
    if (!std::isfinite(r.lo) || !std::isfinite(r.hi)) r = {0.0, 1.0};
    if (std::abs(r.hi - r.lo) < 1.0e-30) r.hi = r.lo + 1.0;
    return r;
}

void WritePPM(const std::string& path, const std::vector<RGB>& pixels, int width, int height)
{
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not write visualizer PPM frame: " + path);
    out << "P6\n" << width << " " << height << "\n255\n";
    for (const RGB& p : pixels)
    {
        out.put(static_cast<char>(p.r));
        out.put(static_cast<char>(p.g));
        out.put(static_cast<char>(p.b));
    }
}

RGB SamplePanel(const Grid& grid, int phi_index, const std::string& mode, const Range& range,
                int px, int py, int width, int height)
{
    const double fx = (width <= 1) ? 0.0 : static_cast<double>(px) / static_cast<double>(width - 1);
    const double fy = (height <= 1) ? 0.0 : static_cast<double>(py) / static_cast<double>(height - 1);
    const int it = std::clamp(static_cast<int>(fx * grid.NTheta()), 0, grid.NTheta() - 1);
    // Top row = outer radius, bottom row = deepest shell.
    const int ir = std::clamp(grid.NR() - 1 - static_cast<int>(fy * grid.NR()), 0, grid.NR() - 1);
    const double v = FieldValue(grid.At(ir, it, phi_index), mode);
    const double q = ScaledValue(v, range.lo, range.hi, UseLogScale(mode));
    return (mode == "velocity") ? VelocityColor(q) : Heat(q);
}

void DrawPanel(std::vector<RGB>& canvas, int canvas_w, int canvas_h,
               const Grid& grid, int phi_index, const std::string& mode, const Range& range,
               int x0, int y0, int panel_w, int panel_h)
{
    (void)canvas_h;
    for (int y = 0; y < panel_h; ++y)
    for (int x = 0; x < panel_w; ++x)
    {
        canvas[static_cast<std::size_t>(y0 + y) * canvas_w + (x0 + x)] =
            SamplePanel(grid, phi_index, mode, range, x, y, panel_w, panel_h);
    }

    // Thin neutral border and a vertical color bar. Text is omitted to keep the
    // writer dependency-free; the sidecar metadata names every panel and scale.
    const RGB border = Border();
    for (int x = -1; x <= panel_w; ++x)
    {
        const int xx = x0 + x;
        if (xx >= 0 && xx < canvas_w)
        {
            if (y0 - 1 >= 0) canvas[static_cast<std::size_t>(y0 - 1) * canvas_w + xx] = border;
            canvas[static_cast<std::size_t>(y0 + panel_h) * canvas_w + xx] = border;
        }
    }
    for (int y = -1; y <= panel_h; ++y)
    {
        const int yy = y0 + y;
        if (yy >= 0 && yy < canvas_h)
        {
            if (x0 - 1 >= 0) canvas[static_cast<std::size_t>(yy) * canvas_w + (x0 - 1)] = border;
            canvas[static_cast<std::size_t>(yy) * canvas_w + (x0 + panel_w)] = border;
        }
    }
}

void WriteModePPM(const Grid& grid, const std::string& output_dir, int step, int phi_index, std::string mode)
{
    if (grid.NPhi() <= 0) return;
    mode = CanonicalMode(mode);
    if (!IsSupportedMode(mode)) throw std::runtime_error("Unsupported visualizer mode: " + mode);

    std::filesystem::create_directories(output_dir + "/frames");
    phi_index = ((phi_index % grid.NPhi()) + grid.NPhi()) % grid.NPhi();

    const int width = std::max(64, grid.NTheta());
    const int height = std::max(64, grid.NR());
    const Range range = ComputeRange(grid, phi_index, mode);
    std::vector<RGB> pixels(static_cast<std::size_t>(width) * height);
    for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x)
        pixels[static_cast<std::size_t>(y) * width + x] = SamplePanel(grid, phi_index, mode, range, x, y, width, height);

    WritePPM(FramePath(output_dir, mode, step), pixels, width, height);

    std::ofstream meta(FramePath(output_dir, mode, step) + ".txt");
    if (meta)
    {
        meta << "mode=" << mode << "\n"
             << "step=" << step << "\n"
             << "phi_index=" << phi_index << "\n"
             << "min=" << range.lo << "\n"
             << "max=" << range.hi << "\n"
             << "scale=" << (UseLogScale(mode) ? "log" : "linear") << "\n"
             << "width_theta=" << width << "\n"
             << "height_radius=" << height << "\n";
    }
}

void WriteDashboardPPM(const Grid& grid, const std::string& output_dir, int step, int phi_index, const VisualizerConfig& cfg)
{
    std::vector<std::string> modes;
    for (const std::string& raw : SplitModes(cfg.modes))
    {
        const std::string mode = CanonicalMode(raw);
        if (!IsSupportedMode(mode)) throw std::runtime_error("Unsupported visualizer mode: " + raw);
        modes.push_back(mode);
    }

    phi_index = ((phi_index % grid.NPhi()) + grid.NPhi()) % grid.NPhi();
    const int panel_w = std::max(64, cfg.panel_width);
    const int panel_h = std::max(64, cfg.panel_height);
    const int margin = std::max(0, cfg.margin);
    const int n = static_cast<int>(modes.size());
    const int width = margin + n * panel_w + (n - 1) * margin + margin;
    const int height = panel_h + 2 * margin;
    std::vector<RGB> canvas(static_cast<std::size_t>(width) * height, Background());

    std::vector<Range> ranges;
    for (const std::string& mode : modes) ranges.push_back(ComputeRange(grid, phi_index, mode));

    int x = margin;
    for (int i = 0; i < n; ++i)
    {
        DrawPanel(canvas, width, height, grid, phi_index, modes[i], ranges[i], x, margin, panel_w, panel_h);
        x += panel_w + margin;
    }

    const std::string path = DashboardPath(output_dir, step);
    WritePPM(path, canvas, width, height);

    // Keep a stable filename that Preview/browser refresh workflows can watch.
    try
    {
        std::filesystem::copy_file(path, output_dir + "/frames/dashboard_latest.ppm",
                                   std::filesystem::copy_options::overwrite_existing);
    }
    catch (...) {}

    std::ofstream meta(path + ".txt");
    if (meta)
    {
        meta << "step=" << step << "\n"
             << "phi_index=" << phi_index << "\n"
             << "dashboard_width=" << width << "\n"
             << "dashboard_height=" << height << "\n";
        for (std::size_t i = 0; i < modes.size(); ++i)
        {
            meta << "panel_" << i << "=" << modes[i]
                 << ",min=" << ranges[i].lo
                 << ",max=" << ranges[i].hi
                 << ",scale=" << (UseLogScale(modes[i]) ? "log" : "linear") << "\n";
        }
    }
}

} // namespace

std::vector<std::string> Visualizer::SupportedModes()
{
    return {"temperature", "velocity", "density", "pressure", "entropy", "opacity_proxy"};
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

    std::filesystem::create_directories(output_dir + "/frames");

    if (cfg.write_individual)
    {
        for (const std::string& raw_mode : SplitModes(cfg.modes))
        {
            const std::string mode = CanonicalMode(raw_mode);
            if (!IsSupportedMode(mode))
                throw std::runtime_error("Unsupported visualizer mode: " + raw_mode);
            WriteModePPM(grid, output_dir, step, cfg.phi_index, mode);
        }
    }

    if (cfg.dashboard)
        WriteDashboardPPM(grid, output_dir, step, cfg.phi_index, cfg);
}

} // namespace Jupiter
