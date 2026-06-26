#include "config/SimConfig.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace jrt {
namespace {

bool hasValue(int i, int argc) { return i + 1 < argc; }

double parseDouble(const char* s, const std::string& name) {
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (!end || *end != '\0') throw std::runtime_error("Bad value for " + name + ": " + s);
    return v;
}

int parseInt(const char* s, const std::string& name) {
    char* end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (!end || *end != '\0') throw std::runtime_error("Bad value for " + name + ": " + s);
    return static_cast<int>(v);
}

} // namespace

SimConfig SimConfig::fromArgs(int argc, char** argv) {
    SimConfig cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const std::string& name) {
            if (!hasValue(i, argc)) throw std::runtime_error("Missing value after " + name);
        };
        if (a == "--nr") { need(a); cfg.nr = static_cast<std::size_t>(parseInt(argv[++i], a)); }
        else if (a == "--ntheta") { need(a); cfg.ntheta = static_cast<std::size_t>(parseInt(argv[++i], a)); }
        else if (a == "--nphi") { need(a); cfg.nphi = static_cast<std::size_t>(parseInt(argv[++i], a)); }
        else if (a == "--ng") { need(a); cfg.ng = static_cast<std::size_t>(parseInt(argv[++i], a)); }
        else if (a == "--steps") { need(a); cfg.steps = parseInt(argv[++i], a); }
        else if (a == "--output-every") { need(a); cfg.output_every = parseInt(argv[++i], a); }
        else if (a == "--max-dt") { need(a); cfg.max_dt_s = parseDouble(argv[++i], a); }
        else if (a == "--cfl") { need(a); cfg.cfl = parseDouble(argv[++i], a); }
        else if (a == "--output-dir") { need(a); cfg.output_dir = argv[++i]; }
        else if (a == "--r-inner-fraction") { need(a); cfg.r_inner_fraction = parseDouble(argv[++i], a); }
        else if (a == "--jet-speed") { need(a); cfg.jet_speed_m_s = parseDouble(argv[++i], a); }
        else if (a == "--help" || a == "-h") {
            std::cout << "JupiterRT Patch 04\n"
                      << "  --nr N --ntheta N --nphi N      grid size excluding ghost cells\n"
                      << "  --steps N                       number of matter timesteps\n"
                      << "  --output-every N                frame/diagnostic cadence\n"
                      << "  --max-dt SECONDS                cap the CFL timestep\n"
                      << "  --output-dir PATH               output directory\n"
                      << "  --r-inner-fraction X            inner radius as fraction of R_J\n"
                      << "  --jet-speed MPS                 initialized alternating jet amplitude\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + a);
        }
    }

    if (cfg.nr < 6 || cfg.ntheta < 6 || cfg.nphi < 8) {
        throw std::runtime_error("Grid too small; use at least --nr 6 --ntheta 6 --nphi 8.");
    }
    if (cfg.ng < 1) throw std::runtime_error("Patch 04 requires at least one ghost cell: --ng 1.");
    if (cfg.output_every < 1) cfg.output_every = 1;
    if (cfg.steps < 0) cfg.steps = 0;
    return cfg;
}

} // namespace jrt
