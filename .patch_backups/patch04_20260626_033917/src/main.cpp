#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "common/Constants.hpp"
#include "config/SimConfig.hpp"
#include "core/Grid.hpp"
#include "cuda/CudaHooks.hpp"
#include "io/Diagnostics.hpp"
#include "io/RuntimeTable.hpp"
#include "io/Visualizer.hpp"
#include "physics/Dynamics.hpp"
#include "physics/RadiativeTransfer.hpp"
#include "physics/TimeStep.hpp"
#include "pipeline/InitialConditions.hpp"

using namespace Jupiter;

struct Args
{
    std::string config_path = "config/default.sim";
    bool cuda_smoke_test = false;

    bool override_nr = false;
    bool override_ntheta = false;
    bool override_nphi = false;
    bool override_steps = false;
    bool override_dt = false;
    bool override_assets = false;
    bool override_use_cuda_rt = false;
    bool override_visualizer_modes = false;

    int nr = 0;
    int ntheta = 0;
    int nphi = 0;
    int steps = 0;
    double dt = 0.0;
    std::string assets;
    bool use_cuda_rt = false;
    std::string visualizer_modes;
};

static Args parse_args(int argc, char** argv)
{
    Args a;
    for (int i = 1; i < argc; ++i)
    {
        const std::string k = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + name);
            return argv[++i];
        };
        if (k == "--config") a.config_path = need("--config");
        else if (k == "--nr") { a.nr = std::stoi(need("--nr")); a.override_nr = true; }
        else if (k == "--ntheta") { a.ntheta = std::stoi(need("--ntheta")); a.override_ntheta = true; }
        else if (k == "--nphi") { a.nphi = std::stoi(need("--nphi")); a.override_nphi = true; }
        else if (k == "--steps") { a.steps = std::stoi(need("--steps")); a.override_steps = true; }
        else if (k == "--dt") { a.dt = std::stod(need("--dt")); a.override_dt = true; }
        else if (k == "--assets") { a.assets = need("--assets"); a.override_assets = true; }
        else if (k == "--cuda-smoke-test") a.cuda_smoke_test = true;
        else if (k == "--use-cuda-rt") { a.use_cuda_rt = true; a.override_use_cuda_rt = true; }
        else if (k == "--viz") { a.visualizer_modes = need("--viz"); a.override_visualizer_modes = true; }
        else if (k == "--help")
        {
            std::cout << "Usage: jupiter [--config config/default.sim] [--nr N] [--ntheta N] [--nphi N]\n"
                         "               [--steps N] [--dt SEC] [--assets data/runtime]\n"
                         "               [--cuda-smoke-test] [--use-cuda-rt] [--viz temperature,velocity]\n";
            std::exit(0);
        }
        else throw std::runtime_error("Unknown argument: " + k);
    }
    return a;
}

static void apply_overrides(SimConfig& cfg, const Args& args)
{
    if (args.override_nr) cfg.nr = args.nr;
    if (args.override_ntheta) cfg.ntheta = args.ntheta;
    if (args.override_nphi) cfg.nphi = args.nphi;
    if (args.override_steps) cfg.steps = args.steps;
    if (args.override_dt) cfg.requested_dt_seconds = args.dt;
    if (args.override_assets) cfg.runtime_asset_dir = args.assets;
    if (args.override_use_cuda_rt) cfg.use_cuda_rt = args.use_cuda_rt;
    if (args.override_visualizer_modes) cfg.visualizer_modes = args.visualizer_modes;
}

static void write_spectrum_csv(const RayResult& ray, const std::string& path)
{
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream out(path);
    out << "wavelength_m,intensity_W_m3_sr\n";
    for (std::size_t i = 0; i < ray.wavelength_m.size(); ++i)
        out << ray.wavelength_m[i] << ',' << ray.intensity_W_m3_sr[i] << '\n';
}

int main(int argc, char** argv)
{
    try
    {
        const Args args = parse_args(argc, argv);
        SimConfig cfg = LoadSimConfig(args.config_path);
        apply_overrides(cfg, args);
        PrintSimConfig(cfg);

        if (args.cuda_smoke_test) RunCudaSmokeTest();

        RuntimeTables tables;
        if (cfg.use_runtime_tables && tables.Load(cfg.runtime_asset_dir))
            std::cout << "Loaded runtime spectral assets from " << cfg.runtime_asset_dir << " with "
                      << tables.Wavelengths().size() << " wavelengths.\n";
        else
            std::cout << "Runtime assets not loaded (" << tables.LastError()
                      << "). Using analytic fallback opacities/emission.\n";

        Grid grid(cfg.nr, cfg.ntheta, cfg.nphi, cfg.RInnerMeters(), cfg.ROuterMeters());
        InitialConditions::FillJupiterLike(grid, cfg.initial);
        Dynamics::UpdatePrimitives(grid);

        std::cout << "Cells = " << static_cast<long long>(grid.Size()) << "\n";
        std::cout << "CUDA built = " << (CudaBuilt() ? "yes" : "no") << "\n";

        double time_seconds = 0.0;
        TimeStepReport dt_report = TimeStep::ComputeCFL(grid, cfg.requested_dt_seconds, cfg.cfl_safety, cfg.use_cfl_dt);

        if (cfg.write_diagnostics) Diagnostics::AppendCSV(grid, cfg.output_dir, 0, time_seconds, dt_report);
        VisualizerConfig viz_cfg;
        viz_cfg.modes = cfg.visualizer_modes;
        viz_cfg.format = cfg.visualizer_format;
        viz_cfg.phi_index = cfg.visualizer_phi_index;

        if (cfg.visualizer_enabled) Visualizer::WriteMeridionalFrames(grid, cfg.output_dir, 0, viz_cfg);

        const DynamicsConfig dyn = cfg.Dynamics();
        for (int s = 1; s <= cfg.steps; ++s)
        {
            dt_report = TimeStep::ComputeCFL(grid, cfg.requested_dt_seconds, cfg.cfl_safety, cfg.use_cfl_dt);
            Dynamics::ApplySourceTerms(grid, dt_report.used_dt, dyn);
            Dynamics::ShiftEpsilon(grid, dt_report.used_dt, dyn);
            time_seconds += dt_report.used_dt;

            if (cfg.write_diagnostics) Diagnostics::AppendCSV(grid, cfg.output_dir, s, time_seconds, dt_report);
            if (cfg.visualizer_enabled && (s % cfg.visualizer_every == 0))
                Visualizer::WriteMeridionalFrames(grid, cfg.output_dir, s, viz_cfg);

            std::cout << "step " << s << " / " << cfg.steps
                      << ", dt=" << dt_report.used_dt
                      << " s, stable_dt=" << dt_report.stable_dt
                      << " s, t=" << time_seconds << " s\n";
        }

        if (cfg.write_spectrum)
        {
            RTConfig rt = cfg.RT();
            const int itheta = (cfg.ray_theta_index < 0) ? (cfg.ntheta / 2) : cfg.ray_theta_index;
            const int iphi = cfg.ray_phi_index;
            const RayResult ray = RadiativeTransfer::IntegrateRadialRay(grid, tables, itheta, iphi, rt);
            const std::string spectrum_path = cfg.output_dir + "/synthetic_radial_ray.csv";
            write_spectrum_csv(ray, spectrum_path);
            std::cout << "Wrote " << spectrum_path << " with " << ray.wavelength_m.size() << " bins.\n";
        }

        if (cfg.visualizer_enabled)
            std::cout << "Visualizer frames: " << cfg.output_dir << "/frames/*.ppm"
                      << " (modes: " << cfg.visualizer_modes << ")\n";
        if (cfg.write_diagnostics)
            std::cout << "Diagnostics: " << cfg.output_dir << "/diagnostics.csv\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
        return 1;
    }
}
