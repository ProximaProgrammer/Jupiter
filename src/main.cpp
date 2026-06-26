#include <fstream>
#include <iostream>
#include <string>

#include "common/Constants.hpp"
#include "core/Grid.hpp"
#include "cuda/CudaHooks.hpp"
#include "io/RuntimeTable.hpp"
#include "physics/Dynamics.hpp"
#include "physics/RadiativeTransfer.hpp"
#include "pipeline/InitialConditions.hpp"

using namespace Jupiter;

struct Args
{
    int nr = 32;
    int ntheta = 24;
    int nphi = 48;
    int steps = 1;
    double dt = 0.05;
    double r_inner = 0.20 * Constants::R_J;
    double r_outer = Constants::R_J;
    std::string assets = "data/runtime";
    bool cuda_smoke_test = false;
    bool use_cuda_rt = false;
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
        if (k == "--nr") a.nr = std::stoi(need("--nr"));
        else if (k == "--ntheta") a.ntheta = std::stoi(need("--ntheta"));
        else if (k == "--nphi") a.nphi = std::stoi(need("--nphi"));
        else if (k == "--steps") a.steps = std::stoi(need("--steps"));
        else if (k == "--dt") a.dt = std::stod(need("--dt"));
        else if (k == "--assets") a.assets = need("--assets");
        else if (k == "--cuda-smoke-test") a.cuda_smoke_test = true;
        else if (k == "--use-cuda-rt") a.use_cuda_rt = true;
        else if (k == "--help")
        {
            std::cout << "Usage: jupiter [--nr N] [--ntheta N] [--nphi N] [--steps N] [--dt SEC] [--assets data/runtime] [--cuda-smoke-test] [--use-cuda-rt]\n";
            std::exit(0);
        }
        else throw std::runtime_error("Unknown argument: " + k);
    }
    return a;
}

static void write_spectrum_csv(const RayResult& ray, const std::string& path)
{
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

        if (args.cuda_smoke_test) RunCudaSmokeTest();

        RuntimeTables tables;
        if (tables.Load(args.assets))
            std::cout << "Loaded runtime spectral assets from " << args.assets << " with "
                      << tables.Wavelengths().size() << " wavelengths.\n";
        else
            std::cout << "Runtime assets not loaded (" << tables.LastError()
                      << "). Using analytic fallback opacities/emission.\n";

        Grid grid(args.nr, args.ntheta, args.nphi, args.r_inner, args.r_outer);
        InitialConditions::FillJupiterLike(grid);
        Dynamics::UpdatePrimitives(grid);

        std::cout << "Cells = " << static_cast<long long>(grid.Size()) << "\n";
        std::cout << "CUDA built = " << (CudaBuilt() ? "yes" : "no") << "\n";

        DynamicsConfig dyn;
        dyn.epsilon_max = 0.01;
        for (int s = 0; s < args.steps; ++s)
        {
            Dynamics::ApplySourceTerms(grid, args.dt, dyn);
            Dynamics::ShiftEpsilon(grid, args.dt, dyn);
        }

        RTConfig rt;
        rt.include_doppler = true;
        rt.include_rotation_velocity = true;
        rt.use_cuda_rt_if_available = args.use_cuda_rt;
        const RayResult ray = RadiativeTransfer::IntegrateRadialRay(grid, tables, args.ntheta / 2, 0, rt);

        system("mkdir -p data/output");
        write_spectrum_csv(ray, "data/output/synthetic_radial_ray.csv");
        std::cout << "Wrote data/output/synthetic_radial_ray.csv with " << ray.wavelength_m.size() << " bins.\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
        return 1;
    }
}
