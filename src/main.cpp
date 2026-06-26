#include "config/SimConfig.h"
#include "core/Grid.h"
#include "io/Diagnostics.h"
#include "io/Visualizer.h"
#include "physics/Dynamics.h"
#include "physics/RadiativeTransfer.h"
#include "physics/TimeStep.h"

#include <exception>
#include <filesystem>
#include <iostream>

using namespace jrt;

int main(int argc, char** argv) {
    try {
        SimConfig cfg = SimConfig::fromArgs(argc, argv);
        std::filesystem::create_directories(cfg.output_dir);

        std::cout << "JupiterRT Patch 04\n"
                  << "Grid: nr=" << cfg.nr << " ntheta=" << cfg.ntheta << " nphi=" << cfg.nphi
                  << " ng=" << cfg.ng << "\n"
                  << "Rendering second-to-outermost radial shell; outermost radial layer is ghost.\n";

        Grid grid(cfg);
        grid.initializeJovianEnvelope();

        Diagnostics diag(cfg.output_dir);
        double time = 0.0;
        StepSummary initialSummary;
        initialSummary.dt = estimateStableDt(grid);
        diag.writeStep(grid, 0, time, initialSummary);
        const std::size_t jProfile = grid.jBegin(); // closest to theta = 0 without entering polar ghost
        const std::size_t kProfile = 0;             // closest to phi = 0
        diag.writeRadialProfile(grid, 0, jProfile, kProfile);
        diag.writeShellCsv(grid, 0);
        diag.writeShellPly(grid, 0);
        diag.writeVelocityObj(grid, 0);
        writeRadialSpectrumCsv(grid, cfg.output_dir + "/synthetic_radial_ray.csv", jProfile, kProfile);
        writeDashboard(grid, 0, time, cfg.output_dir);

        for (int step = 1; step <= cfg.steps; ++step) {
            const double dtRequest = estimateStableDt(grid);
            StepSummary summary = advanceDynamics(grid, dtRequest);
            time += summary.dt;

            if (step % cfg.output_every == 0 || step == cfg.steps) {
                diag.writeStep(grid, step, time, summary);
                diag.writeRadialProfile(grid, step, jProfile, kProfile);
                diag.writeShellCsv(grid, step);
                diag.writeShellPly(grid, step);
                diag.writeVelocityObj(grid, step);
                writeRadialSpectrumCsv(grid, cfg.output_dir + "/synthetic_radial_ray.csv", jProfile, kProfile);
                writeDashboard(grid, step, time, cfg.output_dir);
                std::cout << "step " << step << "/" << cfg.steps
                          << " t=" << time << " s"
                          << " dt=" << summary.dt
                          << " max|u|=" << summary.maxCorotatingSpeed
                          << " max|v_inertial|=" << summary.maxInertialSpeed << "\n";
            }
        }

        std::cout << "Done. Key outputs:\n"
                  << "  " << cfg.output_dir << "/diagnostics.csv\n"
                  << "  " << cfg.output_dir << "/profiles/radial_theta_phi_latest.csv\n"
                  << "  " << cfg.output_dir << "/frames/dashboard_latest.ppm\n"
                  << "  " << cfg.output_dir << "/shells/outer_interior_shell_latest.ply\n"
                  << "  " << cfg.output_dir << "/shells/velocity_shell_latest.obj\n"
                  << "  " << cfg.output_dir << "/synthetic_radial_ray.csv\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "JupiterRT error: " << e.what() << "\n";
        return 1;
    }
}
