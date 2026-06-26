#include "io/Diagnostics.h"
#include "physics/EOS.h"
#include "physics/RadiativeTransfer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace jrt {
namespace {


std::string frameName(const std::string& dir, const std::string& prefix, int step, const std::string& ext) {
    std::ostringstream os;
    os << dir << "/" << prefix << '_' << std::setw(6) << std::setfill('0') << step << ext;
    return os.str();
}

std::array<unsigned char,3> temperatureColor(double T, double Tmin, double Tmax) {
    const double q = std::clamp((T - Tmin) / std::max(1.0, Tmax - Tmin), 0.0, 1.0);
    const auto toByte = [](double x) { return static_cast<unsigned char>(std::clamp(x, 0.0, 255.0)); };
    return {toByte(255.0 * q), toByte(255.0 * std::sqrt(q) * (1.0 - 0.35*q)), toByte(255.0 * (1.0 - q))};
}

} // namespace

Diagnostics::Diagnostics(std::string outputDir) : outDir(std::move(outputDir)) {
    std::filesystem::create_directories(outDir);
    std::filesystem::create_directories(outDir + "/profiles");
    std::filesystem::create_directories(outDir + "/shells");
}

void Diagnostics::writeStep(const Grid& grid, int step, double time_s, const StepSummary& summary) {
    const std::string path = outDir + "/diagnostics.csv";
    std::ofstream f(path, headerWritten ? std::ios::app : std::ios::out);
    if (!f) throw std::runtime_error("Could not write " + path);
    if (!headerWritten) {
        f << "step,time_s,dt_s,T_min_K,T_max_K,rho_min_kg_m-3,rho_max_kg_m-3,p_min_Pa,p_max_Pa,max_u_corot_m_s,max_v_inertial_m_s,max_signal_m_s,shell_index,outer_ghost_index\n";
        headerWritten = true;
    }

    double Tmin = std::numeric_limits<double>::infinity();
    double Tmax = -Tmin;
    double rmin = Tmin, rmax = -Tmin, pmin = Tmin, pmax = -Tmin;
    for (std::size_t i = grid.iBegin(); i < grid.iEnd(); ++i) {
        for (std::size_t j = grid.jBegin(); j < grid.jEnd(); ++j) {
            for (std::size_t k = grid.kBegin(); k < grid.kEnd(); ++k) {
                const Cell& c = grid.at(i,j,k);
                Tmin = std::min(Tmin, c.T); Tmax = std::max(Tmax, c.T);
                rmin = std::min(rmin, c.rho); rmax = std::max(rmax, c.rho);
                pmin = std::min(pmin, c.p); pmax = std::max(pmax, c.p);
            }
        }
    }
    f << std::setprecision(12)
      << step << ',' << time_s << ',' << summary.dt << ','
      << Tmin << ',' << Tmax << ',' << rmin << ',' << rmax << ',' << pmin << ',' << pmax << ','
      << summary.maxCorotatingSpeed << ',' << summary.maxInertialSpeed << ',' << summary.maxSignalSpeed << ','
      << grid.outerInteriorI() << ',' << grid.outerGhostI() << '\n';
}

void Diagnostics::writeRadialProfile(const Grid& grid, int step, std::size_t j, std::size_t k) {
    const std::string path = frameName(outDir + "/profiles", "radial_theta_phi", step, ".csv");
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Could not write " + path);
    f << "# radial profile closest to theta=0, phi=0; j=" << j << ", k=" << k << "\n";
    f << "i,r_m,r_over_RJ,theta_rad,phi_rad,T_K,rho_kg_m-3,p_Pa,u_r_m_s,u_theta_m_s,u_phi_corot_m_s,v_inertial_x_m_s,v_inertial_y_m_s,v_inertial_z_m_s,kappa_4p8um_m2_kg,tau_cell_4p8um\n";
    f << std::setprecision(12);
    for (std::size_t i = grid.iBegin(); i < grid.iEnd(); ++i) {
        const Cell& c = grid.at(i,j,k);
        const auto v = grid.inertialVelocityCartesian(i,j,k);
        const double kappa = spectralOpacityM2PerKg(c, 4.8e-6);
        f << i << ',' << grid.radius(i) << ',' << grid.radius(i)/grid.rOuter << ',' << grid.theta(j) << ',' << grid.phi(k) << ','
          << c.T << ',' << c.rho << ',' << c.p << ',' << c.u.r << ',' << c.u.th << ',' << c.u.ph << ','
          << v[0] << ',' << v[1] << ',' << v[2] << ',' << kappa << ',' << kappa * c.rho * grid.dr << '\n';
    }
    std::filesystem::copy_file(path, outDir + "/profiles/radial_theta_phi_latest.csv", std::filesystem::copy_options::overwrite_existing);
}

void Diagnostics::writeShellCsv(const Grid& grid, int step) {
    const std::size_t i = grid.outerInteriorI();
    const std::string path = frameName(outDir + "/shells", "outer_interior_shell", step, ".csv");
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Could not write " + path);
    f << "# second-to-outermost radial layer; outermost radial layer is ghost boundary\n";
    f << "i,j,k,x_m,y_m,z_m,lat_rad,lon_rad,T_K,rho_kg_m-3,p_Pa,u_r_m_s,u_theta_m_s,u_phi_corot_m_s,vx_inertial_m_s,vy_inertial_m_s,vz_inertial_m_s,kappa_4p8um_m2_kg\n";
    f << std::setprecision(12);
    for (std::size_t j = grid.jBegin(); j < grid.jEnd(); ++j) {
        for (std::size_t k = grid.kBegin(); k < grid.kEnd(); ++k) {
            const Cell& c = grid.at(i,j,k);
            const auto x = grid.cartesianPosition(i,j,k);
            const auto v = grid.inertialVelocityCartesian(i,j,k);
            f << i << ',' << j << ',' << k << ',' << x[0] << ',' << x[1] << ',' << x[2] << ','
              << grid.latitude(j) << ',' << grid.phi(k) << ',' << c.T << ',' << c.rho << ',' << c.p << ','
              << c.u.r << ',' << c.u.th << ',' << c.u.ph << ',' << v[0] << ',' << v[1] << ',' << v[2] << ','
              << spectralOpacityM2PerKg(c, 4.8e-6) << '\n';
        }
    }
    std::filesystem::copy_file(path, outDir + "/shells/outer_interior_shell_latest.csv", std::filesystem::copy_options::overwrite_existing);
}

void Diagnostics::writeShellPly(const Grid& grid, int step) {
    const std::size_t i = grid.outerInteriorI();
    const std::size_t nv = grid.cfg.ntheta * grid.cfg.nphi;
    const std::size_t nf = (grid.cfg.ntheta - 1) * grid.cfg.nphi;
    double Tmin = std::numeric_limits<double>::infinity();
    double Tmax = -Tmin;
    for (std::size_t j = grid.jBegin(); j < grid.jEnd(); ++j) {
        for (std::size_t k = grid.kBegin(); k < grid.kEnd(); ++k) {
            const double T = grid.at(i,j,k).T;
            Tmin = std::min(Tmin, T); Tmax = std::max(Tmax, T);
        }
    }

    const std::string path = frameName(outDir + "/shells", "outer_interior_shell", step, ".ply");
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Could not write " + path);
    f << "ply\nformat ascii 1.0\ncomment JupiterRT Patch 04 second-to-outermost shell; outer radial layer is ghost\n";
    f << "element vertex " << nv << "\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\n";
    f << "element face " << nf << "\nproperty list uchar int vertex_indices\nend_header\n";
    for (std::size_t jj = 0; jj < grid.cfg.ntheta; ++jj) {
        const std::size_t j = grid.jBegin() + jj;
        for (std::size_t k = grid.kBegin(); k < grid.kEnd(); ++k) {
            const auto x = grid.cartesianPosition(i,j,k);
            const auto color = temperatureColor(grid.at(i,j,k).T, Tmin, Tmax);
            f << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << static_cast<int>(color[0]) << ' ' << static_cast<int>(color[1]) << ' ' << static_cast<int>(color[2]) << '\n';
        }
    }
    for (std::size_t jj = 0; jj + 1 < grid.cfg.ntheta; ++jj) {
        for (std::size_t k = 0; k < grid.cfg.nphi; ++k) {
            const std::size_t k2 = (k + 1) % grid.cfg.nphi;
            const std::size_t a = jj * grid.cfg.nphi + k;
            const std::size_t b = jj * grid.cfg.nphi + k2;
            const std::size_t c = (jj + 1) * grid.cfg.nphi + k2;
            const std::size_t d = (jj + 1) * grid.cfg.nphi + k;
            f << "4 " << a << ' ' << b << ' ' << c << ' ' << d << '\n';
        }
    }
    std::filesystem::copy_file(path, outDir + "/shells/outer_interior_shell_latest.ply", std::filesystem::copy_options::overwrite_existing);
}

void Diagnostics::writeVelocityObj(const Grid& grid, int step, std::size_t strideTheta, std::size_t stridePhi) {
    const std::size_t i = grid.outerInteriorI();
    const std::string path = frameName(outDir + "/shells", "velocity_shell", step, ".obj");
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Could not write " + path);
    f << "# JupiterRT Patch 04 velocity arrows on second-to-outermost shell\n";
    std::size_t vertex = 1;
    const double arrowScale = 6000.0; // m per (m/s), purely visual
    for (std::size_t j = grid.jBegin(); j < grid.jEnd(); j += std::max<std::size_t>(1, strideTheta)) {
        for (std::size_t k = grid.kBegin(); k < grid.kEnd(); k += std::max<std::size_t>(1, stridePhi)) {
            const auto x = grid.cartesianPosition(i,j,k);
            const auto v = grid.inertialVelocityCartesian(i,j,k);
            f << "v " << x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
            f << "v " << x[0] + arrowScale * v[0] << ' ' << x[1] + arrowScale * v[1] << ' ' << x[2] + arrowScale * v[2] << '\n';
            f << "l " << vertex << ' ' << vertex + 1 << '\n';
            vertex += 2;
        }
    }
    std::filesystem::copy_file(path, outDir + "/shells/velocity_shell_latest.obj", std::filesystem::copy_options::overwrite_existing);
}

} // namespace jrt
