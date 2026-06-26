#include "physics/RadiativeTransfer.h"
#include "physics/EOS.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace jrt {

namespace {

double gaussian(double x, double mu, double sig) {
    const double z = (x - mu) / sig;
    return std::exp(-0.5 * z * z);
}

} // namespace

double spectralOpacityM2PerKg(const Cell& c, double wavelength_m) {
    const double um = wavelength_m * 1.0e6;
    // Deliberate first physical approximation: grey CIA/cloud continuum plus
    // methane/H3+ structures in the JIRAM 2-5 micron range. Later patches can load
    // the data/runtime opacity tables generated from ExoMol/HITRAN-like assets.
    const double methane = c.x_ch4 * (0.9 * gaussian(um, 3.30, 0.16) + 0.25 * gaussian(um, 2.30, 0.20));
    const double h3p = 0.004 * c.h3p_emissivity * gaussian(um, 3.70, 0.35);
    const double window = 0.35 + 0.65 * (1.0 - gaussian(um, 4.95, 0.40));
    return std::max(1.0e-7, c.grey_opacity * window + methane + h3p);
}

std::vector<SpectrumSample> traceRadialRay(const Grid& grid, std::size_t j, std::size_t k, std::size_t bins) {
    std::vector<SpectrumSample> out;
    out.reserve(bins);
    const double lambda0 = 2.0e-6;
    const double lambda1 = 5.0e-6;

    for (std::size_t b = 0; b < bins; ++b) {
        const double q = (static_cast<double>(b) + 0.5) / static_cast<double>(bins);
        const double lambdaObs = lambda0 + q * (lambda1 - lambda0);
        double I = planckLambda(lambdaObs, grid.at(grid.iBegin(), j, k).T);
        double tauTotal = 0.0;

        for (std::size_t i = grid.iBegin(); i < grid.iEnd(); ++i) {
            const Cell& c = grid.at(i, j, k);
            // A radial ray uses u_r for the line-of-sight Doppler correction.
            // The general spacecraft line-of-sight will use the full inertial velocity.
            const double lambdaComoving = lambdaObs * (1.0 - c.u.r / speedOfLight);
            const double kappa = spectralOpacityM2PerKg(c, lambdaComoving);
            const double tau = std::max(0.0, kappa * c.rho * grid.dr);
            const double e = std::exp(-std::min(80.0, tau));
            const double S = planckLambda(lambdaComoving, c.T);
            I = I * e + S * (1.0 - e);
            tauTotal += tau;
        }
        out.push_back({lambdaObs, I, tauTotal});
    }
    return out;
}

void writeRadialSpectrumCsv(const Grid& grid, const std::string& path, std::size_t j, std::size_t k) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    const auto spec = traceRadialRay(grid, j, k, 336);
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Could not write " + path);
    f << "wavelength_m,intensity_W_m-3_sr-1,tau_total\n";
    f << std::setprecision(12);
    for (const auto& s : spec) {
        f << s.wavelength_m << ',' << s.intensity << ',' << s.tau_total << '\n';
    }
}

} // namespace jrt
