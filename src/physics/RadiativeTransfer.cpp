#include "RadiativeTransfer.hpp"
#include "EOS.hpp"
#include "../common/Constants.hpp"
#include "../cuda/CudaHooks.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace Jupiter
{

static double number_density_from_mass_fraction(const Cell& cell, Species s)
{
    const double X = cell.primitive.mass_fraction[SIndex(s)];
    return cell.primitive.density * X / SpeciesMolarMassKgPerMol[SIndex(s)] * 6.02214076e23;
}

static double pressure_bar(const Cell& cell)
{
    return cell.primitive.pressure * 1.0e-5;
}

static double analytic_fallback_sigma(const char* species, double wavelength_m, double T)
{
    const double um = wavelength_m * 1.0e6;
    auto band = [&](double c, double w, double amp)
    {
        const double x = (std::log(um / c) / w);
        return amp * std::exp(-0.5 * x * x);
    };
    double sigma = 0.0;
    if (std::string(species) == "CH4") sigma += band(3.3, 0.08, 1.0e-25) + band(7.7, 0.10, 3.0e-26);
    if (std::string(species) == "NH3") sigma += band(10.5, 0.12, 5.0e-26);
    if (std::string(species) == "PH3") sigma += band(4.3, 0.10, 2.0e-26);
    if (std::string(species) == "H2O") sigma += band(2.7, 0.12, 2.0e-26) + band(6.3, 0.12, 2.0e-26);
    // Smoothly weaken fallback at extreme temperatures so it does not dominate deep layers.
    return sigma * std::clamp(1200.0 / std::max(T, 100.0), 0.05, 10.0);
}

double RadiativeTransfer::ComovingWavelength(double lab_wavelength_m,
                                             const Cell& cell,
                                             const Vec3& photon_direction_cart,
                                             const RTConfig& cfg)
{
    if (!cfg.include_doppler) return lab_wavelength_m;
    const auto& g = cell.geometry;
    const auto& p = cell.primitive;
    Vec3 v = SphericalVelocityToCartesian(g.theta_center, g.phi_center,
                                          p.velocity_r, p.velocity_theta, p.velocity_phi);
    if (cfg.include_rotation_velocity)
    {
        const Vec3 omega{0.0, 0.0, Constants::Omega_J};
        const Vec3 r = PositionCartesian(g.r_center, g.theta_center, g.phi_center);
        v += Cross(omega, r); // v_rot = Omega_J x r.
    }
    const double beta = std::clamp(Dot(Unit(photon_direction_cart), v) / Constants::c, -0.2, 0.2);
    const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    return lab_wavelength_m / (gamma * (1.0 - beta));
}

double RadiativeTransfer::CellAbsorptionPerMeter(const Cell& cell,
                                                 const RuntimeTables& tables,
                                                 double wavelength_m,
                                                 const RTConfig& cfg)
{
    const double T = cell.primitive.temperature;
    const double pbar = pressure_bar(cell);
    double alpha = 0.0;

    const auto add_species = [&](Species s, const char* name)
    {
        const double n = number_density_from_mass_fraction(cell, s);
        double sigma = 0.0;
        if (cfg.use_runtime_tables && !tables.Empty()) sigma = tables.InterpolateOpacity(name, pbar, T, wavelength_m);
        if (sigma <= 0.0) sigma = analytic_fallback_sigma(name, wavelength_m, T);
        alpha += n * sigma;
    };

    add_species(Species::CH4, "CH4");
    add_species(Species::NH3, "NH3");
    add_species(Species::PH3, "PH3");
    add_species(Species::H2O, "H2O");

    // Collision-induced absorption continuum. If runtime CIA tables exist, use them as
    // alpha = k_cia n_A n_B with k_cia in m^5 molecule^-2. Otherwise use a smooth fallback.
    const double nH2 = number_density_from_mass_fraction(cell, Species::H2);
    const double nHe = number_density_from_mass_fraction(cell, Species::He);
    double cia_alpha = 0.0;
    if (cfg.use_runtime_tables && !tables.Empty())
    {
        cia_alpha += tables.InterpolateCIA("H2-H2", T, wavelength_m) * nH2 * nH2;
        cia_alpha += tables.InterpolateCIA("H2-He", T, wavelength_m) * nH2 * nHe;
    }
    if (cia_alpha <= 0.0)
    {
        const double um = wavelength_m * 1e6;
        const double cia_shape = std::exp(-0.5 * std::pow(std::log(um / 2.2) / 0.35, 2)) +
                                 0.35 * std::exp(-0.5 * std::pow(std::log(um / 4.0) / 0.45, 2));
        cia_alpha = 1.0e-56 * (nH2*nH2 + 0.25*nH2*nHe) * cia_shape;
    }
    alpha += cia_alpha;

    return std::max(0.0, cfg.opacity_scale * alpha);
}

double RadiativeTransfer::CellSourceFunction(const Cell& cell,
                                             const RuntimeTables& tables,
                                             double wavelength_m,
                                             const RTConfig& cfg)
{
    const double B = EOS::PlanckLambda(wavelength_m, cell.primitive.temperature);
    double extra = 0.0;
    const double nH3 = number_density_from_mass_fraction(cell, Species::H3plus);
    if (nH3 > 0.0)
    {
        const double h3_shape = (!tables.Empty()) ? tables.InterpolateH3Emission(cell.primitive.temperature, wavelength_m) : 0.0;
        if (h3_shape > 0.0) extra += cfg.h3plus_emission_scale * nH3 * h3_shape;
    }
    return B + extra;
}

RayResult RadiativeTransfer::IntegrateRadialRay(const Grid& grid,
                                                const RuntimeTables& tables,
                                                int itheta,
                                                int iphi,
                                                const RTConfig& cfg)
{
    RayResult rr;
    if (!tables.Empty()) rr.wavelength_m = tables.Wavelengths();
    else
    {
        const int n = 1024;
        rr.wavelength_m.resize(n);
        const double lo = std::log(0.8e-6), hi = std::log(50.0e-6);
        for (int i = 0; i < n; ++i) rr.wavelength_m[i] = std::exp(lo + (hi-lo) * i / (n-1));
    }
    rr.intensity_W_m3_sr.assign(rr.wavelength_m.size(), 0.0);

    itheta = std::clamp(itheta, 0, grid.NTheta() - 1);
    iphi = ((iphi % grid.NPhi()) + grid.NPhi()) % grid.NPhi();

    // Outward photon direction for this radial column.
    const Cell& surf = grid.At(grid.NR()-1, itheta, iphi);
    const Vec3 n = SphericalBasisR(surf.geometry.theta_center, surf.geometry.phi_center);

    std::vector<double> alpha_layer(rr.wavelength_m.size());
    std::vector<double> source_layer(rr.wavelength_m.size());

    for (int ir = 0; ir < grid.NR(); ++ir)
    {
        const Cell& cell = grid.At(ir, itheta, iphi);
        const double ds = cell.geometry.length_r;

        for (std::size_t l = 0; l < rr.wavelength_m.size(); ++l)
        {
            const double lambda_lab = rr.wavelength_m[l];
            const double lambda_cell = ComovingWavelength(lambda_lab, cell, n, cfg);
            alpha_layer[l] = CellAbsorptionPerMeter(cell, tables, lambda_cell, cfg);
            source_layer[l] = CellSourceFunction(cell, tables, lambda_cell, cfg);
        }

        // On CPU builds this is a CPU fallback. On RTX builds with --use-cuda-rt it runs the spectral
        // layer update in one CUDA kernel launch. Table interpolation stays CPU-side for now, which keeps
        // the code simple while making the spectral recurrence GPU-ready.
        if (cfg.use_cuda_rt_if_available && CudaBuilt())
        {
            CudaApplyRTLayer(rr.intensity_W_m3_sr, alpha_layer, source_layer, ds);
        }
        else
        {
            for (std::size_t l = 0; l < rr.wavelength_m.size(); ++l)
            {
                const double tau = std::clamp(alpha_layer[l] * ds, 0.0, 700.0);
                const double atten = std::exp(-tau);
                rr.intensity_W_m3_sr[l] = rr.intensity_W_m3_sr[l] * atten + source_layer[l] * (1.0 - atten);
            }
        }
    }
    return rr;
}

} // namespace Jupiter
