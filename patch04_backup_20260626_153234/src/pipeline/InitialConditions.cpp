#include "InitialConditions.hpp"
#include "../common/Constants.hpp"
#include "../physics/EOS.hpp"
#include "../physics/Layer.hpp"

#include <algorithm>
#include <cmath>

namespace Jupiter
{

namespace
{

double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

// Deterministic, smooth angular seed; no random library needed and identical
// across machines. The radial envelope keeps perturbations away from hard
// boundaries while giving the molecular/supercritical transition something to do.
double seed_pattern(const CellGeometry& g, int modes)
{
    const double r = g.r_center / Constants::R_J;
    const double theta = g.theta_center;
    const double phi = g.phi_center;
    const double envelope = std::sin(Constants::pi * clamp01((r - 0.22) / 0.76));
    const double m = std::max(1, modes);
    const double a = std::sin((m + 1) * theta + 1.7 * std::cos(3.0 * phi));
    const double b = std::cos((m + 2) * phi - 0.5 * std::sin(5.0 * theta));
    const double c = std::sin(2.0 * Constants::pi * r * (0.5 * m + 1.0));
    return envelope * (0.55 * a * b + 0.45 * c * std::cos(phi));
}

} // namespace

void InitialConditions::FillJupiterLike(Grid& grid, const InitialConditionConfig& cfg)
{
    for (Cell& c : grid.Cells())
    {
        const double x = clamp01(c.geometry.r_center / Constants::R_J);
        const double q = 1.0 - x;
        const LayerKind layer = LayerForRadius(c.geometry.r_center);
        const LayerParameters lp = ParametersForLayer(layer);

        PrimitiveState p;
        p.density = cfg.surface_density_kg_m3 * std::pow(cfg.deep_density_kg_m3 / cfg.surface_density_kg_m3, std::pow(q, 1.42));
        p.temperature = cfg.surface_temperature_K + (cfg.deep_temperature_K - cfg.surface_temperature_K) * std::pow(q, 1.16);
        p.pressure = cfg.surface_pressure_Pa * std::pow(cfg.deep_pressure_Pa / cfg.surface_pressure_Pa, std::pow(q, 1.28));

        const double seed = seed_pattern(c.geometry, cfg.seed_modes);
        p.temperature *= std::max(0.2, 1.0 + cfg.thermal_perturbation_amplitude * seed);
        p.pressure *= std::max(0.2, 1.0 + 0.35 * cfg.thermal_perturbation_amplitude * seed);

        p.velocity_r = cfg.velocity_perturbation_m_s * 0.45 * seed;
        p.velocity_theta = cfg.velocity_perturbation_m_s * 0.25 * std::cos(c.geometry.theta_center * (cfg.seed_modes + 1)) * std::sin(2.0 * c.geometry.phi_center);
        p.velocity_phi = 0.0;
        p.mean_molar_mass = lp.mean_molar_mass;
        p.gamma_eff = lp.gamma;
        p.mass_fraction.fill(0.0);

        // Rough Jovian composition by mass, not number fraction. Trace absorber
        // fractions are intentionally tiny. The H3+ seed lives only high in the
        // atmosphere; the chemistry source term then keeps it time-dependent.
        p.mass_fraction[SIndex(Species::H2)] = 0.735;
        p.mass_fraction[SIndex(Species::He)] = 0.245;
        p.mass_fraction[SIndex(Species::CH4)] = 2.0e-3 * std::exp(-std::max(0.0, p.temperature - 1800.0) / 2200.0);
        p.mass_fraction[SIndex(Species::NH3)] = 4.0e-4 * std::exp(-std::max(0.0, p.temperature - 700.0) / 800.0);
        p.mass_fraction[SIndex(Species::PH3)] = 2.0e-6;
        p.mass_fraction[SIndex(Species::H2O)] = 1.0e-3;
        p.mass_fraction[SIndex(Species::H3plus)] = (x > 0.94) ? 1.0e-12 * (1.0 + 0.5 * std::max(0.0, seed)) : 0.0;

        // Zonal flow relative to the rotating grid. The explicit RT layer later
        // adds v_rot = Omega x r before applying Doppler shifts.
        const double lat = Constants::pi/2.0 - c.geometry.theta_center;
        p.velocity_phi += 70.0 * std::cos(lat) * std::exp(-std::pow((1.0 - x)/0.085, 2));
        p.velocity_phi += 18.0 * std::sin(6.0 * lat) * std::exp(-std::pow((1.0 - x)/0.13, 2));

        // Renormalize only the dominant background fractions so trace species do
        // not accidentally make the total mass fraction exceed one.
        double trace = 0.0;
        for (std::size_t i = SIndex(Species::CH4); i < NumSpecies; ++i) trace += p.mass_fraction[i];
        const double background = std::max(0.0, 1.0 - trace);
        const double h2he = 0.735 + 0.245;
        p.mass_fraction[SIndex(Species::H2)] = background * 0.735 / h2he;
        p.mass_fraction[SIndex(Species::He)] = background * 0.245 / h2he;

        c.conserved = EOS::PrimitiveToConserved(p, c.geometry);
        c.primitive = EOS::ConservedToPrimitive(c.conserved, c.geometry).primitive;
    }
}

} // namespace Jupiter
