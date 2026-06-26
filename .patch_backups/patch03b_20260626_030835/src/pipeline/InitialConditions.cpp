#include "InitialConditions.hpp"
#include "../common/Constants.hpp"
#include "../physics/EOS.hpp"

#include <algorithm>
#include <cmath>

namespace Jupiter
{

void InitialConditions::FillJupiterLike(Grid& grid, const InitialConditionConfig& cfg)
{
    for (Cell& c : grid.Cells())
    {
        const double x = std::clamp(c.geometry.r_center / Constants::R_J, 0.0, 1.0);
        const double q = 1.0 - x;

        PrimitiveState p;
        p.density = cfg.surface_density_kg_m3 * std::pow(cfg.deep_density_kg_m3 / cfg.surface_density_kg_m3, std::pow(q, 1.4));
        p.temperature = cfg.surface_temperature_K + (cfg.deep_temperature_K - cfg.surface_temperature_K) * std::pow(q, 1.15);
        p.pressure = cfg.surface_pressure_Pa * std::pow(cfg.deep_pressure_Pa / cfg.surface_pressure_Pa, std::pow(q, 1.25));
        p.velocity_r = 0.0;
        p.velocity_theta = 0.0;
        p.velocity_phi = 0.0;
        p.mean_molar_mass = 2.30e-3;
        p.mass_fraction.fill(0.0);

        // Rough Jovian composition by mass, not number fraction. Trace absorber fractions are intentionally tiny.
        p.mass_fraction[SIndex(Species::H2)] = 0.735;
        p.mass_fraction[SIndex(Species::He)] = 0.245;
        p.mass_fraction[SIndex(Species::CH4)] = 2.0e-3;
        p.mass_fraction[SIndex(Species::NH3)] = 4.0e-4;
        p.mass_fraction[SIndex(Species::PH3)] = 2.0e-6;
        p.mass_fraction[SIndex(Species::H2O)] = 1.0e-3;
        p.mass_fraction[SIndex(Species::H3plus)] = (x > 0.94) ? 1.0e-12 : 0.0;

        // Add a small zonal flow perturbation relative to the rotating grid.
        const double lat = Constants::pi/2.0 - c.geometry.theta_center;
        p.velocity_phi = 60.0 * std::cos(lat) * std::exp(-std::pow((1.0 - x)/0.08, 2));

        c.conserved = EOS::PrimitiveToConserved(p, c.geometry);
        c.primitive = EOS::ConservedToPrimitive(c.conserved, c.geometry).primitive;
    }
}

} // namespace Jupiter
