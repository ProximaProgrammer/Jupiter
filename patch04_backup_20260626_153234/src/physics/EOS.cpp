#include "EOS.hpp"
#include "../common/Constants.hpp"

#include <algorithm>
#include <cmath>

namespace Jupiter
{

static double sqr(double x) { return x * x; }

static double safe_temperature(double p_thermal, double rho, double mu)
{
    if (rho <= 0.0) return 2.7;
    return std::clamp(p_thermal * mu / (rho * Constants::Rgas), 2.7, 1.0e6);
}

double EOS::MeanMolarMassFromFractions(const ConservedState& u)
{
    if (u.mass <= 0.0) return 2.30e-3;
    double denom = 0.0;
    for (std::size_t i = 0; i < NumSpecies; ++i)
    {
        const double X = std::max(0.0, u.species_mass[i] / u.mass);
        denom += X / SpeciesMolarMassKgPerMol[i];
    }
    return denom > 0.0 ? 1.0 / denom : 2.30e-3;
}

double EOS::ColdPressure(double rho, const LayerParameters& lp)
{
    if (lp.bulk_modulus <= 0.0 || rho <= lp.rho0) return 0.0;
    const double x = rho / lp.rho0;
    return lp.bulk_modulus / lp.cold_gamma * (std::pow(x, lp.cold_gamma) - 1.0);
}

double EOS::ElectronDegeneracyPressure(double rho, double mean_molar_mass)
{
    if (rho <= 0.0) return 0.0;
    // Non-relativistic electron degeneracy: P = (hbar^2/5m_e)(3pi^2)^(2/3) n_e^(5/3).
    // mu_e is estimated from H/He mixture. This is not a full ionization EOS; it supplies the cold metallic floor.
    const double mu_e = std::max(1.0, mean_molar_mass / (1.0e-3));
    const double n_e = rho / (mu_e * Constants::m_u);
    const double prefactor = (Constants::hbar * Constants::hbar / (5.0 * Constants::m_e)) * std::pow(3.0 * Constants::pi * Constants::pi, 2.0 / 3.0);
    return prefactor * std::pow(n_e, 5.0 / 3.0);
}

EOSResult EOS::ConservedToPrimitive(const ConservedState& u, const CellGeometry& g)
{
    EOSResult out;
    PrimitiveState& p = out.primitive;
    if (g.volume <= 0.0 || u.mass <= 0.0)
    {
        p.temperature = 2.7;
        return out;
    }

    const double rho = std::max(1.0e-12, u.mass / g.volume);
    p.density = rho;
    p.velocity_r = u.momentum_r / std::max(u.mass, 1.0e-300);
    p.velocity_theta = u.momentum_theta / std::max(u.mass, 1.0e-300);
    p.velocity_phi = u.momentum_phi / std::max(u.mass, 1.0e-300);
    p.magnetic_r = u.magnetic_r;
    p.magnetic_theta = u.magnetic_theta;
    p.magnetic_phi = u.magnetic_phi;

    for (std::size_t i = 0; i < NumSpecies; ++i)
        p.mass_fraction[i] = std::max(0.0, u.species_mass[i] / std::max(u.mass, 1.0e-300));

    const LayerParameters lp = ParametersForLayer(LayerForRadius(g.r_center));
    p.gamma_eff = lp.gamma;
    p.mean_molar_mass = MeanMolarMassFromFractions(u);
    if (!(p.mean_molar_mass > 0.0) || !std::isfinite(p.mean_molar_mass)) p.mean_molar_mass = lp.mean_molar_mass;

    const double kinetic = 0.5 * u.mass * (sqr(p.velocity_r) + sqr(p.velocity_theta) + sqr(p.velocity_phi));
    const double magnetic = 0.5 * (sqr(u.magnetic_r) + sqr(u.magnetic_theta) + sqr(u.magnetic_phi)) * g.volume;
    const double e_int_total = std::max(1.0e-30, u.total_energy - kinetic - magnetic);
    const double e_int_density = e_int_total / g.volume;
    out.internal_energy_density = e_int_density;

    p.cold_pressure = ColdPressure(rho, lp);
    p.electron_degeneracy_pressure = lp.electron_degenerate ? ElectronDegeneracyPressure(rho, p.mean_molar_mass) : 0.0;
    const double p_thermal = std::max(1.0e-20, (lp.gamma - 1.0) * e_int_density);
    p.pressure = p_thermal + p.cold_pressure + p.electron_degeneracy_pressure;
    p.temperature = safe_temperature(p_thermal, rho, p.mean_molar_mass);
    p.sound_speed = std::sqrt(std::max(0.0, lp.gamma * (p.pressure + p.cold_pressure) / rho));
    return out;
}

ConservedState EOS::PrimitiveToConserved(const PrimitiveState& p, const CellGeometry& g)
{
    ConservedState u;
    u.mass = std::max(0.0, p.density * g.volume);
    u.momentum_r = u.mass * p.velocity_r;
    u.momentum_theta = u.mass * p.velocity_theta;
    u.momentum_phi = u.mass * p.velocity_phi;
    u.magnetic_r = p.magnetic_r;
    u.magnetic_theta = p.magnetic_theta;
    u.magnetic_phi = p.magnetic_phi;

    const LayerParameters lp = ParametersForLayer(LayerForRadius(g.r_center));
    const double cold = ColdPressure(p.density, lp);
    const double degen = lp.electron_degenerate ? ElectronDegeneracyPressure(p.density, p.mean_molar_mass) : 0.0;
    const double p_thermal = std::max(1.0e-20, p.pressure - cold - degen);
    const double e_int_density = p_thermal / std::max(1.0e-12, lp.gamma - 1.0);
    const double kinetic = 0.5 * u.mass * (sqr(p.velocity_r) + sqr(p.velocity_theta) + sqr(p.velocity_phi));
    const double magnetic = 0.5 * (sqr(p.magnetic_r) + sqr(p.magnetic_theta) + sqr(p.magnetic_phi)) * g.volume;
    u.total_energy = e_int_density * g.volume + kinetic + magnetic;
    for (std::size_t i = 0; i < NumSpecies; ++i)
        u.species_mass[i] = u.mass * std::max(0.0, p.mass_fraction[i]);
    return u;
}

double EOS::PlanckLambda(double wavelength_m, double T)
{
    const double lam = std::max(wavelength_m, 1.0e-20);
    const double x = Constants::h * Constants::c / (lam * Constants::kB * std::max(T, 2.7));
    if (x > 700.0) return 0.0;
    const double denom = std::expm1(x);
    if (denom <= 0.0) return 0.0;
    return 2.0 * Constants::h * Constants::c * Constants::c / (std::pow(lam, 5) * denom);
}

} // namespace Jupiter
