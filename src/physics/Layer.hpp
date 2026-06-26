#pragma once

#include <string_view>
#include "../common/Constants.hpp"

namespace Jupiter
{

enum class LayerKind
{
    UpperAtmosphere,
    MolecularHydrogen,
    SupercriticalHydrogen,
    MetallicHydrogen
};

struct LayerParameters
{
    LayerKind kind = LayerKind::UpperAtmosphere;
    double gamma = 1.4;
    double mean_molar_mass = 2.30e-3;
    double rho0 = 1.0;
    double bulk_modulus = 0.0;
    double cold_gamma = 2.0;
    double gruneisen = 1.0;
    bool electron_degenerate = false;
};

inline LayerKind LayerForRadius(double r)
{
    const double x = r / Constants::R_J;
    if (x > 0.96) return LayerKind::UpperAtmosphere;
    if (x > 0.78) return LayerKind::MolecularHydrogen;
    if (x > 0.55) return LayerKind::SupercriticalHydrogen;
    return LayerKind::MetallicHydrogen;
}

inline std::string_view LayerName(LayerKind k)
{
    switch (k)
    {
        case LayerKind::UpperAtmosphere: return "upper_atmosphere";
        case LayerKind::MolecularHydrogen: return "molecular_hydrogen";
        case LayerKind::SupercriticalHydrogen: return "supercritical_hydrogen";
        case LayerKind::MetallicHydrogen: return "metallic_hydrogen";
    }
    return "unknown";
}

inline LayerParameters ParametersForLayer(LayerKind k)
{
    switch (k)
    {
        case LayerKind::UpperAtmosphere:
            // H2/He ideal gas with trace molecules.
            return {k, 1.40, 2.30e-3, 0.16, 0.0, 2.0, 1.0, false};
        case LayerKind::MolecularHydrogen:
            // Dense molecular H2/He: still mostly molecular, but gamma softened by rotational/vibrational modes.
            return {k, 1.34, 2.27e-3, 20.0, 0.0, 2.0, 1.0, false};
        case LayerKind::SupercriticalHydrogen:
            // Analytic Mie-Gruneisen-like closure. Replace with table later through the same EOS interface.
            return {k, 1.25, 2.15e-3, 350.0, 2.0e11, 2.35, 1.25, false};
        case LayerKind::MetallicHydrogen:
            // Partially ionized metallic H/He: cold degeneracy + thermal ions/electrons.
            return {k, 1.18, 1.45e-3, 900.0, 8.0e11, 2.60, 1.50, true};
    }
    return {};
}

} // namespace Jupiter
