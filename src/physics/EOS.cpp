#include "physics/EOS.h"
#include "core/Grid.h"

#include <algorithm>
#include <cmath>

namespace jrt {

double specificGasConstant(double meanMolecularMassKg) {
    return kBoltzmann / meanMolecularMassKg;
}

double pressureIdealGas(double rho, double T, double meanMolecularMassKg) {
    return std::max(1.0e-30, rho) * specificGasConstant(meanMolecularMassKg) * std::max(1.0, T);
}

double densityIdealGas(double p, double T, double meanMolecularMassKg) {
    return std::max(1.0e-30, p) / (specificGasConstant(meanMolecularMassKg) * std::max(1.0, T));
}

double soundSpeedIdealGas(double gamma, double p, double rho) {
    return std::sqrt(std::max(1.0, gamma * std::max(1.0e-30, p) / std::max(1.0e-30, rho)));
}

double internalEnergyPerMass(double gamma, double p, double rho) {
    return std::max(0.0, p) / ((gamma - 1.0) * std::max(1.0e-30, rho));
}

double planckLambda(double wavelength_m, double temperature_k) {
    const double lambda = std::max(1.0e-12, wavelength_m);
    const double T = std::max(1.0, temperature_k);
    const double h = 6.62607015e-34;
    const double c = speedOfLight;
    const double k = kBoltzmann;
    const double x = h * c / (lambda * k * T);
    if (x > 700.0) return 0.0;
    return (2.0 * h * c * c) / std::pow(lambda, 5.0) / std::expm1(x);
}

double equilibriumTemperatureAtRadius(double r, double r_inner, double r_outer, double top_T, double deep_T) {
    const double eta = std::clamp((r - r_inner) / (r_outer - r_inner), 0.0, 1.0);
    return top_T + (deep_T - top_T) * std::pow(1.0 - eta, 1.2);
}

} // namespace jrt
