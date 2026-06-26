#pragma once

namespace jrt {

double specificGasConstant(double meanMolecularMassKg);
double pressureIdealGas(double rho, double T, double meanMolecularMassKg);
double densityIdealGas(double p, double T, double meanMolecularMassKg);
double soundSpeedIdealGas(double gamma, double p, double rho);
double internalEnergyPerMass(double gamma, double p, double rho);

double planckLambda(double wavelength_m, double temperature_k);
double equilibriumTemperatureAtRadius(double r, double r_inner, double r_outer, double top_T, double deep_T);

} // namespace jrt
