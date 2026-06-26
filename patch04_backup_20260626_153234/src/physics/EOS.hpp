#pragma once

#include "Layer.hpp"
#include "../core/Cell.hpp"
#include "../state/ConservedState.hpp"
#include "../state/PrimitiveState.hpp"

namespace Jupiter
{

struct EOSResult
{
    PrimitiveState primitive;
    double internal_energy_density = 0.0;
};

class EOS
{
public:
    static EOSResult ConservedToPrimitive(const ConservedState& u, const CellGeometry& g);
    static ConservedState PrimitiveToConserved(const PrimitiveState& p, const CellGeometry& g);

    static double ColdPressure(double rho, const LayerParameters& lp);
    static double ElectronDegeneracyPressure(double rho, double mean_molar_mass);
    static double MeanMolarMassFromFractions(const ConservedState& u);
    static double PlanckLambda(double wavelength_m, double T);
};

} // namespace Jupiter
