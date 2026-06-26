#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace Jupiter
{

struct SpeciesOpacityTable
{
    std::string name;
    std::vector<double> pressures_bar;
    std::vector<double> temperatures_K;
    // sigma[p][T][lambda], flattened ((p*nT + t)*nLambda + l), m^2 molecule^-1.
    std::vector<float> sigma_m2_per_molecule;
};

struct CIATable
{
    std::string pair_name; // e.g. H2-H2 or H2-He
    std::vector<double> temperatures_K;
    // k_cia[T][lambda], m^5 molecule^-2. Absorption alpha = k_cia n_A n_B.
    std::vector<float> coeff_m5_per_molecule2;
};

struct EmissionTable
{
    std::string name;
    std::vector<double> temperatures_K;
    // F[T][lambda], source-native H3+ emissivity shape.
    std::vector<float> F_lambda;
};

class RuntimeTables
{
public:
    bool Load(const std::string& asset_dir);
    bool Empty() const { return wavelengths_m.empty(); }

    double InterpolateOpacity(const std::string& species, double pressure_bar, double temperature_K, double wavelength_m) const;
    double InterpolateCIA(const std::string& pair_name, double temperature_K, double wavelength_m) const;
    double InterpolateH3Emission(double temperature_K, double wavelength_m) const;

    const std::vector<double>& Wavelengths() const { return wavelengths_m; }
    const std::string& LastError() const { return last_error_; }
    std::size_t OpacityTableCount() const { return opacities_.size(); }
    std::size_t CIATableCount() const { return cia_.size(); }

private:
    std::vector<double> wavelengths_m;
    std::vector<SpeciesOpacityTable> opacities_;
    std::vector<CIATable> cia_;
    EmissionTable h3plus_;
    std::string last_error_;

    std::size_t LowerIndex(const std::vector<double>& axis, double x) const;
    double InterpLambda(const std::vector<float>& values, std::size_t stride_offset, std::size_t nLambda, double wavelength_m) const;
};

} // namespace Jupiter
