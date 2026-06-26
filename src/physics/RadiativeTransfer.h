#pragma once

#include "core/Grid.h"

#include <cstddef>
#include <string>
#include <vector>

namespace jrt {

struct SpectrumSample {
    double wavelength_m = 0.0;
    double intensity = 0.0;
    double tau_total = 0.0;
};

double spectralOpacityM2PerKg(const Cell& c, double wavelength_m);
std::vector<SpectrumSample> traceRadialRay(const Grid& grid, std::size_t j, std::size_t k, std::size_t bins = 336);
void writeRadialSpectrumCsv(const Grid& grid, const std::string& path, std::size_t j, std::size_t k);

} // namespace jrt
