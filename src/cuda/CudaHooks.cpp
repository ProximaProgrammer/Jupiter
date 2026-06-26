#include "CudaHooks.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace Jupiter
{

#ifndef JUPITER_ENABLE_CUDA
bool CudaBuilt() { return false; }

void RunCudaSmokeTest()
{
    std::cout << "CUDA was not enabled in this build. Reconfigure on the RTX host with -DJUPITER_ENABLE_CUDA=ON.\n";
}

void CudaApplyRTLayer(std::vector<double>& intensity,
                      const std::vector<double>& alpha_per_m,
                      const std::vector<double>& source_function,
                      double ds_m)
{
    if (intensity.size() != alpha_per_m.size() || intensity.size() != source_function.size())
        throw std::runtime_error("CudaApplyRTLayer CPU fallback received mismatched vector sizes");

    for (std::size_t i = 0; i < intensity.size(); ++i)
    {
        const double tau = std::clamp(alpha_per_m[i] * ds_m, 0.0, 700.0);
        const double atten = std::exp(-tau);
        intensity[i] = intensity[i] * atten + source_function[i] * (1.0 - atten);
    }
}
#endif

} // namespace Jupiter
