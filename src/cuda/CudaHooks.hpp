#pragma once

#include <vector>

namespace Jupiter
{

bool CudaBuilt();
void RunCudaSmokeTest();

// Applies one RT layer to a spectral vector on the GPU when CUDA is enabled:
// I_lambda <- I_lambda exp(-alpha_lambda ds) + S_lambda (1 - exp(-alpha_lambda ds)).
// The CPU version is intentionally present so physics code can call it unconditionally.
void CudaApplyRTLayer(std::vector<double>& intensity,
                      const std::vector<double>& alpha_per_m,
                      const std::vector<double>& source_function,
                      double ds_m);

} // namespace Jupiter
