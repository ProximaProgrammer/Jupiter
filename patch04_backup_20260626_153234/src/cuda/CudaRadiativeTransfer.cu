#ifdef JUPITER_ENABLE_CUDA
#include "CudaHooks.hpp"

#include <cuda_runtime.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Jupiter
{

static void cuda_check(cudaError_t e, const char* where)
{
    if (e != cudaSuccess)
        throw std::runtime_error(std::string(where) + ": " + cudaGetErrorString(e));
}

__global__ void rt_step_kernel(double* I, const double* alpha, const double* source, double ds, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double tau_raw = alpha[i] * ds;
    const double tau = fmin(fmax(tau_raw, 0.0), 700.0);
    const double atten = exp(-tau);
    I[i] = I[i] * atten + source[i] * (1.0 - atten);
}

bool CudaBuilt() { return true; }

void CudaApplyRTLayer(std::vector<double>& intensity,
                      const std::vector<double>& alpha_per_m,
                      const std::vector<double>& source_function,
                      double ds_m)
{
    if (intensity.size() != alpha_per_m.size() || intensity.size() != source_function.size())
        throw std::runtime_error("CudaApplyRTLayer received mismatched vector sizes");
    const int n = static_cast<int>(intensity.size());
    if (n == 0) return;

    double* dI = nullptr;
    double* da = nullptr;
    double* dS = nullptr;
    const std::size_t bytes = static_cast<std::size_t>(n) * sizeof(double);

    cuda_check(cudaMalloc(&dI, bytes), "cudaMalloc I");
    cuda_check(cudaMalloc(&da, bytes), "cudaMalloc alpha");
    cuda_check(cudaMalloc(&dS, bytes), "cudaMalloc source");
    cuda_check(cudaMemcpy(dI, intensity.data(), bytes, cudaMemcpyHostToDevice), "copy I to GPU");
    cuda_check(cudaMemcpy(da, alpha_per_m.data(), bytes, cudaMemcpyHostToDevice), "copy alpha to GPU");
    cuda_check(cudaMemcpy(dS, source_function.data(), bytes, cudaMemcpyHostToDevice), "copy source to GPU");

    const int block = 256;
    const int grid = (n + block - 1) / block;
    rt_step_kernel<<<grid, block>>>(dI, da, dS, ds_m, n);
    cuda_check(cudaGetLastError(), "launch rt_step_kernel");
    cuda_check(cudaDeviceSynchronize(), "sync rt_step_kernel");
    cuda_check(cudaMemcpy(intensity.data(), dI, bytes, cudaMemcpyDeviceToHost), "copy I to CPU");

    cudaFree(dI);
    cudaFree(da);
    cudaFree(dS);
}

void RunCudaSmokeTest()
{
    int device = -1;
    cuda_check(cudaGetDevice(&device), "cudaGetDevice");
    cudaDeviceProp prop{};
    cuda_check(cudaGetDeviceProperties(&prop, device), "cudaGetDeviceProperties");
    std::cout << "CUDA device " << device << ": " << prop.name
              << ", SM " << prop.major << "." << prop.minor
              << ", global mem " << static_cast<double>(prop.totalGlobalMem)/(1024.0*1024.0*1024.0)
              << " GiB\n";

    constexpr int n = 1 << 20;
    std::vector<double> I(n, 0.0), alpha(n, 1.0e-7), source(n, 42.0);
    CudaApplyRTLayer(I, alpha, source, 1000.0);
    std::cout << "CUDA RT smoke test complete. I[0]=" << I[0] << "\n";
}

} // namespace Jupiter
#endif
