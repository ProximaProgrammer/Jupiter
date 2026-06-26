#pragma once

#include <vector>
#include "../common/Vec3.hpp"
#include "../core/Grid.hpp"
#include "../io/RuntimeTable.hpp"

namespace Jupiter
{

struct RayResult
{
    std::vector<double> wavelength_m;
    std::vector<double> intensity_W_m3_sr;
};

struct RTConfig
{
    bool use_runtime_tables = true;
    bool include_doppler = true;
    bool include_rotation_velocity = true;
    bool use_cuda_rt_if_available = false;
    double opacity_scale = 1.0;
    double h3plus_emission_scale = 1.0;
};

class RadiativeTransfer
{
public:
    static RayResult IntegrateRadialRay(const Grid& grid,
                                        const RuntimeTables& tables,
                                        int itheta,
                                        int iphi,
                                        const RTConfig& cfg = {});

    static double CellAbsorptionPerMeter(const Cell& cell,
                                         const RuntimeTables& tables,
                                         double wavelength_m,
                                         const RTConfig& cfg);

    static double CellSourceFunction(const Cell& cell,
                                     const RuntimeTables& tables,
                                     double wavelength_m,
                                     const RTConfig& cfg);

    static double ComovingWavelength(double lab_wavelength_m,
                                     const Cell& cell,
                                     const Vec3& photon_direction_cart,
                                     const RTConfig& cfg);
};

} // namespace Jupiter
