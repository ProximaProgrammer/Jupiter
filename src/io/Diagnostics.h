#pragma once

#include "core/Grid.h"
#include "physics/Dynamics.h"

#include <cstddef>
#include <string>

namespace jrt {

class Diagnostics {
public:
    explicit Diagnostics(std::string outputDir);
    void writeStep(const Grid& grid, int step, double time_s, const StepSummary& summary);
    void writeRadialProfile(const Grid& grid, int step, std::size_t j, std::size_t k);
    void writeShellCsv(const Grid& grid, int step);
    void writeShellPly(const Grid& grid, int step);
    void writeVelocityObj(const Grid& grid, int step, std::size_t strideTheta = 3, std::size_t stridePhi = 3);
private:
    std::string outDir;
    bool headerWritten = false;
};

} // namespace jrt
