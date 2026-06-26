#pragma once

#include "core/Grid.h"

namespace jrt {

struct StepSummary {
    double dt = 0.0;
    double maxSignalSpeed = 0.0;
    double maxCorotatingSpeed = 0.0;
    double maxInertialSpeed = 0.0;
};

// Advances the rotating-frame primitive variables by one explicit CFL step.
// The PDEs being discretized are documented in docs/physics_pdes_patch04.md.
StepSummary advanceDynamics(Grid& grid, double requestedDt);

} // namespace jrt
