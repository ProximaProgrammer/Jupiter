#include "physics/TimeStep.h"
#include "physics/EOS.h"

#include <algorithm>
#include <cmath>

namespace jrt {

double estimateStableDt(const Grid& g) {
    double dt = g.cfg.max_dt_s;
    for (std::size_t i = g.iBegin(); i < g.iEnd(); ++i) {
        for (std::size_t j = g.jBegin(); j < g.jEnd(); ++j) {
            const double r = g.radius(i);
            const double s = std::max(0.08, std::sin(g.theta(j)));
            const double dx = std::min({g.dr, r*g.dtheta, r*s*g.dphi});
            for (std::size_t k = g.kBegin(); k < g.kEnd(); ++k) {
                const Cell& c = g.at(i,j,k);
                const double cs = soundSpeedIdealGas(g.cfg.gamma, c.p, c.rho);
                const double u = std::sqrt(c.u.r*c.u.r + c.u.th*c.u.th + c.u.ph*c.u.ph);
                dt = std::min(dt, g.cfg.cfl * dx / std::max(1.0, cs + u));
            }
        }
    }
    return std::clamp(dt, g.cfg.min_dt_s, g.cfg.max_dt_s);
}

} // namespace jrt
