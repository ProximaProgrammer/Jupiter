#include "Dynamics.hpp"
#include "EOS.hpp"
#include "Layer.hpp"
#include "../common/Constants.hpp"
#include "../common/Vec3.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace Jupiter
{

namespace
{

double sqr(double x) { return x * x; }


int wrap_phi(int ip, int nphi)
{
    if (nphi <= 0) return 0;
    while (ip < 0) ip += nphi;
    while (ip >= nphi) ip -= nphi;
    return ip;
}

const Cell& clamp_cell(const Grid& grid, int ir, int it, int ip)
{
    ir = std::clamp(ir, 0, grid.NR() - 1);
    it = std::clamp(it, 0, grid.NTheta() - 1);
    ip = wrap_phi(ip, grid.NPhi());
    return grid.At(ir, it, ip);
}

double pressure_at(const Grid& grid, int ir, int it, int ip)
{
    return clamp_cell(grid, ir, it, ip).primitive.pressure;
}

double temperature_at(const Grid& grid, int ir, int it, int ip)
{
    return clamp_cell(grid, ir, it, ip).primitive.temperature;
}

double velocity_component_at(const Grid& grid, int ir, int it, int ip, int component)
{
    const PrimitiveState& p = clamp_cell(grid, ir, it, ip).primitive;
    if (component == 0) return p.velocity_r;
    if (component == 1) return p.velocity_theta;
    return p.velocity_phi;
}

double radial_gradient_pressure(const Grid& grid, int ir, int it, int ip)
{
    const int irm = std::max(0, ir - 1);
    const int irp = std::min(grid.NR() - 1, ir + 1);
    const double r0 = grid.At(irm, it, ip).geometry.r_center;
    const double r1 = grid.At(irp, it, ip).geometry.r_center;
    return (pressure_at(grid, irp, it, ip) - pressure_at(grid, irm, it, ip)) / std::max(r1 - r0, 1.0);
}

double radial_gradient_log_entropy(const Grid& grid, int ir, int it, int ip, double gamma)
{
    const int irm = std::max(0, ir - 1);
    const int irp = std::min(grid.NR() - 1, ir + 1);
    const Cell& a = grid.At(irm, it, ip);
    const Cell& b = grid.At(irp, it, ip);
    const double sa = std::log(std::max(a.primitive.pressure, 1.0e-300)) - gamma * std::log(std::max(a.primitive.density, 1.0e-300));
    const double sb = std::log(std::max(b.primitive.pressure, 1.0e-300)) - gamma * std::log(std::max(b.primitive.density, 1.0e-300));
    return (sb - sa) / std::max(b.geometry.r_center - a.geometry.r_center, 1.0);
}

double enclosed_mass(double r)
{
    // Smooth interior mass model. It avoids treating the shell start at 0.2 R_J
    // as if the whole planet were a point mass while still approaching M_J at
    // the visible radius. Replace this with a real interior density integral later.
    const double x = std::clamp(r / Constants::R_J, 1.0e-6, 1.0);
    const double core = 0.030;
    const double y = std::pow(x, 2.35);
    const double y1 = 1.0 / (1.0 + core);
    return Constants::M_J * (y / (y + core)) / y1;
}

double gravity_magnitude(double r)
{
    return Constants::G * enclosed_mass(r) / std::max(r*r, 1.0);
}

double equilibrium_temperature(double r, double theta, double phi)
{
    const double x = std::clamp(r / Constants::R_J, 0.0, 1.0);
    const double q = 1.0 - x;
    const double deep = 25000.0;
    const double surface = 165.0;
    double Teq = surface + (deep - surface) * std::pow(q, 1.16);

    // A very small nonaxisymmetric upper-atmosphere forcing prevents the RT/
    // cooling source from exactly preserving a one-dimensional state.
    const double top = std::exp(-std::pow((1.0 - x) / 0.075, 2));
    const double lat = Constants::pi / 2.0 - theta;
    Teq *= 1.0 + top * 0.015 * std::cos(2.0 * lat) * std::cos(phi);
    return std::max(2.7, Teq);
}

double radiative_timescale(const Cell& c, const DynamicsConfig& cfg)
{
    const double x = std::clamp(c.geometry.r_center / Constants::R_J, 0.0, 1.0);
    const double w_top = std::exp(-std::pow((1.0 - x) / 0.055, 2));
    const double inv = w_top / std::max(cfg.top_radiative_timescale_s, 1.0) +
                       (1.0 - w_top) / std::max(cfg.deep_radiative_timescale_s, 1.0);
    return 1.0 / std::max(inv, 1.0e-30);
}

double laplacian_scalar(const Grid& grid, int ir, int it, int ip,
                        double center_value,
                        double (*sampler)(const Grid&, int, int, int))
{
    const Cell& c = grid.At(ir, it, ip);
    const auto& g = c.geometry;
    const double dr2 = sqr(std::max(g.length_r, 1.0));
    const double dt2 = sqr(std::max(g.length_theta, 1.0));
    const double dp2 = sqr(std::max(g.length_phi, 1.0));

    const double vrp = sampler(grid, std::min(grid.NR()-1, ir+1), it, ip);
    const double vrm = sampler(grid, std::max(0, ir-1), it, ip);
    const double vtp = sampler(grid, ir, std::min(grid.NTheta()-1, it+1), ip);
    const double vtm = sampler(grid, ir, std::max(0, it-1), ip);
    const double vpp = sampler(grid, ir, it, wrap_phi(ip+1, grid.NPhi()));
    const double vpm = sampler(grid, ir, it, wrap_phi(ip-1, grid.NPhi()));

    return (vrp - 2.0*center_value + vrm) / dr2 +
           (vtp - 2.0*center_value + vtm) / dt2 +
           (vpp - 2.0*center_value + vpm) / dp2;
}

struct FluxMove
{
    int from_ir = 0, from_it = 0, from_ip = 0;
    int to_ir = 0, to_it = 0, to_ip = 0;
    double eps = 0.0;
};

void add_move(Grid& grid, std::vector<ConservedState>& delta, const FluxMove& m)
{
    if (m.eps <= 0.0) return;
    if (!grid.InBounds(m.from_ir, m.from_it, m.from_ip)) return;
    if (!grid.InBounds(m.to_ir, m.to_it, m.to_ip)) return;
    const std::size_t a = grid.Flatten(m.from_ir, m.from_it, m.from_ip);
    const std::size_t b = grid.Flatten(m.to_ir, m.to_it, m.to_ip);
    const ConservedState moved = grid.Cells()[a].conserved * m.eps;
    delta[a] -= moved;
    delta[b] += moved;
}

void advect_face(Grid& grid, std::vector<ConservedState>& delta,
                 int a_ir, int a_it, int a_ip,
                 int b_ir, int b_it, int b_ip,
                 double face_velocity, double face_area, double dt, double eps_max)
{
    if (face_area <= 0.0) return;
    if (face_velocity > 0.0)
    {
        const Cell& donor = grid.At(a_ir, a_it, a_ip);
        const double eps = std::clamp(face_velocity * face_area * dt / std::max(donor.geometry.volume, 1.0), 0.0, eps_max);
        add_move(grid, delta, {a_ir, a_it, a_ip, b_ir, b_it, b_ip, eps});
    }
    else if (face_velocity < 0.0)
    {
        const Cell& donor = grid.At(b_ir, b_it, b_ip);
        const double eps = std::clamp((-face_velocity) * face_area * dt / std::max(donor.geometry.volume, 1.0), 0.0, eps_max);
        add_move(grid, delta, {b_ir, b_it, b_ip, a_ir, a_it, a_ip, eps});
    }
}

void renormalize_species(ConservedState& u)
{
    if (u.mass <= 0.0) return;
    double sum = 0.0;
    for (double m : u.species_mass) sum += std::max(0.0, m);
    if (sum <= 0.0)
    {
        u.species_mass[SIndex(Species::H2)] = 0.75 * u.mass;
        u.species_mass[SIndex(Species::He)] = 0.25 * u.mass;
        return;
    }
    const double scale = u.mass / sum;
    for (double& m : u.species_mass) m = std::max(0.0, m) * scale;
}

} // namespace

void Dynamics::UpdatePrimitives(Grid& grid)
{
    for (Cell& c : grid.Cells()) c.primitive = EOS::ConservedToPrimitive(c.conserved, c.geometry).primitive;
}

void Dynamics::ApplySourceTerms(Grid& grid, double dt, const DynamicsConfig& cfg)
{
    UpdatePrimitives(grid);
    std::vector<ConservedState> delta(grid.Size());

    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const Cell& c = grid.At(ir, it, ip);
        const auto& g = c.geometry;
        const auto& p = c.primitive;
        if (c.conserved.mass <= 0.0) continue;

        double ar = 0.0, at = 0.0, ap = 0.0;

        if (cfg.include_gravity)
            ar += -gravity_magnitude(g.r_center);

        if (cfg.include_pressure_gradient)
        {
            const int itm = std::max(0, it - 1);
            const int itp = std::min(grid.NTheta() - 1, it + 1);
            const int ipm = wrap_phi(ip - 1, grid.NPhi());
            const int ipp = wrap_phi(ip + 1, grid.NPhi());

            const double dpdr = radial_gradient_pressure(grid, ir, it, ip);
            const double dpdt = (pressure_at(grid, ir, itp, ip) - pressure_at(grid, ir, itm, ip)) /
                std::max(grid.At(ir, itp, ip).geometry.theta_center - grid.At(ir, itm, ip).geometry.theta_center, 1.0e-12);
            const double dpdp = (pressure_at(grid, ir, it, ipp) - pressure_at(grid, ir, it, ipm)) /
                std::max(2.0 * g.dphi, 1.0e-12);

            ar += -dpdr / std::max(p.density, 1.0e-30);
            at += -(1.0 / std::max(g.r_center, 1.0)) * dpdt / std::max(p.density, 1.0e-30);
            ap += -(1.0 / std::max(g.r_center * std::sin(g.theta_center), 1.0)) * dpdp / std::max(p.density, 1.0e-30);
        }

        if (cfg.include_convective_buoyancy)
        {
            const double dlnSdr = radial_gradient_log_entropy(grid, ir, it, ip, p.gamma_eff);
            const double pressure_scale_height = std::abs(p.pressure / std::max(std::abs(radial_gradient_pressure(grid, ir, it, ip)), 1.0e-30));
            const double excess = std::clamp(-dlnSdr * pressure_scale_height, 0.0, 4.0);
            if (excess > 0.0)
            {
                const double pattern = 0.65 + 0.35 * std::sin(3.0*g.theta_center + 2.0*std::cos(5.0*g.phi_center));
                ar += cfg.convective_buoyancy_strength * gravity_magnitude(g.r_center) * excess * pattern;
            }
        }

        if (cfg.rotating_frame_sources)
        {
            // Centrifugal acceleration: Omega^2 r sin(theta) [sin(theta)e_r + cos(theta)e_theta].
            const double st = std::sin(g.theta_center), ct = std::cos(g.theta_center);
            ar += Constants::Omega_J * Constants::Omega_J * g.r_center * st * st;
            at += Constants::Omega_J * Constants::Omega_J * g.r_center * st * ct;

            // Coriolis acceleration in local basis from -2 Omega x u.
            const Vec3 u_cart = SphericalVelocityToCartesian(g.theta_center, g.phi_center,
                                                            p.velocity_r, p.velocity_theta, p.velocity_phi);
            const Vec3 omega{0.0, 0.0, Constants::Omega_J};
            const Vec3 a_cor = -2.0 * Cross(omega, u_cart);
            ar += Dot(a_cor, SphericalBasisR(g.theta_center, g.phi_center));
            at += Dot(a_cor, SphericalBasisTheta(g.theta_center, g.phi_center));
            ap += Dot(a_cor, SphericalBasisPhi(g.phi_center));
        }

        if (cfg.velocity_drag_timescale_s > 0.0)
        {
            const double drag = 1.0 / cfg.velocity_drag_timescale_s;
            ar += -p.velocity_r * drag;
            at += -p.velocity_theta * drag;
            ap += -p.velocity_phi * drag;
        }

        ConservedState d;
        d.momentum_r = c.conserved.mass * ar * dt;
        d.momentum_theta = c.conserved.mass * at * dt;
        d.momentum_phi = c.conserved.mass * ap * dt;
        const double power = c.conserved.mass * (p.velocity_r * ar + p.velocity_theta * at + p.velocity_phi * ap);
        d.total_energy = power * dt;

        // Radiative relaxation: a grey stand-in for the much more expensive full
        // multi-bin RT update. The spectrum writer still performs formal I_out =
        // I_in exp(-tau) + S(1-exp(-tau)).
        if (cfg.include_radiative_relaxation)
        {
            const double tau = radiative_timescale(c, cfg);
            const double Teq = equilibrium_temperature(g.r_center, g.theta_center, g.phi_center);
            const double dT = (Teq - p.temperature) * std::clamp(dt / tau, -0.20, 0.20);
            const double cv_rho = p.density * Constants::Rgas / std::max(p.mean_molar_mass * (p.gamma_eff - 1.0), 1.0e-30);
            d.total_energy += cv_rho * dT * g.volume;
        }

        delta[grid.Flatten(ir, it, ip)] += d;
    }

    for (std::size_t i = 0; i < grid.Size(); ++i)
    {
        grid.Cells()[i].conserved += delta[i];
        ClipConserved(grid.Cells()[i].conserved);
        renormalize_species(grid.Cells()[i].conserved);
    }
    UpdatePrimitives(grid);

    if (cfg.include_thermal_diffusion || cfg.include_viscosity || cfg.include_chemistry)
    {
        std::vector<ConservedState> delta2(grid.Size());
        UpdatePrimitives(grid);

        for (int ir = 0; ir < grid.NR(); ++ir)
        for (int it = 0; it < grid.NTheta(); ++it)
        for (int ip = 0; ip < grid.NPhi(); ++ip)
        {
            const Cell& c = grid.At(ir, it, ip);
            const auto& g = c.geometry;
            const auto& p = c.primitive;
            if (c.conserved.mass <= 0.0) continue;
            ConservedState d;

            if (cfg.include_thermal_diffusion && cfg.thermal_diffusivity_m2_s > 0.0)
            {
                const double lapT = laplacian_scalar(grid, ir, it, ip, p.temperature, temperature_at);
                const double dT = std::clamp(cfg.thermal_diffusivity_m2_s * lapT * dt, -0.05*p.temperature, 0.05*p.temperature);
                const double cv_rho = p.density * Constants::Rgas / std::max(p.mean_molar_mass * (p.gamma_eff - 1.0), 1.0e-30);
                d.total_energy += cv_rho * dT * g.volume;
            }

            if (cfg.include_viscosity && cfg.kinematic_viscosity_m2_s > 0.0)
            {
                for (int component = 0; component < 3; ++component)
                {
                    auto sampler = [](const Grid&, int, int, int) { return 0.0; };
                    (void)sampler;
                    const double v0 = (component == 0) ? p.velocity_r : (component == 1 ? p.velocity_theta : p.velocity_phi);
                    const double lapV =
                        (velocity_component_at(grid, std::min(grid.NR()-1, ir+1), it, ip, component) - 2.0*v0 + velocity_component_at(grid, std::max(0, ir-1), it, ip, component)) / sqr(std::max(g.length_r, 1.0)) +
                        (velocity_component_at(grid, ir, std::min(grid.NTheta()-1, it+1), ip, component) - 2.0*v0 + velocity_component_at(grid, ir, std::max(0, it-1), ip, component)) / sqr(std::max(g.length_theta, 1.0)) +
                        (velocity_component_at(grid, ir, it, wrap_phi(ip+1, grid.NPhi()), component) - 2.0*v0 + velocity_component_at(grid, ir, it, wrap_phi(ip-1, grid.NPhi()), component)) / sqr(std::max(g.length_phi, 1.0));
                    const double dv = std::clamp(cfg.kinematic_viscosity_m2_s * lapV * dt, -0.15*std::max(1.0, std::abs(v0)), 0.15*std::max(1.0, std::abs(v0)));
                    if (component == 0) d.momentum_r += c.conserved.mass * dv;
                    else if (component == 1) d.momentum_theta += c.conserved.mass * dv;
                    else d.momentum_phi += c.conserved.mass * dv;
                }
            }

            if (cfg.include_chemistry && cfg.chemistry_timescale_s > 0.0)
            {
                const double x = g.r_center / Constants::R_J;
                const double ion_layer = std::exp(-std::pow((x - 0.982) / 0.035, 2));
                const double hot = 1.0 / (1.0 + std::exp(-(p.temperature - 700.0) / 180.0));
                const double target_X = 2.0e-12 * ion_layer * hot;
                const double current = p.mass_fraction[SIndex(Species::H3plus)];
                const double frac = std::clamp(dt / cfg.chemistry_timescale_s, 0.0, 0.25);
                const double dX = (target_X - current) * frac;
                d.species_mass[SIndex(Species::H3plus)] += c.conserved.mass * dX;
                d.species_mass[SIndex(Species::H2)] -= c.conserved.mass * dX;
            }

            delta2[grid.Flatten(ir, it, ip)] += d;
        }

        for (std::size_t i = 0; i < grid.Size(); ++i)
        {
            grid.Cells()[i].conserved += delta2[i];
            ClipConserved(grid.Cells()[i].conserved);
            renormalize_species(grid.Cells()[i].conserved);
        }
        UpdatePrimitives(grid);
    }

    if (cfg.enforce_boundaries) ApplyBoundaryConditions(grid, cfg);
}

void Dynamics::AdvectFiniteVolume(Grid& grid, double dt, const DynamicsConfig& cfg)
{
    UpdatePrimitives(grid);
    std::vector<ConservedState> delta(grid.Size());

    // Radial faces: positive v_r moves outward, ir -> ir+1.
    for (int ir = 0; ir < grid.NR() - 1; ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const Cell& a = grid.At(ir, it, ip);
        const Cell& b = grid.At(ir + 1, it, ip);
        const double vface = 0.5 * (a.primitive.velocity_r + b.primitive.velocity_r);
        advect_face(grid, delta, ir, it, ip, ir + 1, it, ip,
                    vface, 0.5*(a.geometry.area_r_plus + b.geometry.area_r_minus), dt, cfg.epsilon_max);
    }

    // Theta faces: positive v_theta moves toward increasing theta.
    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta() - 1; ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const Cell& a = grid.At(ir, it, ip);
        const Cell& b = grid.At(ir, it + 1, ip);
        const double vface = 0.5 * (a.primitive.velocity_theta + b.primitive.velocity_theta);
        advect_face(grid, delta, ir, it, ip, ir, it + 1, ip,
                    vface, 0.5*(a.geometry.area_theta_plus + b.geometry.area_theta_minus), dt, cfg.epsilon_max);
    }

    // Phi faces: periodic. Only visit ip -> ip+1, including the wrap face.
    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const int jp = wrap_phi(ip + 1, grid.NPhi());
        const Cell& a = grid.At(ir, it, ip);
        const Cell& b = grid.At(ir, it, jp);
        const double vface = 0.5 * (a.primitive.velocity_phi + b.primitive.velocity_phi);
        advect_face(grid, delta, ir, it, ip, ir, it, jp,
                    vface, 0.5*(a.geometry.area_phi_plus + b.geometry.area_phi_minus), dt, cfg.epsilon_max);
    }

    for (std::size_t i = 0; i < grid.Size(); ++i)
    {
        grid.Cells()[i].conserved += delta[i];
        ClipConserved(grid.Cells()[i].conserved);
        renormalize_species(grid.Cells()[i].conserved);
    }
    UpdatePrimitives(grid);
    if (cfg.enforce_boundaries) ApplyBoundaryConditions(grid, cfg);
}

static double transfer_fraction(double v, double dt, double length, double eps_max)
{
    if (v <= 0.0) return 0.0;
    return std::clamp(v * dt / std::max(length, 1.0), 0.0, eps_max);
}

static void transfer(Grid& grid, std::vector<ConservedState>& delta,
                     int ir0, int it0, int ip0, int ir1, int it1, int ip1, double eps)
{
    if (eps <= 0.0 || !grid.InBounds(ir1, it1, ip1)) return;
    const std::size_t a = grid.Flatten(ir0, it0, ip0);
    const std::size_t b = grid.Flatten(ir1, it1, ip1);
    const ConservedState moved = grid.Cells()[a].conserved * eps;
    delta[a] -= moved;
    delta[b] += moved;
}

void Dynamics::ShiftEpsilon(Grid& grid, double dt, const DynamicsConfig& cfg)
{
    UpdatePrimitives(grid);
    std::vector<ConservedState> delta(grid.Size());

    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        const Cell& c = grid.At(ir, it, ip);
        const auto& p = c.primitive;
        const auto& g = c.geometry;
        if (c.conserved.mass <= 0.0) continue;

        const double erp = transfer_fraction( p.velocity_r, dt, g.length_r, cfg.epsilon_max);
        const double erm = transfer_fraction(-p.velocity_r, dt, g.length_r, cfg.epsilon_max);
        transfer(grid, delta, ir, it, ip, ir + 1, it, ip, erp);
        transfer(grid, delta, ir, it, ip, ir - 1, it, ip, erm);

        const double etp = transfer_fraction( p.velocity_theta, dt, g.length_theta, cfg.epsilon_max);
        const double etm = transfer_fraction(-p.velocity_theta, dt, g.length_theta, cfg.epsilon_max);
        transfer(grid, delta, ir, it, ip, ir, it + 1, ip, etp);
        transfer(grid, delta, ir, it, ip, ir, it - 1, ip, etm);

        const int ipp = wrap_phi(ip + 1, grid.NPhi());
        const int ipm = wrap_phi(ip - 1, grid.NPhi());
        const double epp = transfer_fraction( p.velocity_phi, dt, g.length_phi, cfg.epsilon_max);
        const double epm = transfer_fraction(-p.velocity_phi, dt, g.length_phi, cfg.epsilon_max);
        transfer(grid, delta, ir, it, ip, ir, it, ipp, epp);
        transfer(grid, delta, ir, it, ip, ir, it, ipm, epm);
    }

    for (std::size_t i = 0; i < grid.Size(); ++i)
    {
        grid.Cells()[i].conserved += delta[i];
        ClipConserved(grid.Cells()[i].conserved);
        renormalize_species(grid.Cells()[i].conserved);
    }
    UpdatePrimitives(grid);
    if (cfg.enforce_boundaries) ApplyBoundaryConditions(grid, cfg);
}

void Dynamics::ApplyBoundaryConditions(Grid& grid, const DynamicsConfig& cfg)
{
    (void)cfg;
    UpdatePrimitives(grid);

    for (int it = 0; it < grid.NTheta(); ++it)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        Cell& inner = grid.At(0, it, ip);
        Cell& outer = grid.At(grid.NR() - 1, it, ip);
        if (inner.primitive.velocity_r < 0.0) inner.conserved.momentum_r = 0.0;
        if (outer.primitive.velocity_r > 0.0) outer.conserved.momentum_r = 0.0;
    }

    // At the polar caps, damp theta velocity to avoid coordinate singularity noise.
    for (int ir = 0; ir < grid.NR(); ++ir)
    for (int ip = 0; ip < grid.NPhi(); ++ip)
    {
        grid.At(ir, 0, ip).conserved.momentum_theta = 0.0;
        grid.At(ir, grid.NTheta()-1, ip).conserved.momentum_theta = 0.0;
    }

    for (Cell& c : grid.Cells())
    {
        ClipConserved(c.conserved);
        renormalize_species(c.conserved);
    }
    UpdatePrimitives(grid);
}

} // namespace Jupiter
