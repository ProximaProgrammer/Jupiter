# JupiterRT Patch 04 physics map

Patch 04 is the first patch where the display is driven by PDE-evolved fields rather than by a mostly static initialization. It is still intentionally explicit and conservative in scope: it evolves the outer Jovian envelope, exports the second-to-outermost shell, and keeps the outermost radial layer as a ghost boundary.

## Rotating-frame state

The grid is Eulerian and co-rotates with Jupiter. The evolved velocity `u` is the co-rotating velocity. Whenever a Doppler-relevant or observer velocity is needed, the code uses

\[
\vec v_\mathrm{inertial}=\vec u+\vec\Omega_J\times\vec r .
\]

That fixes the issue from the previous discussion where the visualizer mostly showed the trivial rotation field.

## PDEs represented in `src/physics/Dynamics.cpp`

The target continuum equations are the compressible rotating-frame hydrodynamic equations, with source terms left explicit so later patches can replace closures with better tables:

\[
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\vec u)=0,
\]

\[
\frac{D\vec u}{Dt}=-\frac{1}{\rho}\nabla p-\nabla\Phi
-2\vec\Omega\times\vec u
-\vec\Omega\times(\vec\Omega\times\vec r)
+\nu\nabla^2\vec u-\frac{\vec u}{\tau_\mathrm{drag}},
\]

where \(\Phi=-GM_J/r\). In spherical components the code applies the centrifugal acceleration

\[
\vec a_\mathrm{cent}=\Omega^2 r\sin\theta(\sin\theta\,\hat e_r+\cos\theta\,\hat e_\theta)
\]

and the Coriolis components for \(\vec\Omega=\Omega(\cos\theta\,\hat e_r-\sin\theta\,\hat e_\theta)\).

The temperature equation is the ideal-gas energy equation in primitive form with explicit turbulent diffusion and radiative relaxation:

\[
\frac{DT}{Dt}=-(\gamma-1)T\nabla\cdot\vec u
+\chi\nabla^2T
-\frac{T-T_\mathrm{eq}(r)}{\tau_\mathrm{rad}(r)}.
\]

The EOS used in this patch is

\[
p=\rho R_\mathrm{spec}T,\quad R_\mathrm{spec}=k_B/\bar m.
\]

That is correct for the current outer molecular envelope approximation. The code keeps the EOS isolated in `src/physics/EOS.cpp` so a later patch can switch the deeper layers to SCvH/ab-initio tables without rewriting the dynamics loop.

## Radiative transfer represented in `src/physics/RadiativeTransfer.cpp`

The radial ray solver uses the cell update you specified:

\[
I_\lambda^\mathrm{out}=I_\lambda^\mathrm{in}e^{-\tau_\lambda}+S_\lambda(1-e^{-\tau_\lambda}),
\quad S_\lambda\approx B_\lambda(T),
\quad \tau_\lambda=\kappa_\lambda\rho\Delta s.
\]

For now \(\kappa_\lambda\) is a physically labeled approximation: grey CIA/cloud opacity plus methane and H3+ shaped features in the JIRAM 2-5 micron range. The function boundary is deliberately clean so Patch 05 can replace this with the generated runtime tables in `data/runtime/*.bin`.

## Boundary conditions

- Radial inner ghost: reflecting radial velocity, copied thermodynamic state.
- Radial outer ghost: open/radiative copy with mild top-temperature relaxation.
- Polar ghosts: reflected meridional velocity.
- Longitude: periodic indexing.

The rendered/exported shell is always `outerInteriorI() = ng + nr - 1`. The true outermost radial layer `outerGhostI() = ng + nr` is not plotted because it is the boundary-condition ghost layer.

## Outputs

- `data/output/diagnostics.csv`: global min/max and CFL summaries.
- `data/output/profiles/radial_theta_phi_latest.csv`: graph-ready radial profile near \(\theta=0,\phi=0\).
- `data/output/synthetic_radial_ray.csv`: 336-bin JIRAM-like radial spectrum from the same angular column.
- `data/output/frames/dashboard_latest.ppm`: labeled 2D dashboard.
- `data/output/shells/outer_interior_shell_latest.csv`: 3D shell points with temperature and velocity.
- `data/output/shells/outer_interior_shell_latest.ply`: colorized 3D mesh.
- `data/output/shells/velocity_shell_latest.obj`: sparse velocity arrows.

