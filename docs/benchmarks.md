# Benchmark Roadmap

The verification suite should grow from simple, analytically constrained flows toward increasingly demanding validation cases.

## Current Executables

| Executable | Source | Purpose |
| --- | --- | --- |
| `lbm_periodic_smoke.exe` | `examples/periodic_smoke.cpp` | Minimal periodic streaming and mass-conservation smoke test |
| `lbm_poiseuille.exe` | `examples/poiseuille_flow.cpp` | Pressure-gradient-like channel flow with no-slip walls |
| `lbm_couette.exe` | `examples/couette_flow.cpp` | Moving-wall channel flow with a linear analytical profile |
| `lbm_lid_driven_cavity.exe` | `examples/lid_driven_cavity.cpp` | Recirculating cavity flow with a moving lid |
| `lbm_taylor_green_vortex.exe` | `examples/taylor_green_vortex.cpp` | Transient viscous vortex decay |

Current benchmark executables use D3Q19. D3Q27 is defined in the core lattice layer and will become a selectable benchmark option once the command-line dispatch is added.

## Tier 1: Numerical Foundation

| Case | Purpose | Primary Checks | Status |
| --- | --- | --- | --- |
| Periodic uniform flow | Streaming, mass conservation | Total mass, uniform velocity preservation | Smoke example |
| Planar Poiseuille flow | No-slip walls, forcing, viscosity | Parabolic profile, convergence history, mass conservation | First benchmark |
| Couette flow | Moving wall boundary condition | Linear profile, wall velocity enforcement | First benchmark |
| Taylor-Green vortex | Transient viscous decay | Kinetic-energy decay, grid convergence | First transient benchmark |

## Tier 2: Canonical CFD Benchmarks

| Case | Purpose | Primary Checks | Status |
| --- | --- | --- | --- |
| Lid-driven cavity | Recirculation and corner vortices | Centerline velocity profiles vs Ghia et al. | First benchmark |
| Backward-facing step | Separation and reattachment | Reattachment length vs benchmark data | Planned |
| Von Karman vortex street | Unsteady bluff-body shedding | Strouhal number, drag/lift history | Planned next transient external-flow case |
| Taylor-Couette flow | Rotating-wall flow and instability onset | Azimuthal velocity profile, torque | Planned |

## Tier 3: Application-Oriented Extensions

| Case | Purpose | Primary Checks | Status |
| --- | --- | --- | --- |
| Porous channel | Complex solid geometry | Permeability, pressure drop | Planned |
| Simplified straw | High-speed internal flow precursor for cavitation studies | Pressure minimum, acceleration regions | Planned |
| Bubble in quiescent liquid | Multiphase foundation | Spurious currents, Laplace law | Future physics |
| Gas-evolving electrode | Electrochemical multiphase target | Bubble growth and departure trends | Future physics |

Each benchmark should eventually provide:

- A standalone executable or input deck
- VTK/HDF5 visualization output
- CSV convergence and diagnostic history
- SI scaling metadata with `dx_m`, `dt_s`, `rho0_kg_m3`, and derived physical units
- A reference comparison against analytical, experimental, DNS, or literature data
- A regression test with a short runtime

Current benchmark output folders include the solver dimension and lattice, for example:

```text
outputs\YYYYMMDD_HHMMSS_taylor_green_vortex_2d_d3q19\
```

The standard output contract is:

- `metadata.json` for run parameters, lattice choice, solver dimension, physical dimension, and SI scaling
- `*_history.csv` for time history and convergence diagnostics
- `*_profile.csv` when an analytical or reference profile is available
- `*_step_XXXXXX.vtk` for transient visualization snapshots
- `*_final.vtk` for the final field
- `*_report.html` when exported from the research viewer

The viewer should keep color ranges fixed across VTK snapshots for a selected field so that transient decay, growth, and approach to steady state remain visually meaningful.
