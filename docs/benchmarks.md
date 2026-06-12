# Benchmark Roadmap

The verification suite should grow from simple, analytically constrained flows toward increasingly demanding validation cases.

## Tier 1: Numerical Foundation

| Case | Purpose | Primary Checks | Status |
| --- | --- | --- | --- |
| Periodic uniform flow | Streaming, mass conservation | Total mass, uniform velocity preservation | Smoke example |
| Planar Poiseuille flow | No-slip walls, forcing, viscosity | Parabolic profile, convergence history, mass conservation | First benchmark |
| Couette flow | Moving wall boundary condition | Linear profile, wall velocity enforcement | Planned |
| Taylor-Green vortex | Transient viscous decay | Kinetic-energy decay, grid convergence | Planned |

## Tier 2: Canonical CFD Benchmarks

| Case | Purpose | Primary Checks | Status |
| --- | --- | --- | --- |
| Lid-driven cavity | Recirculation and corner vortices | Centerline velocity profiles vs Ghia et al. | Planned |
| Backward-facing step | Separation and reattachment | Reattachment length vs benchmark data | Planned |
| Von Karman vortex street | Unsteady bluff-body shedding | Strouhal number, drag/lift history | Planned |
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
- A reference comparison against analytical, experimental, DNS, or literature data
- A regression test with a short runtime

