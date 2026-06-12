# Architecture Notes

The framework is organized around a small set of separable concepts:

- **Lattice definitions** describe discrete velocities, weights, and opposite directions.
- **Population fields** own distribution functions using `f[q][cell]` Structure-of-Arrays storage.
- **Geometry fields** identify fluid and solid cells for boundary handling.
- **Boundary velocity fields** store moving-wall velocities independently from solid/fluid geometry.
- **Macroscopic operators** compute density and velocity from populations.
- **Collision operators** relax populations toward equilibrium.
- **Forcing operators** add body forces using a Guo forcing term.
- **Streaming operators** move populations between cells.
- **I/O helpers** write simulation fields for visualization and regression checks.

The current implementation focuses on a single-process CPU kernel. The data layout and interfaces are chosen so that OpenMP, MPI, and GPU backends can be added without changing the public numerical model.
