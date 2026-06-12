# Architecture Notes

The framework is organized around a small set of separable concepts:

- **Lattice definitions** describe discrete velocities, weights, and opposite directions.
- **Population fields** own distribution functions using `f[q][cell]` Structure-of-Arrays storage.
- **Macroscopic operators** compute density and velocity from populations.
- **Collision operators** relax populations toward equilibrium.
- **Streaming operators** move populations between cells.
- **I/O helpers** write simulation fields for visualization and regression checks.

The current implementation focuses on a single-process CPU kernel. The data layout and interfaces are chosen so that OpenMP, MPI, and GPU backends can be added without changing the public numerical model.

