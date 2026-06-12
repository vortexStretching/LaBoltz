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
- **I/O helpers** create timestamped output folders, write SI scaling metadata, and export simulation fields for visualization and regression checks.

The current implementation focuses on a single-process CPU kernel. The data layout and interfaces are chosen so that OpenMP, MPI, and GPU backends can be added without changing the public numerical model.

Benchmark programs are deliberately standalone. Each one owns its case-specific setup, command-line options, diagnostics, and output schedule while reusing the same core lattice, population, collision, streaming, boundary, and VTK helpers. This keeps early verification cases easy to read and modify.

The Python viewer is independent from the C++ solver. It reads generated CSV, VTK, and JSON metadata files from `outputs/`, which preserves the headless execution path needed for future HPC runs.
