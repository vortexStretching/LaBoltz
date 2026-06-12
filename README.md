# LaBoltz

LaBoltz is a general-purpose Lattice Boltzmann Method framework for fluid and energy research. The immediate goal is a verified incompressible-flow foundation that can grow into 3D, transient, transport, multiphase, and electrochemical simulations.

This repository is intentionally small at the moment. That is a strength: every benchmark we add should teach us something about the numerical engine.

## Current Status

The current codebase includes:

- C++20 core solver components
- D3Q19 and D3Q27 lattice definitions
- Structure-of-Arrays population storage, `f[q][cell]`
- BGK collision
- Periodic streaming
- Solid-cell bounce-back walls
- Guo body forcing
- Moving-wall support for Couette and lid-driven cavity flow
- ASCII VTK output for visualization
- CSV diagnostics for convergence and validation checks
- SI scaling metadata for benchmark runs
- A lightweight Python/Tk research viewer
- CTest smoke tests

The benchmark executables currently run with D3Q19. The folder names include `d3q19` because the code records which lattice was used. A command-line selector such as `--lattice d3q27` is planned, but it is not wired in yet.

## Repository Map

```text
include/lbm/        Core header-only numerical components
examples/           Standalone benchmark programs
tests/              CTest/unit-style checks
gui/                Lightweight research viewer
scripts/            Convenience scripts for building and running
docs/               Architecture and benchmark notes
outputs/            Generated simulation output, ignored by Git
build/              Generated CMake build tree, ignored by Git
```

## Source Files vs Executables

The `.cpp` files are source code. They are meant for humans and compilers.

For example:

```text
examples/taylor_green_vortex.cpp
```

CMake and the C++ compiler turn that source file into a Windows executable:

```text
build/Release/lbm_taylor_green_vortex.exe
```

You run the `.exe`, not the `.cpp`.

## Requirements

On Windows, the recommended setup is:

- Visual Studio 2022 Build Tools with the C++ workload
- CMake 3.24 or newer
- Python 3 for the viewer

The helper scripts initialize the Visual Studio compiler environment for you, so `cmake` and `cl` do not need to be permanently available on your normal PowerShell `PATH`.

## Build

From the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1
```

This builds the default `Debug` configuration and runs the tests.

For actual simulations, use `Release` because it is optimized and much faster:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1 -Configuration Release
```

The executable files will be created under:

```text
build/Debug/
build/Release/
```

To list all built executables:

```powershell
Get-ChildItem .\build -Recurse -Filter *.exe
```

## Run Benchmarks

Build in `Release` mode first:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1 -Configuration Release
```

Then run one of the benchmark executables.

Planar Poiseuille flow:

```powershell
.\build\Release\lbm_poiseuille.exe --steps 5000 --report 1000 --vtk-interval 5000
```

Planar Couette flow:

```powershell
.\build\Release\lbm_couette.exe --steps 6000 --report 1000 --vtk-interval 6000
```

Lid-driven cavity flow:

```powershell
.\build\Release\lbm_lid_driven_cavity.exe --steps 8000 --report 1000 --vtk-interval 4000
```

Taylor-Green vortex:

```powershell
.\build\Release\lbm_taylor_green_vortex.exe --steps 2000 --report 250 --vtk-interval 500
```

Use `--help` on any executable to see its options:

```powershell
.\build\Release\lbm_taylor_green_vortex.exe --help
```

## Common Run Options

Most benchmarks support:

| Option | Meaning |
| --- | --- |
| `--steps N` | Number of time steps |
| `--report N` | How often diagnostics are printed and written |
| `--vtk-interval N` | How often VTK snapshots are written |
| `--nx N`, `--ny N`, `--nz N` | Domain size in lattice nodes |
| `--tau VALUE` | BGK relaxation time in lattice units |
| `--dx-m VALUE` | Physical lattice spacing in meters |
| `--dt-s VALUE` | Physical time step in seconds |
| `--rho0-kg-m3 VALUE` | Reference density in kg/m^3 |

Case-specific options include:

| Case | Extra option |
| --- | --- |
| Poiseuille | `--acceleration VALUE` |
| Couette | `--wall-velocity VALUE` |
| Lid-driven cavity | `--lid-velocity VALUE` |
| Taylor-Green vortex | `--amplitude VALUE` |

## SI Units

The solver advances in lattice units. SI values are currently added through scaling metadata:

```powershell
.\build\Release\lbm_taylor_green_vortex.exe --dx-m 0.001 --dt-s 0.0001 --rho0-kg-m3 1000
```

This means:

- one lattice spacing is `0.001 m`
- one time step is `0.0001 s`
- the reference density is `1000 kg/m^3`
- the velocity scale is `dx_m / dt_s`
- the viscosity scale is `dx_m^2 / dt_s`
- the mass scale is `rho0_kg_m3 * dx_m^3`

Each run writes these values to `metadata.json` and adds SI columns to the benchmark history CSV where appropriate.

## Output Files

Each benchmark creates a new timestamped output folder:

```text
outputs/YYYYMMDD_HHMMSS_case_name_2d_d3q19/
```

Example:

```text
outputs/20260612_161218_taylor_green_vortex_2d_d3q19/
```

Typical files are:

```text
metadata.json                         Run settings and SI scaling
case_history.csv                      Convergence and diagnostic history
case_profile.csv                      Profile comparison when available
case_step_000500.vtk                  Saved field snapshot
case_final.vtk                        Final field snapshot
case_report.html                      Optional exported viewer report
```

The `2d`/`3d` and `d3q19`/`d3q27` parts are deliberate. They make future comparisons easier to organize.

## Research Viewer

Launch the viewer with:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_gui.ps1
```

The viewer can:

- detect the latest output folder
- switch between benchmark cases
- show convergence history
- show profile comparisons
- inspect VTK field slices
- browse saved transient VTK snapshots
- use a fixed field color range across snapshots
- display SI metadata
- export an HTML report

To generate a report directly from the terminal:

```powershell
python gui\laboltz_viewer.py outputs --case taylor_green_vortex --export-report
```

To check that a case can be loaded:

```powershell
python gui\laboltz_viewer.py outputs --case taylor_green_vortex --check
```

## Changing an Existing Benchmark

The benchmark source files live in `examples/`.

For example, to change Taylor-Green vortex defaults, edit:

```text
examples/taylor_green_vortex.cpp
```

Look for the `Config` struct near the top of the file. That is where default values such as `steps`, `nx`, `ny`, `tau`, and the case-specific velocity scale are set.

After changing C++ code, rebuild:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1 -Configuration Release
```

Then rerun the executable.

## Adding a New Example

Create a new source file under `examples/`, for example:

```text
examples/backward_facing_step.cpp
```

Then register it in `CMakeLists.txt`:

```cmake
add_executable(lbm_backward_facing_step examples/backward_facing_step.cpp)
target_link_libraries(lbm_backward_facing_step PRIVATE lbm::core)
```

The first name passed to `add_executable` becomes the `.exe` name:

```text
build/Release/lbm_backward_facing_step.exe
```

The `lbm_` prefix is only a naming convention. It keeps the project executables easy to recognize.

## Direct CMake Commands

If CMake and a compiler are already on `PATH`, the manual equivalent of the build script is:

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```
