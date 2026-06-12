# CFD_LBM

General-purpose high-performance Lattice Boltzmann framework for fluid and energy research.

This repository is starting with a deliberately small numerical kernel:

- C++20 core solver layer
- D3Q19 and D3Q27 lattice definitions
- Structure-of-Arrays population storage
- BGK collision
- Periodic streaming
- ASCII VTK output for early visualization
- Lightweight CTest smoke tests

## Build

Requirements:

- CMake 3.24 or newer
- A C++20 compiler

On this Windows workstation, the easiest command is:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1
```

That script initializes the Visual Studio Build Tools environment, configures the project, builds it, and runs the test suite.

For longer benchmark runs, use a Release build:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build.ps1 -Configuration Release
```

If CMake and a compiler are already on `PATH`, you can also use the direct commands:

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build -C Debug --output-on-failure
```

## First Example

```powershell
.\build\Debug\lbm_periodic_smoke.exe
```

The example runs a tiny periodic D3Q19 domain and writes `periodic_smoke.vtk`.

## Current Numerical Capabilities

- Periodic streaming
- Solid-cell bounce-back for no-slip walls
- Guo body-force term for pressure-gradient-like forcing
- BGK collision
- D3Q19 and D3Q27 lattices

## Benchmarks

The first verification benchmark is planar Poiseuille flow:

```powershell
.\build\Release\lbm_poiseuille.exe --steps 5000 --report 1000 --vtk-interval 5000
```

It writes VTK and CSV diagnostics under `outputs/`.

The growing benchmark ladder is tracked in [docs/benchmarks.md](docs/benchmarks.md).

## Research Viewer

The current lightweight GUI uses only Python's standard library, so it works even before the future PySide6 interface is installed:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_gui.ps1
```

It loads `outputs/poiseuille_history.csv`, `outputs/poiseuille_profile.csv`, and `outputs/poiseuille_final.vtk`, then shows convergence plots, profile comparison, and a 2D field-slice view. The viewer can also export an HTML report for presentation.

To generate the report directly from the terminal:

```powershell
python gui\laboltz_viewer.py outputs --export-report
```
