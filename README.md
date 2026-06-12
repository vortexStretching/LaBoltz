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
- BGK collision
- D3Q19 and D3Q27 lattices
