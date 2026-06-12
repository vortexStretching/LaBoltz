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

```powershell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## First Example

```powershell
.\build\Debug\lbm_periodic_smoke.exe
```

The example runs a tiny periodic D3Q19 domain and writes `periodic_smoke.vtk`.

