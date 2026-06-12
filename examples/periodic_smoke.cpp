#include "lbm/initialization.hpp"
#include "lbm/lattice.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/population.hpp"
#include "lbm/solver_step.hpp"
#include "lbm/vtk.hpp"

#include <array>
#include <cmath>
#include <iostream>

int main() {
  using Lattice = lbm::D3Q19;

  const lbm::Extent3D extent{16, 8, 8};
  lbm::PopulationField<Lattice> current(extent);
  lbm::PopulationField<Lattice> scratch(extent);

  constexpr double rho = 1.0;
  constexpr std::array<double, 3> velocity{0.04, 0.0, 0.0};
  constexpr double omega = 1.0;

  lbm::initialize_equilibrium(current, rho, velocity);

  const double initial_mass = lbm::total_mass(current);
  for (int step = 0; step < 25; ++step) {
    lbm::collide_and_stream_periodic(current, scratch, omega);
  }
  const double final_mass = lbm::total_mass(current);

  lbm::write_legacy_vtk("periodic_smoke.vtk", current, "Periodic D3Q19 smoke test");

  std::cout << "Initial mass: " << initial_mass << '\n';
  std::cout << "Final mass:   " << final_mass << '\n';
  std::cout << "Mass error:   " << std::abs(final_mass - initial_mass) << '\n';
  std::cout << "Wrote periodic_smoke.vtk\n";

  return 0;
}

