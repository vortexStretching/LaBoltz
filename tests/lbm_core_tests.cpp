#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/initialization.hpp"
#include "lbm/lattice.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/population.hpp"
#include "lbm/streaming.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>

namespace {

constexpr double tolerance = 1.0e-12;

void test_lattice_weights_and_opposites() {
  static_assert(lbm::lattice_has_valid_opposites<lbm::D3Q19>());
  static_assert(lbm::lattice_has_valid_opposites<lbm::D3Q27>());

  assert(std::abs(lbm::lattice_weight_sum<lbm::D3Q19>() - 1.0) < tolerance);
  assert(std::abs(lbm::lattice_weight_sum<lbm::D3Q27>() - 1.0) < tolerance);
}

template <typename Lattice>
void test_equilibrium_recovers_macroscopic_values() {
  const lbm::Extent3D extent{2, 2, 2};
  lbm::PopulationField<Lattice> field(extent);

  constexpr double rho = 1.17;
  constexpr std::array<double, 3> velocity{0.03, -0.02, 0.01};
  lbm::initialize_equilibrium(field, rho, velocity);

  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    const auto macro = lbm::macroscopic_at(field, cell);
    assert(std::abs(macro.rho - rho) < tolerance);
    assert(std::abs(macro.velocity[0] - velocity[0]) < tolerance);
    assert(std::abs(macro.velocity[1] - velocity[1]) < tolerance);
    assert(std::abs(macro.velocity[2] - velocity[2]) < tolerance);
  }
}

void test_streaming_is_periodic() {
  using Lattice = lbm::D3Q19;

  const lbm::Extent3D extent{4, 3, 2};
  lbm::PopulationField<Lattice> source(extent);
  lbm::PopulationField<Lattice> destination(extent);

  for (std::size_t q = 0; q < Lattice::q; ++q) {
    for (std::size_t cell = 0; cell < source.cell_count(); ++cell) {
      source(q, cell) = static_cast<double>(1000 * q + cell);
    }
  }

  lbm::stream_periodic(source, destination);

  const std::size_t source_cell = lbm::cell_index(extent, 0, 0, 0);
  const std::size_t wrapped_target = lbm::cell_index(extent, extent.nx - 1, 0, 0);
  assert(destination(2, wrapped_target) == source(2, source_cell));
}

void test_bgk_preserves_mass_for_periodic_step() {
  using Lattice = lbm::D3Q19;

  const lbm::Extent3D extent{6, 4, 3};
  lbm::PopulationField<Lattice> current(extent);
  lbm::PopulationField<Lattice> scratch(extent);

  lbm::initialize_equilibrium(current, 1.0, {0.02, 0.01, 0.0});
  const double initial_mass = lbm::total_mass(current);

  lbm::collide_bgk(current, 1.0);
  lbm::stream_periodic(current, scratch);

  assert(std::abs(lbm::total_mass(scratch) - initial_mass) < tolerance);
}

}  // namespace

int main() {
  test_lattice_weights_and_opposites();
  test_equilibrium_recovers_macroscopic_values<lbm::D3Q19>();
  test_equilibrium_recovers_macroscopic_values<lbm::D3Q27>();
  test_streaming_is_periodic();
  test_bgk_preserves_mass_for_periodic_step();

  return 0;
}

