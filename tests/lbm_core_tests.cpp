#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/geometry.hpp"
#include "lbm/initialization.hpp"
#include "lbm/lattice.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/population.hpp"
#include "lbm/streaming.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace {

constexpr double tolerance = 1.0e-12;

void require(const bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void test_lattice_weights_and_opposites() {
  static_assert(lbm::lattice_has_valid_opposites<lbm::D3Q19>());
  static_assert(lbm::lattice_has_valid_opposites<lbm::D3Q27>());

  require(std::abs(lbm::lattice_weight_sum<lbm::D3Q19>() - 1.0) < tolerance, "D3Q19 weights do not sum to one");
  require(std::abs(lbm::lattice_weight_sum<lbm::D3Q27>() - 1.0) < tolerance, "D3Q27 weights do not sum to one");
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
    require(std::abs(macro.rho - rho) < tolerance, "Equilibrium density moment is incorrect");
    require(std::abs(macro.velocity[0] - velocity[0]) < tolerance, "Equilibrium ux moment is incorrect");
    require(std::abs(macro.velocity[1] - velocity[1]) < tolerance, "Equilibrium uy moment is incorrect");
    require(std::abs(macro.velocity[2] - velocity[2]) < tolerance, "Equilibrium uz moment is incorrect");
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
  require(destination(2, wrapped_target) == source(2, source_cell), "Periodic streaming did not wrap correctly");
}

void test_streaming_bounce_back_reflects_from_solid_cell() {
  using Lattice = lbm::D3Q19;

  const lbm::Extent3D extent{3, 1, 1};
  lbm::PopulationField<Lattice> source(extent);
  lbm::PopulationField<Lattice> destination(extent);
  lbm::GeometryField geometry(extent);

  const std::size_t left = lbm::cell_index(extent, 0, 0, 0);
  const std::size_t center = lbm::cell_index(extent, 1, 0, 0);
  const std::size_t right = lbm::cell_index(extent, 2, 0, 0);

  geometry(right) = lbm::CellType::Solid;
  source(1, center) = 2.5;

  lbm::stream_periodic_bounce_back(source, destination, geometry);

  require(destination(2, center) == 2.5, "Bounce-back did not reflect into the opposite direction");
  require(destination(1, right) == 0.0, "Bounce-back streamed into a solid cell");
  require(destination(1, left) == 0.0, "Bounce-back affected an unrelated cell");
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

  require(std::abs(lbm::total_mass(scratch) - initial_mass) < tolerance, "Periodic BGK step does not conserve mass");
}

void test_geometry_aware_collision_skips_solid_cells() {
  using Lattice = lbm::D3Q19;

  const lbm::Extent3D extent{2, 1, 1};
  lbm::PopulationField<Lattice> field(extent);
  lbm::GeometryField geometry(extent);

  geometry(1) = lbm::CellType::Solid;
  lbm::initialize_equilibrium(field, 1.0, {0.0, 0.0, 0.0});
  for (std::size_t q = 0; q < Lattice::q; ++q) {
    field(q, 1) = 0.0;
  }

  lbm::collide_bgk(field, geometry, 1.0);

  const auto macro = lbm::macroscopic_at(field, 0);
  require(std::abs(macro.rho - 1.0) < tolerance, "Geometry-aware collision changed fluid density");
}

void test_guo_forcing_preserves_mass_and_adds_momentum() {
  using Lattice = lbm::D3Q19;

  const lbm::Extent3D extent{3, 1, 1};
  lbm::PopulationField<Lattice> field(extent);
  lbm::GeometryField geometry(extent);

  lbm::initialize_equilibrium(field, 1.0, {0.0, 0.0, 0.0});
  const double initial_mass = lbm::total_mass(field);

  lbm::collide_bgk_guo(field, geometry, 1.0, {1.0e-6, 0.0, 0.0});

  require(std::abs(lbm::total_mass(field) - initial_mass) < 1.0e-10, "Guo forcing does not conserve mass");
  const auto macro = lbm::macroscopic_at_with_acceleration(field, 1, {1.0e-6, 0.0, 0.0});
  require(macro.velocity[0] > 0.0, "Guo forcing did not add positive x momentum");
}

}  // namespace

int main() {
  test_lattice_weights_and_opposites();
  test_equilibrium_recovers_macroscopic_values<lbm::D3Q19>();
  test_equilibrium_recovers_macroscopic_values<lbm::D3Q27>();
  test_streaming_is_periodic();
  test_streaming_bounce_back_reflects_from_solid_cell();
  test_bgk_preserves_mass_for_periodic_step();
  test_geometry_aware_collision_skips_solid_cells();
  test_guo_forcing_preserves_mass_and_adds_momentum();

  return 0;
}
