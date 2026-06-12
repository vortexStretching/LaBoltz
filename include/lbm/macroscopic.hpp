#pragma once

#include "lbm/population.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>

namespace lbm {

struct Macroscopic {
  double rho{};
  std::array<double, 3> velocity{};
};

template <typename Lattice>
[[nodiscard]] Macroscopic macroscopic_at(
    const PopulationField<Lattice>& field,
    const std::size_t cell) {
  if (cell >= field.cell_count()) {
    throw std::out_of_range("Cell index is outside the population field");
  }

  Macroscopic result{};

  for (std::size_t q = 0; q < Lattice::q; ++q) {
    const double f = field(q, cell);
    result.rho += f;
    result.velocity[0] += f * static_cast<double>(Lattice::velocities[q][0]);
    result.velocity[1] += f * static_cast<double>(Lattice::velocities[q][1]);
    result.velocity[2] += f * static_cast<double>(Lattice::velocities[q][2]);
  }

  if (result.rho == 0.0) {
    throw std::runtime_error("Cannot compute velocity for zero-density cell");
  }

  result.velocity[0] /= result.rho;
  result.velocity[1] /= result.rho;
  result.velocity[2] /= result.rho;

  return result;
}

template <typename Lattice>
[[nodiscard]] Macroscopic macroscopic_at_with_acceleration(
    const PopulationField<Lattice>& field,
    const std::size_t cell,
    const std::array<double, 3>& acceleration) {
  Macroscopic result = macroscopic_at(field, cell);

  result.velocity[0] += 0.5 * acceleration[0];
  result.velocity[1] += 0.5 * acceleration[1];
  result.velocity[2] += 0.5 * acceleration[2];

  return result;
}

template <typename Lattice>
[[nodiscard]] double total_mass(const PopulationField<Lattice>& field) noexcept {
  double mass = 0.0;
  for (double value : field.data()) {
    mass += value;
  }
  return mass;
}

}  // namespace lbm
