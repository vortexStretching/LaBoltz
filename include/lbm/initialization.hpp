#pragma once

#include "lbm/equilibrium.hpp"
#include "lbm/population.hpp"

#include <array>

namespace lbm {

template <typename Lattice>
void initialize_equilibrium(
    PopulationField<Lattice>& field,
    const double rho,
    const std::array<double, 3>& velocity) {
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    for (std::size_t q = 0; q < Lattice::q; ++q) {
      field(q, cell) = equilibrium<Lattice>(q, rho, velocity);
    }
  }
}

}  // namespace lbm

