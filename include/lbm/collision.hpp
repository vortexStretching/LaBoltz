#pragma once

#include "lbm/equilibrium.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/population.hpp"

#include <cstddef>
#include <stdexcept>

namespace lbm {

template <typename Lattice>
void collide_bgk(PopulationField<Lattice>& field, const double omega) {
  if (omega <= 0.0 || omega >= 2.0) {
    throw std::invalid_argument("BGK relaxation frequency omega must be in (0, 2)");
  }

  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    const Macroscopic macro = macroscopic_at(field, cell);
    for (std::size_t q = 0; q < Lattice::q; ++q) {
      const double feq = equilibrium<Lattice>(q, macro.rho, macro.velocity);
      field(q, cell) += omega * (feq - field(q, cell));
    }
  }
}

}  // namespace lbm

