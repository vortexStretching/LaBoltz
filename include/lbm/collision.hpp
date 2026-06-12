#pragma once

#include "lbm/equilibrium.hpp"
#include "lbm/geometry.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/population.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>

namespace lbm {

inline void validate_omega(const double omega) {
  if (omega <= 0.0 || omega >= 2.0) {
    throw std::invalid_argument("BGK relaxation frequency omega must be in (0, 2)");
  }
}

template <typename Lattice>
void collide_bgk(PopulationField<Lattice>& field, const double omega) {
  validate_omega(omega);
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    const Macroscopic macro = macroscopic_at(field, cell);
    for (std::size_t q = 0; q < Lattice::q; ++q) {
      const double feq = equilibrium<Lattice>(q, macro.rho, macro.velocity);
      field(q, cell) += omega * (feq - field(q, cell));
    }
  }
}

template <typename Lattice>
void collide_bgk(
    PopulationField<Lattice>& field,
    const GeometryField& geometry,
    const double omega) {
  validate_omega(omega);
  validate_matching_extents(field.extent(), geometry.extent());

  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    if (geometry.is_solid(cell)) {
      continue;
    }

    const Macroscopic macro = macroscopic_at(field, cell);
    for (std::size_t q = 0; q < Lattice::q; ++q) {
      const double feq = equilibrium<Lattice>(q, macro.rho, macro.velocity);
      field(q, cell) += omega * (feq - field(q, cell));
    }
  }
}

template <typename Lattice>
void collide_bgk_guo(
    PopulationField<Lattice>& field,
    const GeometryField& geometry,
    const double omega,
    const std::array<double, 3>& acceleration) {
  validate_omega(omega);
  validate_matching_extents(field.extent(), geometry.extent());

  constexpr double cs2 = Lattice::cs2;
  constexpr double cs4 = Lattice::cs2 * Lattice::cs2;
  const double force_prefactor = 1.0 - 0.5 * omega;

  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    if (geometry.is_solid(cell)) {
      continue;
    }

    const Macroscopic macro = macroscopic_at_with_acceleration(field, cell, acceleration);
    const std::array<double, 3> force_density{
        macro.rho * acceleration[0],
        macro.rho * acceleration[1],
        macro.rho * acceleration[2],
    };

    for (std::size_t q = 0; q < Lattice::q; ++q) {
      const auto& c = Lattice::velocities[q];
      const double cx = static_cast<double>(c[0]);
      const double cy = static_cast<double>(c[1]);
      const double cz = static_cast<double>(c[2]);
      const double cu = cx * macro.velocity[0] + cy * macro.velocity[1] + cz * macro.velocity[2];

      const double forcing =
          Lattice::weights[q] * force_prefactor *
          ((((cx - macro.velocity[0]) / cs2) + (cu * cx / cs4)) * force_density[0] +
           (((cy - macro.velocity[1]) / cs2) + (cu * cy / cs4)) * force_density[1] +
           (((cz - macro.velocity[2]) / cs2) + (cu * cz / cs4)) * force_density[2]);

      const double feq = equilibrium<Lattice>(q, macro.rho, macro.velocity);
      field(q, cell) += omega * (feq - field(q, cell)) + forcing;
    }
  }
}

}  // namespace lbm
