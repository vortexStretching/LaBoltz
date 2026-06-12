#pragma once

#include "lbm/boundary.hpp"
#include "lbm/domain.hpp"
#include "lbm/geometry.hpp"
#include "lbm/population.hpp"

#include <cstddef>
#include <stdexcept>

namespace lbm {

[[nodiscard]] inline std::size_t periodic_offset(
    const std::size_t coordinate,
    const int delta,
    const std::size_t extent) noexcept {
  const auto shifted = static_cast<long long>(coordinate) + static_cast<long long>(delta);
  const auto wrapped = (shifted % static_cast<long long>(extent) + static_cast<long long>(extent)) %
                       static_cast<long long>(extent);
  return static_cast<std::size_t>(wrapped);
}

template <typename Lattice>
void stream_periodic(
    const PopulationField<Lattice>& source,
    PopulationField<Lattice>& destination) {
  validate_matching_extents(source.extent(), destination.extent());

  const Extent3D& extent = source.extent();

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t source_cell = cell_index(extent, x, y, z);

        for (std::size_t q = 0; q < Lattice::q; ++q) {
          const auto& c = Lattice::velocities[q];
          const std::size_t target_x = periodic_offset(x, c[0], extent.nx);
          const std::size_t target_y = periodic_offset(y, c[1], extent.ny);
          const std::size_t target_z = periodic_offset(z, c[2], extent.nz);
          const std::size_t target_cell = cell_index(extent, target_x, target_y, target_z);

          destination(q, target_cell) = source(q, source_cell);
        }
      }
    }
  }
}

template <typename Lattice>
void stream_periodic_bounce_back(
    const PopulationField<Lattice>& source,
    PopulationField<Lattice>& destination,
    const GeometryField& geometry) {
  validate_matching_extents(source.extent(), destination.extent());
  validate_matching_extents(source.extent(), geometry.extent());

  const Extent3D& extent = source.extent();
  destination.fill(0.0);

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t source_cell = cell_index(extent, x, y, z);

        if (geometry.is_solid(source_cell)) {
          continue;
        }

        for (std::size_t q = 0; q < Lattice::q; ++q) {
          const auto& c = Lattice::velocities[q];
          const std::size_t target_x = periodic_offset(x, c[0], extent.nx);
          const std::size_t target_y = periodic_offset(y, c[1], extent.ny);
          const std::size_t target_z = periodic_offset(z, c[2], extent.nz);
          const std::size_t target_cell = cell_index(extent, target_x, target_y, target_z);

          if (geometry.is_fluid(target_cell)) {
            destination(q, target_cell) = source(q, source_cell);
          } else {
            destination(Lattice::opposite[q], source_cell) = source(q, source_cell);
          }
        }
      }
    }
  }
}

template <typename Lattice>
void stream_periodic_moving_bounce_back(
    const PopulationField<Lattice>& source,
    PopulationField<Lattice>& destination,
    const GeometryField& geometry,
    const BoundaryVelocityField& wall_velocity,
    const double wall_density = 1.0) {
  validate_matching_extents(source.extent(), destination.extent());
  validate_matching_extents(source.extent(), geometry.extent());
  validate_matching_extents(source.extent(), wall_velocity.extent());

  const Extent3D& extent = source.extent();
  destination.fill(0.0);

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t source_cell = cell_index(extent, x, y, z);

        if (geometry.is_solid(source_cell)) {
          continue;
        }

        for (std::size_t q = 0; q < Lattice::q; ++q) {
          const auto& c = Lattice::velocities[q];
          const std::size_t target_x = periodic_offset(x, c[0], extent.nx);
          const std::size_t target_y = periodic_offset(y, c[1], extent.ny);
          const std::size_t target_z = periodic_offset(z, c[2], extent.nz);
          const std::size_t target_cell = cell_index(extent, target_x, target_y, target_z);

          if (geometry.is_fluid(target_cell)) {
            destination(q, target_cell) = source(q, source_cell);
          } else {
            const auto velocity = wall_velocity(target_cell);
            const double c_dot_wall_velocity =
                static_cast<double>(c[0]) * velocity[0] +
                static_cast<double>(c[1]) * velocity[1] +
                static_cast<double>(c[2]) * velocity[2];
            const double moving_wall_correction =
                2.0 * Lattice::weights[q] * wall_density * c_dot_wall_velocity /
                Lattice::cs2;
            destination(Lattice::opposite[q], source_cell) =
                source(q, source_cell) - moving_wall_correction;
          }
        }
      }
    }
  }
}

}  // namespace lbm
