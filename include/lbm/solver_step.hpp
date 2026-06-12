#pragma once

#include "lbm/boundary.hpp"
#include "lbm/collision.hpp"
#include "lbm/geometry.hpp"
#include "lbm/population.hpp"
#include "lbm/streaming.hpp"

#include <array>

namespace lbm {

template <typename Lattice>
void collide_and_stream_periodic(
    PopulationField<Lattice>& current,
    PopulationField<Lattice>& scratch,
    const double omega) {
  collide_bgk(current, omega);
  stream_periodic(current, scratch);
  current = scratch;
}

template <typename Lattice>
void collide_and_stream_bounce_back(
    PopulationField<Lattice>& current,
    PopulationField<Lattice>& scratch,
    const GeometryField& geometry,
    const double omega) {
  collide_bgk(current, geometry, omega);
  stream_periodic_bounce_back(current, scratch, geometry);
  current = scratch;
}

template <typename Lattice>
void collide_and_stream_forced_bounce_back(
    PopulationField<Lattice>& current,
    PopulationField<Lattice>& scratch,
    const GeometryField& geometry,
    const double omega,
    const std::array<double, 3>& acceleration) {
  collide_bgk_guo(current, geometry, omega, acceleration);
  stream_periodic_bounce_back(current, scratch, geometry);
  current = scratch;
}

template <typename Lattice>
void collide_and_stream_moving_bounce_back(
    PopulationField<Lattice>& current,
    PopulationField<Lattice>& scratch,
    const GeometryField& geometry,
    const BoundaryVelocityField& wall_velocity,
    const double omega,
    const double wall_density = 1.0) {
  collide_bgk(current, geometry, omega);
  stream_periodic_moving_bounce_back(current, scratch, geometry, wall_velocity, wall_density);
  current = scratch;
}

}  // namespace lbm
