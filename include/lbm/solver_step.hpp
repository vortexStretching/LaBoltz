#pragma once

#include "lbm/collision.hpp"
#include "lbm/geometry.hpp"
#include "lbm/population.hpp"
#include "lbm/streaming.hpp"

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
  collide_bgk(current, omega);
  stream_periodic_bounce_back(current, scratch, geometry);
  current = scratch;
}

}  // namespace lbm
