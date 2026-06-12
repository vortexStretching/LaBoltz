#pragma once

#include <array>
#include <cstddef>

namespace lbm {

template <typename Lattice>
[[nodiscard]] constexpr double equilibrium(
    const std::size_t q,
    const double rho,
    const std::array<double, 3>& velocity) noexcept {
  const auto& c = Lattice::velocities[q];
  const double cu = static_cast<double>(c[0]) * velocity[0] +
                    static_cast<double>(c[1]) * velocity[1] +
                    static_cast<double>(c[2]) * velocity[2];
  const double u2 = velocity[0] * velocity[0] +
                    velocity[1] * velocity[1] +
                    velocity[2] * velocity[2];

  return Lattice::weights[q] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
}

}  // namespace lbm

