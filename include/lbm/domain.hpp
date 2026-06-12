#pragma once

#include <cstddef>
#include <stdexcept>

namespace lbm {

struct Extent3D {
  std::size_t nx{};
  std::size_t ny{};
  std::size_t nz{};

  [[nodiscard]] constexpr std::size_t cell_count() const noexcept {
    return nx * ny * nz;
  }

  [[nodiscard]] constexpr bool empty() const noexcept {
    return nx == 0 || ny == 0 || nz == 0;
  }
};

inline void validate_extent(const Extent3D& extent) {
  if (extent.empty()) {
    throw std::invalid_argument("LBM domain extents must be non-zero");
  }
}

[[nodiscard]] constexpr std::size_t cell_index(
    const Extent3D& extent,
    const std::size_t x,
    const std::size_t y,
    const std::size_t z) noexcept {
  return (z * extent.ny + y) * extent.nx + x;
}

}  // namespace lbm

