#pragma once

#include "lbm/domain.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace lbm {

class BoundaryVelocityField {
 public:
  explicit BoundaryVelocityField(const Extent3D& extent)
      : extent_(extent), velocities_(checked_cell_count(extent), {0.0, 0.0, 0.0}) {}

  [[nodiscard]] const Extent3D& extent() const noexcept {
    return extent_;
  }

  [[nodiscard]] std::size_t cell_count() const noexcept {
    return extent_.cell_count();
  }

  [[nodiscard]] std::array<double, 3>& operator()(const std::size_t cell) {
    return velocities_.at(cell);
  }

  [[nodiscard]] std::array<double, 3> operator()(const std::size_t cell) const {
    return velocities_.at(cell);
  }

  [[nodiscard]] std::array<double, 3>& operator()(
      const std::size_t x,
      const std::size_t y,
      const std::size_t z) {
    return velocities_.at(cell_index(extent_, x, y, z));
  }

  [[nodiscard]] std::array<double, 3> operator()(
      const std::size_t x,
      const std::size_t y,
      const std::size_t z) const {
    return velocities_.at(cell_index(extent_, x, y, z));
  }

  [[nodiscard]] std::span<const std::array<double, 3>> velocities() const noexcept {
    return std::span<const std::array<double, 3>>(velocities_.data(), velocities_.size());
  }

 private:
  [[nodiscard]] static std::size_t checked_cell_count(const Extent3D& extent) {
    validate_extent(extent);
    return extent.cell_count();
  }

  Extent3D extent_{};
  std::vector<std::array<double, 3>> velocities_;
};

}  // namespace lbm

