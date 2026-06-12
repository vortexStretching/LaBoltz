#pragma once

#include "lbm/domain.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace lbm {

enum class CellType {
  Fluid,
  Solid,
};

class GeometryField {
 public:
  explicit GeometryField(const Extent3D& extent)
      : extent_(extent), cells_(checked_cell_count(extent), CellType::Fluid) {}

  [[nodiscard]] const Extent3D& extent() const noexcept {
    return extent_;
  }

  [[nodiscard]] std::size_t cell_count() const noexcept {
    return extent_.cell_count();
  }

  [[nodiscard]] CellType& operator()(const std::size_t cell) {
    return cells_.at(cell);
  }

  [[nodiscard]] CellType operator()(const std::size_t cell) const {
    return cells_.at(cell);
  }

  [[nodiscard]] CellType& operator()(
      const std::size_t x,
      const std::size_t y,
      const std::size_t z) {
    return cells_.at(cell_index(extent_, x, y, z));
  }

  [[nodiscard]] CellType operator()(
      const std::size_t x,
      const std::size_t y,
      const std::size_t z) const {
    return cells_.at(cell_index(extent_, x, y, z));
  }

  [[nodiscard]] bool is_fluid(const std::size_t cell) const {
    return cells_.at(cell) == CellType::Fluid;
  }

  [[nodiscard]] bool is_solid(const std::size_t cell) const {
    return cells_.at(cell) == CellType::Solid;
  }

  [[nodiscard]] std::span<const CellType> cells() const noexcept {
    return std::span<const CellType>(cells_.data(), cells_.size());
  }

 private:
  [[nodiscard]] static std::size_t checked_cell_count(const Extent3D& extent) {
    validate_extent(extent);
    return extent.cell_count();
  }

  Extent3D extent_{};
  std::vector<CellType> cells_;
};

inline void validate_matching_extents(const Extent3D& lhs, const Extent3D& rhs) {
  if (lhs.nx != rhs.nx || lhs.ny != rhs.ny || lhs.nz != rhs.nz) {
    throw std::invalid_argument("Fields must have identical extents");
  }
}

}  // namespace lbm

