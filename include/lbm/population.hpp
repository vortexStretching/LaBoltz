#pragma once

#include "lbm/domain.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace lbm {

template <typename Lattice>
class PopulationField {
 public:
  explicit PopulationField(const Extent3D& extent)
      : extent_(extent), data_(Lattice::q * checked_cell_count(extent), 0.0) {}

  [[nodiscard]] const Extent3D& extent() const noexcept {
    return extent_;
  }

  [[nodiscard]] std::size_t cell_count() const noexcept {
    return extent_.cell_count();
  }

  [[nodiscard]] std::size_t size() const noexcept {
    return data_.size();
  }

  [[nodiscard]] double& operator()(const std::size_t q, const std::size_t cell) {
    return data_.at(offset(q, cell));
  }

  [[nodiscard]] const double& operator()(const std::size_t q, const std::size_t cell) const {
    return data_.at(offset(q, cell));
  }

  [[nodiscard]] std::span<double> direction(const std::size_t q) {
    if (q >= Lattice::q) {
      throw std::out_of_range("Lattice direction is out of range");
    }
    return std::span<double>(data_.data() + q * cell_count(), cell_count());
  }

  [[nodiscard]] std::span<const double> direction(const std::size_t q) const {
    if (q >= Lattice::q) {
      throw std::out_of_range("Lattice direction is out of range");
    }
    return std::span<const double>(data_.data() + q * cell_count(), cell_count());
  }

  [[nodiscard]] std::span<double> data() noexcept {
    return std::span<double>(data_.data(), data_.size());
  }

  [[nodiscard]] std::span<const double> data() const noexcept {
    return std::span<const double>(data_.data(), data_.size());
  }

  void fill(const double value) {
    std::fill(data_.begin(), data_.end(), value);
  }

 private:
  [[nodiscard]] static std::size_t checked_cell_count(const Extent3D& extent) {
    validate_extent(extent);
    return extent.cell_count();
  }

  [[nodiscard]] std::size_t offset(const std::size_t q, const std::size_t cell) const {
    if (q >= Lattice::q) {
      throw std::out_of_range("Lattice direction is out of range");
    }
    if (cell >= cell_count()) {
      throw std::out_of_range("Cell index is outside the population field");
    }
    return q * cell_count() + cell;
  }

  Extent3D extent_{};
  std::vector<double> data_;
};

}  // namespace lbm

