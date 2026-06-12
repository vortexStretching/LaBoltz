#pragma once

#include <array>
#include <cstddef>

namespace lbm {

struct D3Q19 {
  static constexpr std::size_t dimensions = 3;
  static constexpr std::size_t q = 19;
  static constexpr double cs2 = 1.0 / 3.0;

  inline static constexpr std::array<std::array<int, 3>, q> velocities{{
      {{0, 0, 0}},
      {{1, 0, 0}},
      {{-1, 0, 0}},
      {{0, 1, 0}},
      {{0, -1, 0}},
      {{0, 0, 1}},
      {{0, 0, -1}},
      {{1, 1, 0}},
      {{-1, -1, 0}},
      {{1, -1, 0}},
      {{-1, 1, 0}},
      {{1, 0, 1}},
      {{-1, 0, -1}},
      {{1, 0, -1}},
      {{-1, 0, 1}},
      {{0, 1, 1}},
      {{0, -1, -1}},
      {{0, 1, -1}},
      {{0, -1, 1}},
  }};

  inline static constexpr std::array<double, q> weights{
      1.0 / 3.0,
      1.0 / 18.0,
      1.0 / 18.0,
      1.0 / 18.0,
      1.0 / 18.0,
      1.0 / 18.0,
      1.0 / 18.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
      1.0 / 36.0,
  };

  inline static constexpr std::array<std::size_t, q> opposite{
      0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
};

struct D3Q27 {
  static constexpr std::size_t dimensions = 3;
  static constexpr std::size_t q = 27;
  static constexpr double cs2 = 1.0 / 3.0;

  inline static constexpr std::array<std::array<int, 3>, q> velocities{{
      {{0, 0, 0}},
      {{1, 0, 0}},
      {{-1, 0, 0}},
      {{0, 1, 0}},
      {{0, -1, 0}},
      {{0, 0, 1}},
      {{0, 0, -1}},
      {{1, 1, 0}},
      {{-1, -1, 0}},
      {{1, -1, 0}},
      {{-1, 1, 0}},
      {{1, 0, 1}},
      {{-1, 0, -1}},
      {{1, 0, -1}},
      {{-1, 0, 1}},
      {{0, 1, 1}},
      {{0, -1, -1}},
      {{0, 1, -1}},
      {{0, -1, 1}},
      {{1, 1, 1}},
      {{-1, -1, -1}},
      {{1, 1, -1}},
      {{-1, -1, 1}},
      {{1, -1, 1}},
      {{-1, 1, -1}},
      {{-1, 1, 1}},
      {{1, -1, -1}},
  }};

  inline static constexpr std::array<double, q> weights{
      8.0 / 27.0,
      2.0 / 27.0,
      2.0 / 27.0,
      2.0 / 27.0,
      2.0 / 27.0,
      2.0 / 27.0,
      2.0 / 27.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 54.0,
      1.0 / 216.0,
      1.0 / 216.0,
      1.0 / 216.0,
      1.0 / 216.0,
      1.0 / 216.0,
      1.0 / 216.0,
      1.0 / 216.0,
      1.0 / 216.0,
  };

  inline static constexpr std::array<std::size_t, q> opposite{
      0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14,
      13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
};

template <typename Lattice>
[[nodiscard]] constexpr bool lattice_has_valid_opposites() noexcept {
  for (std::size_t q = 0; q < Lattice::q; ++q) {
    const auto opposite = Lattice::opposite[q];
    if (opposite >= Lattice::q || Lattice::opposite[opposite] != q) {
      return false;
    }

    const auto& c = Lattice::velocities[q];
    const auto& co = Lattice::velocities[opposite];
    if (c[0] + co[0] != 0 || c[1] + co[1] != 0 || c[2] + co[2] != 0) {
      return false;
    }
  }

  return true;
}

template <typename Lattice>
[[nodiscard]] constexpr double lattice_weight_sum() noexcept {
  double sum = 0.0;
  for (const double weight : Lattice::weights) {
    sum += weight;
  }
  return sum;
}

template <typename Lattice>
[[nodiscard]] constexpr const char* lattice_name() noexcept {
  return "unknown_lattice";
}

template <>
[[nodiscard]] constexpr const char* lattice_name<D3Q19>() noexcept {
  return "d3q19";
}

template <>
[[nodiscard]] constexpr const char* lattice_name<D3Q27>() noexcept {
  return "d3q27";
}

static_assert(lattice_has_valid_opposites<D3Q19>());
static_assert(lattice_has_valid_opposites<D3Q27>());

}  // namespace lbm
