#pragma once

#include "lbm/domain.hpp"
#include "lbm/geometry.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/population.hpp"

#include <array>
#include <fstream>
#include <stdexcept>
#include <string>

namespace lbm {

template <typename Lattice>
void write_legacy_vtk(
    const std::string& path,
    const PopulationField<Lattice>& field,
    const std::string& title = "LBM output") {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open VTK output file: " + path);
  }

  const Extent3D& extent = field.extent();
  out << "# vtk DataFile Version 3.0\n";
  out << title << "\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_POINTS\n";
  out << "DIMENSIONS " << extent.nx << ' ' << extent.ny << ' ' << extent.nz << "\n";
  out << "ORIGIN 0 0 0\n";
  out << "SPACING 1 1 1\n";
  out << "POINT_DATA " << extent.cell_count() << "\n";

  out << "SCALARS density double 1\n";
  out << "LOOKUP_TABLE default\n";
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    out << macroscopic_at(field, cell).rho << "\n";
  }

  out << "VECTORS velocity double\n";
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    const auto macro = macroscopic_at(field, cell);
    out << macro.velocity[0] << ' ' << macro.velocity[1] << ' ' << macro.velocity[2] << "\n";
  }
}

template <typename Lattice>
void write_legacy_vtk(
    const std::string& path,
    const PopulationField<Lattice>& field,
    const GeometryField& geometry,
    const std::string& title,
    const std::array<double, 3>& acceleration = {0.0, 0.0, 0.0}) {
  validate_matching_extents(field.extent(), geometry.extent());

  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open VTK output file: " + path);
  }

  const Extent3D& extent = field.extent();
  out << "# vtk DataFile Version 3.0\n";
  out << title << "\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_POINTS\n";
  out << "DIMENSIONS " << extent.nx << ' ' << extent.ny << ' ' << extent.nz << "\n";
  out << "ORIGIN 0 0 0\n";
  out << "SPACING 1 1 1\n";
  out << "POINT_DATA " << extent.cell_count() << "\n";

  out << "SCALARS density double 1\n";
  out << "LOOKUP_TABLE default\n";
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    if (geometry.is_solid(cell)) {
      out << "0\n";
      continue;
    }
    out << macroscopic_at_with_acceleration(field, cell, acceleration).rho << "\n";
  }

  out << "VECTORS velocity double\n";
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    if (geometry.is_solid(cell)) {
      out << "0 0 0\n";
      continue;
    }
    const auto macro = macroscopic_at_with_acceleration(field, cell, acceleration);
    out << macro.velocity[0] << ' ' << macro.velocity[1] << ' ' << macro.velocity[2] << "\n";
  }

  out << "SCALARS cell_type int 1\n";
  out << "LOOKUP_TABLE default\n";
  for (std::size_t cell = 0; cell < field.cell_count(); ++cell) {
    out << (geometry.is_solid(cell) ? 1 : 0) << "\n";
  }
}

}  // namespace lbm
