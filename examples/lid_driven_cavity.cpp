#include "lbm/boundary.hpp"
#include "lbm/geometry.hpp"
#include "lbm/initialization.hpp"
#include "lbm/lattice.hpp"
#include "lbm/macroscopic.hpp"
#include "lbm/output.hpp"
#include "lbm/population.hpp"
#include "lbm/solver_step.hpp"
#include "lbm/vtk.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Config {
  std::size_t nx = 66;
  std::size_t ny = 66;
  std::size_t nz = 4;
  int steps = 8000;
  int report_interval = 1000;
  int vtk_interval = 4000;
  double tau = 0.8;
  double lid_velocity = 0.05;
};

struct Stats {
  double max_speed{};
  double average_speed{};
  double center_ux{};
  double center_uy{};
  double kinetic_energy{};
  double fluid_mass{};
};

void print_usage(const char* executable) {
  std::cout << "Usage: " << executable
            << " [--steps N] [--report N] [--vtk-interval N] [--nx N] [--ny N] [--nz N]"
               " [--tau VALUE] [--lid-velocity VALUE]\n";
}

Config parse_arguments(const int argc, char** argv) {
  Config config;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto require_value = [&](const std::string& name) {
      if (i + 1 >= argc) {
        throw std::invalid_argument("Missing value for " + name);
      }
      return std::string(argv[++i]);
    };

    if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    } else if (arg == "--steps") {
      config.steps = std::stoi(require_value(arg));
    } else if (arg == "--report") {
      config.report_interval = std::stoi(require_value(arg));
    } else if (arg == "--vtk-interval") {
      config.vtk_interval = std::stoi(require_value(arg));
    } else if (arg == "--nx") {
      config.nx = static_cast<std::size_t>(std::stoull(require_value(arg)));
    } else if (arg == "--ny") {
      config.ny = static_cast<std::size_t>(std::stoull(require_value(arg)));
    } else if (arg == "--nz") {
      config.nz = static_cast<std::size_t>(std::stoull(require_value(arg)));
    } else if (arg == "--tau") {
      config.tau = std::stod(require_value(arg));
    } else if (arg == "--lid-velocity") {
      config.lid_velocity = std::stod(require_value(arg));
    } else {
      throw std::invalid_argument("Unknown argument: " + arg);
    }
  }

  if (config.nx < 8 || config.ny < 8 || config.nz == 0) {
    throw std::invalid_argument("Lid-driven cavity requires nx >= 8, ny >= 8, and nz > 0");
  }
  if (config.nx != config.ny) {
    throw std::invalid_argument("The first lid-driven cavity benchmark expects a square nx == ny domain");
  }
  if (config.steps < 0 || config.report_interval <= 0 || config.vtk_interval <= 0) {
    throw std::invalid_argument("Step and interval values must be positive");
  }
  if (config.tau <= 0.5) {
    throw std::invalid_argument("Tau must be greater than 0.5 for positive viscosity");
  }

  return config;
}

template <typename Lattice>
Stats compute_stats(
    const lbm::PopulationField<Lattice>& field,
    const lbm::GeometryField& geometry) {
  Stats stats{};
  double sum_speed = 0.0;
  std::size_t fluid_count = 0;
  const lbm::Extent3D& extent = field.extent();

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t cell = lbm::cell_index(extent, x, y, z);
        if (geometry.is_solid(cell)) {
          continue;
        }

        const auto macro = lbm::macroscopic_at(field, cell);
        const double speed = std::sqrt(
            macro.velocity[0] * macro.velocity[0] +
            macro.velocity[1] * macro.velocity[1] +
            macro.velocity[2] * macro.velocity[2]);

        stats.max_speed = std::max(stats.max_speed, speed);
        sum_speed += speed;
        stats.fluid_mass += macro.rho;
        stats.kinetic_energy += 0.5 * macro.rho * speed * speed;
        ++fluid_count;
      }
    }
  }

  const std::size_t center_x = extent.nx / 2;
  const std::size_t center_y = extent.ny / 2;
  for (std::size_t z = 0; z < extent.nz; ++z) {
    const auto macro = lbm::macroscopic_at(field, lbm::cell_index(extent, center_x, center_y, z));
    stats.center_ux += macro.velocity[0];
    stats.center_uy += macro.velocity[1];
  }
  stats.center_ux /= static_cast<double>(extent.nz);
  stats.center_uy /= static_cast<double>(extent.nz);

  if (fluid_count > 0) {
    stats.average_speed = sum_speed / static_cast<double>(fluid_count);
    stats.kinetic_energy /= static_cast<double>(fluid_count);
  }

  return stats;
}

template <typename Lattice>
void write_profile_csv(
    const std::filesystem::path& path,
    const lbm::PopulationField<Lattice>& field) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open profile CSV: " + path.string());
  }

  const lbm::Extent3D& extent = field.extent();
  const std::size_t center_x = extent.nx / 2;
  const std::size_t center_y = extent.ny / 2;
  const std::size_t sample_count = extent.nx - 2;

  out << "coordinate,vertical_centerline_ux,horizontal_centerline_uy\n";
  for (std::size_t index = 0; index < sample_count; ++index) {
    const std::size_t fluid_coordinate = index + 1;
    double vertical_ux = 0.0;
    double horizontal_uy = 0.0;

    for (std::size_t z = 0; z < extent.nz; ++z) {
      const auto vertical_macro =
          lbm::macroscopic_at(field, lbm::cell_index(extent, center_x, fluid_coordinate, z));
      const auto horizontal_macro =
          lbm::macroscopic_at(field, lbm::cell_index(extent, fluid_coordinate, center_y, z));
      vertical_ux += vertical_macro.velocity[0];
      horizontal_uy += horizontal_macro.velocity[1];
    }

    vertical_ux /= static_cast<double>(extent.nz);
    horizontal_uy /= static_cast<double>(extent.nz);

    const double normalized_coordinate =
        (static_cast<double>(fluid_coordinate) - 0.5) / static_cast<double>(extent.nx - 2);
    out << normalized_coordinate << ',' << vertical_ux << ',' << horizontal_uy << '\n';
  }
}

std::string step_tag(const int step) {
  std::ostringstream stream;
  stream << std::setw(6) << std::setfill('0') << step;
  return stream.str();
}

}  // namespace

int main(int argc, char** argv) {
  using Lattice = lbm::D3Q19;

  try {
    const Config config = parse_arguments(argc, argv);
    const lbm::Extent3D extent{config.nx, config.ny, config.nz};
    const double omega = 1.0 / config.tau;
    const double viscosity = Lattice::cs2 * (config.tau - 0.5);
    const double characteristic_length = static_cast<double>(extent.ny - 2);
    const double reynolds = config.lid_velocity * characteristic_length / viscosity;
    const std::filesystem::path output_directory =
        lbm::create_simulation_output_directory("outputs", "lid_driven_cavity");

    lbm::PopulationField<Lattice> current(extent);
    lbm::PopulationField<Lattice> scratch(extent);
    lbm::GeometryField geometry(extent);
    lbm::BoundaryVelocityField wall_velocity(extent);

    lbm::initialize_equilibrium(current, 1.0, {0.0, 0.0, 0.0});

    for (std::size_t z = 0; z < extent.nz; ++z) {
      for (std::size_t y = 0; y < extent.ny; ++y) {
        geometry(0, y, z) = lbm::CellType::Solid;
        geometry(extent.nx - 1, y, z) = lbm::CellType::Solid;
      }
      for (std::size_t x = 0; x < extent.nx; ++x) {
        geometry(x, 0, z) = lbm::CellType::Solid;
        geometry(x, extent.ny - 1, z) = lbm::CellType::Solid;
        wall_velocity(x, extent.ny - 1, z) = {config.lid_velocity, 0.0, 0.0};
      }
    }

    for (std::size_t cell = 0; cell < current.cell_count(); ++cell) {
      if (!geometry.is_solid(cell)) {
        continue;
      }
      for (std::size_t q = 0; q < Lattice::q; ++q) {
        current(q, cell) = 0.0;
        scratch(q, cell) = 0.0;
      }
    }

    const auto history_path = output_directory / "lid_driven_cavity_history.csv";
    std::ofstream history(history_path);
    if (!history) {
      throw std::runtime_error("Could not open convergence history CSV");
    }
    history << "step,max_speed,average_speed,center_ux,center_uy,kinetic_energy,fluid_mass,relative_max_change\n";

    std::cout << "Lid-driven cavity benchmark\n";
    std::cout << "Domain: " << extent.nx << " x " << extent.ny << " x " << extent.nz << '\n';
    std::cout << "tau=" << config.tau << " omega=" << omega << " viscosity=" << viscosity
              << " lid_velocity=" << config.lid_velocity << " Re=" << reynolds << '\n';
    std::cout << "Output directory: " << output_directory.string() << "\n\n";
    std::cout << std::setw(8) << "step" << std::setw(16) << "max_speed" << std::setw(16)
              << "avg_speed" << std::setw(16) << "center_ux" << std::setw(16)
              << "center_uy" << std::setw(16) << "rel_change" << '\n';

    double previous_max_speed = 0.0;
    auto log_stats = [&](const int step) {
      const Stats stats = compute_stats(current, geometry);
      const double relative_change =
          (std::abs(stats.max_speed) > 0.0)
              ? std::abs(stats.max_speed - previous_max_speed) / std::abs(stats.max_speed)
              : 0.0;
      previous_max_speed = stats.max_speed;

      history << step << ',' << stats.max_speed << ',' << stats.average_speed << ','
              << stats.center_ux << ',' << stats.center_uy << ',' << stats.kinetic_energy
              << ',' << stats.fluid_mass << ',' << relative_change << '\n';

      std::cout << std::setw(8) << step << std::setw(16) << std::scientific
                << stats.max_speed << std::setw(16) << stats.average_speed
                << std::setw(16) << stats.center_ux << std::setw(16)
                << stats.center_uy << std::setw(16) << relative_change << '\n';
    };

    log_stats(0);
    lbm::write_legacy_vtk(
        (output_directory / "lid_driven_cavity_step_000000.vtk").string(),
        current,
        geometry,
        "Lid-driven cavity step 0");

    for (int step = 1; step <= config.steps; ++step) {
      lbm::collide_and_stream_moving_bounce_back(
          current, scratch, geometry, wall_velocity, omega);

      if (step % config.report_interval == 0 || step == 1 || step == config.steps) {
        log_stats(step);
      }

      if (step % config.vtk_interval == 0) {
        lbm::write_legacy_vtk(
            (output_directory / ("lid_driven_cavity_step_" + step_tag(step) + ".vtk")).string(),
            current,
            geometry,
            "Lid-driven cavity");
      }
    }

    lbm::write_legacy_vtk(
        (output_directory / "lid_driven_cavity_final.vtk").string(),
        current,
        geometry,
        "Lid-driven cavity final");
    write_profile_csv(output_directory / "lid_driven_cavity_profile.csv", current);

    std::cout << "\nWrote:\n";
    std::cout << "  " << (output_directory / "lid_driven_cavity_history.csv").string() << '\n';
    std::cout << "  " << (output_directory / "lid_driven_cavity_profile.csv").string() << '\n';
    std::cout << "  " << (output_directory / "lid_driven_cavity_final.vtk").string() << '\n';
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    print_usage(argv[0]);
    return 1;
  }

  return 0;
}
