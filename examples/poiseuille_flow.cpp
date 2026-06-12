#include "lbm/geometry.hpp"
#include "lbm/initialization.hpp"
#include "lbm/lattice.hpp"
#include "lbm/macroscopic.hpp"
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
  std::size_t nx = 64;
  std::size_t ny = 24;
  std::size_t nz = 4;
  int steps = 5000;
  int report_interval = 250;
  int vtk_interval = 2500;
  double tau = 0.8;
  double acceleration = 1.0e-6;
};

struct Stats {
  double max_ux{};
  double average_ux{};
  double relative_l2_error{};
  double fluid_mass{};
};

void print_usage(const char* executable) {
  std::cout << "Usage: " << executable
            << " [--steps N] [--report N] [--vtk-interval N] [--nx N] [--ny N] [--nz N]"
               " [--tau VALUE] [--acceleration VALUE]\n";
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
    } else if (arg == "--acceleration") {
      config.acceleration = std::stod(require_value(arg));
    } else {
      throw std::invalid_argument("Unknown argument: " + arg);
    }
  }

  if (config.nx == 0 || config.ny < 4 || config.nz == 0) {
    throw std::invalid_argument("Poiseuille domain requires nx > 0, ny >= 4, and nz > 0");
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
double analytical_ux(
    const std::size_t y,
    const lbm::Extent3D& extent,
    const double viscosity,
    const double acceleration) {
  const double channel_height = static_cast<double>(extent.ny) - 2.0;
  const double y_from_lower_wall = static_cast<double>(y) - 0.5;
  return acceleration * y_from_lower_wall * (channel_height - y_from_lower_wall) /
         (2.0 * viscosity);
}

template <typename Lattice>
Stats compute_stats(
    const lbm::PopulationField<Lattice>& field,
    const lbm::GeometryField& geometry,
    const std::array<double, 3>& acceleration,
    const double viscosity) {
  Stats stats{};
  double sum_ux = 0.0;
  double error_norm = 0.0;
  double reference_norm = 0.0;
  std::size_t fluid_count = 0;
  const lbm::Extent3D& extent = field.extent();

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t cell = lbm::cell_index(extent, x, y, z);
        if (geometry.is_solid(cell)) {
          continue;
        }

        const auto macro = lbm::macroscopic_at_with_acceleration(field, cell, acceleration);
        const double ux = macro.velocity[0];
        const double reference = analytical_ux<Lattice>(y, extent, viscosity, acceleration[0]);
        const double error = ux - reference;

        stats.max_ux = std::max(stats.max_ux, ux);
        sum_ux += ux;
        stats.fluid_mass += macro.rho;
        error_norm += error * error;
        reference_norm += reference * reference;
        ++fluid_count;
      }
    }
  }

  if (fluid_count > 0) {
    stats.average_ux = sum_ux / static_cast<double>(fluid_count);
  }
  if (reference_norm > 0.0) {
    stats.relative_l2_error = std::sqrt(error_norm / reference_norm);
  }

  return stats;
}

template <typename Lattice>
void write_profile_csv(
    const std::filesystem::path& path,
    const lbm::PopulationField<Lattice>& field,
    const lbm::GeometryField& geometry,
    const std::array<double, 3>& acceleration,
    const double viscosity) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open profile CSV: " + path.string());
  }

  const lbm::Extent3D& extent = field.extent();
  std::vector<double> sum_ux(extent.ny, 0.0);
  std::vector<std::size_t> count(extent.ny, 0);

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t cell = lbm::cell_index(extent, x, y, z);
        if (geometry.is_solid(cell)) {
          continue;
        }
        const auto macro = lbm::macroscopic_at_with_acceleration(field, cell, acceleration);
        sum_ux[y] += macro.velocity[0];
        ++count[y];
      }
    }
  }

  out << "y,y_from_lower_wall,ux_mean,ux_analytical\n";
  for (std::size_t y = 0; y < extent.ny; ++y) {
    if (count[y] == 0) {
      continue;
    }
    const double ux_mean = sum_ux[y] / static_cast<double>(count[y]);
    out << y << ',' << (static_cast<double>(y) - 0.5) << ',' << ux_mean << ','
        << analytical_ux<Lattice>(y, extent, viscosity, acceleration[0]) << '\n';
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
    const double channel_height = static_cast<double>(extent.ny) - 2.0;
    const double analytical_max_ux =
        config.acceleration * channel_height * channel_height / (8.0 * viscosity);
    const std::array<double, 3> acceleration{config.acceleration, 0.0, 0.0};
    const std::filesystem::path output_directory{"outputs"};

    std::filesystem::create_directories(output_directory);

    lbm::PopulationField<Lattice> current(extent);
    lbm::PopulationField<Lattice> scratch(extent);
    lbm::GeometryField geometry(extent);

    lbm::initialize_equilibrium(current, 1.0, {0.0, 0.0, 0.0});

    for (std::size_t z = 0; z < extent.nz; ++z) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        geometry(x, 0, z) = lbm::CellType::Solid;
        geometry(x, extent.ny - 1, z) = lbm::CellType::Solid;
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

    const auto history_path = output_directory / "poiseuille_history.csv";
    std::ofstream history(history_path);
    if (!history) {
      throw std::runtime_error("Could not open convergence history CSV");
    }
    history << "step,max_ux,average_ux,relative_l2_error,fluid_mass,relative_max_change\n";

    std::cout << "Poiseuille flow benchmark\n";
    std::cout << "Domain: " << extent.nx << " x " << extent.ny << " x " << extent.nz << '\n';
    std::cout << "tau=" << config.tau << " omega=" << omega << " viscosity=" << viscosity
              << " acceleration=" << config.acceleration << '\n';
    std::cout << "Analytical max ux: " << analytical_max_ux << "\n\n";
    std::cout << std::setw(8) << "step" << std::setw(16) << "max_ux" << std::setw(16)
              << "avg_ux" << std::setw(16) << "L2_error" << std::setw(16) << "rel_change"
              << '\n';

    double previous_max_ux = 0.0;
    auto log_stats = [&](const int step) {
      const Stats stats = compute_stats(current, geometry, acceleration, viscosity);
      const double relative_change =
          (std::abs(stats.max_ux) > 0.0)
              ? std::abs(stats.max_ux - previous_max_ux) / std::abs(stats.max_ux)
              : 0.0;
      previous_max_ux = stats.max_ux;

      history << step << ',' << stats.max_ux << ',' << stats.average_ux << ','
              << stats.relative_l2_error << ',' << stats.fluid_mass << ',' << relative_change
              << '\n';

      std::cout << std::setw(8) << step << std::setw(16) << std::scientific << stats.max_ux
                << std::setw(16) << stats.average_ux << std::setw(16)
                << stats.relative_l2_error << std::setw(16) << relative_change << '\n';
    };

    log_stats(0);
    lbm::write_legacy_vtk(
        (output_directory / "poiseuille_step_000000.vtk").string(),
        current,
        geometry,
        "Poiseuille flow step 0",
        acceleration);

    for (int step = 1; step <= config.steps; ++step) {
      lbm::collide_and_stream_forced_bounce_back(
          current, scratch, geometry, omega, acceleration);

      if (step % config.report_interval == 0 || step == 1 || step == config.steps) {
        log_stats(step);
      }

      if (step % config.vtk_interval == 0) {
        lbm::write_legacy_vtk(
            (output_directory / ("poiseuille_step_" + step_tag(step) + ".vtk")).string(),
            current,
            geometry,
            "Poiseuille flow",
            acceleration);
      }
    }

    lbm::write_legacy_vtk(
        (output_directory / "poiseuille_final.vtk").string(),
        current,
        geometry,
        "Poiseuille flow final",
        acceleration);
    write_profile_csv(
        output_directory / "poiseuille_profile.csv",
        current,
        geometry,
        acceleration,
        viscosity);

    std::cout << "\nWrote:\n";
    std::cout << "  " << (output_directory / "poiseuille_history.csv").string() << '\n';
    std::cout << "  " << (output_directory / "poiseuille_profile.csv").string() << '\n';
    std::cout << "  " << (output_directory / "poiseuille_final.vtk").string() << '\n';
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    print_usage(argv[0]);
    return 1;
  }

  return 0;
}
