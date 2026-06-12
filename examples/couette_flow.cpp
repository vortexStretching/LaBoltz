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
  std::size_t nx = 64;
  std::size_t ny = 24;
  std::size_t nz = 4;
  int steps = 6000;
  int report_interval = 1000;
  int vtk_interval = 3000;
  double tau = 0.8;
  double wall_velocity = 0.05;
  double dx_m = 1.0;
  double dt_s = 1.0;
  double rho0_kg_m3 = 1.0;
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
               " [--tau VALUE] [--wall-velocity VALUE] [--dx-m VALUE] [--dt-s VALUE]"
               " [--rho0-kg-m3 VALUE]\n";
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
    } else if (arg == "--wall-velocity") {
      config.wall_velocity = std::stod(require_value(arg));
    } else if (arg == "--dx-m") {
      config.dx_m = std::stod(require_value(arg));
    } else if (arg == "--dt-s") {
      config.dt_s = std::stod(require_value(arg));
    } else if (arg == "--rho0-kg-m3") {
      config.rho0_kg_m3 = std::stod(require_value(arg));
    } else {
      throw std::invalid_argument("Unknown argument: " + arg);
    }
  }

  if (config.nx == 0 || config.ny < 4 || config.nz == 0) {
    throw std::invalid_argument("Couette domain requires nx > 0, ny >= 4, and nz > 0");
  }
  if (config.steps < 0 || config.report_interval <= 0 || config.vtk_interval <= 0) {
    throw std::invalid_argument("Step and interval values must be positive");
  }
  if (config.tau <= 0.5) {
    throw std::invalid_argument("Tau must be greater than 0.5 for positive viscosity");
  }
  if (config.dx_m <= 0.0 || config.dt_s <= 0.0 || config.rho0_kg_m3 <= 0.0) {
    throw std::invalid_argument("SI scaling values dx, dt, and rho0 must be positive");
  }

  return config;
}

double analytical_ux(
    const std::size_t y,
    const lbm::Extent3D& extent,
    const double wall_velocity) {
  const double channel_height = static_cast<double>(extent.ny) - 2.0;
  const double y_from_lower_wall = static_cast<double>(y) - 0.5;
  return wall_velocity * y_from_lower_wall / channel_height;
}

template <typename Lattice>
Stats compute_stats(
    const lbm::PopulationField<Lattice>& field,
    const lbm::GeometryField& geometry,
    const double wall_velocity) {
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

        const auto macro = lbm::macroscopic_at(field, cell);
        const double ux = macro.velocity[0];
        const double reference = analytical_ux(y, extent, wall_velocity);
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
    const double wall_velocity) {
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
        const auto macro = lbm::macroscopic_at(field, cell);
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
        << analytical_ux(y, extent, wall_velocity) << '\n';
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
    const lbm::UnitScale units{config.dx_m, config.dt_s, config.rho0_kg_m3};
    const double omega = 1.0 / config.tau;
    const double viscosity = Lattice::cs2 * (config.tau - 0.5);
    const double channel_height = static_cast<double>(extent.ny) - 2.0;
    const double reynolds = config.wall_velocity * channel_height / viscosity;
    const std::filesystem::path output_directory =
        lbm::create_simulation_output_directory(
            "outputs", "couette", 2, lbm::lattice_name<Lattice>());

    lbm::write_metadata_json(
        output_directory / "metadata.json",
        {
            {"case", lbm::json_string("couette")},
            {"physical_dimension", lbm::json_integer(2)},
            {"solver_dimension", lbm::json_integer(3)},
            {"lattice", lbm::json_string(lbm::lattice_name<Lattice>())},
            {"nx", lbm::json_integer(extent.nx)},
            {"ny", lbm::json_integer(extent.ny)},
            {"nz", lbm::json_integer(extent.nz)},
            {"steps", lbm::json_integer(static_cast<std::size_t>(config.steps))},
            {"tau_lattice", lbm::json_number(config.tau)},
            {"omega_lattice", lbm::json_number(omega)},
            {"kinematic_viscosity_lattice", lbm::json_number(viscosity)},
            {"kinematic_viscosity_m2_s", lbm::json_number(viscosity * units.viscosity_scale_m2_s())},
            {"wall_velocity_lattice", lbm::json_number(config.wall_velocity)},
            {"wall_velocity_m_s", lbm::json_number(config.wall_velocity * units.velocity_scale_m_s())},
            {"reynolds_lattice", lbm::json_number(reynolds)},
            {"dx_m", lbm::json_number(units.dx_m)},
            {"dt_s", lbm::json_number(units.dt_s)},
            {"rho0_kg_m3", lbm::json_number(units.rho0_kg_m3)},
            {"velocity_scale_m_s", lbm::json_number(units.velocity_scale_m_s())},
            {"mass_scale_kg", lbm::json_number(units.mass_scale_kg())},
            {"vtk_field_units", lbm::json_string("lattice units")},
        });

    lbm::PopulationField<Lattice> current(extent);
    lbm::PopulationField<Lattice> scratch(extent);
    lbm::GeometryField geometry(extent);
    lbm::BoundaryVelocityField wall_velocity(extent);

    lbm::initialize_equilibrium(current, 1.0, {0.0, 0.0, 0.0});

    for (std::size_t z = 0; z < extent.nz; ++z) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        geometry(x, 0, z) = lbm::CellType::Solid;
        geometry(x, extent.ny - 1, z) = lbm::CellType::Solid;
        wall_velocity(x, extent.ny - 1, z) = {config.wall_velocity, 0.0, 0.0};
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

    const auto history_path = output_directory / "couette_history.csv";
    std::ofstream history(history_path);
    if (!history) {
      throw std::runtime_error("Could not open convergence history CSV");
    }
    history << "step,time_s,max_ux_lattice,max_ux_m_s,average_ux_lattice,"
               "average_ux_m_s,relative_l2_error,fluid_mass_lattice,fluid_mass_kg,"
               "relative_max_change\n";

    std::cout << "Couette flow benchmark\n";
    std::cout << "Domain: " << extent.nx << " x " << extent.ny << " x " << extent.nz << '\n';
    std::cout << "tau=" << config.tau << " omega=" << omega << " viscosity=" << viscosity
              << " wall_velocity=" << config.wall_velocity << " Re=" << reynolds << "\n\n";
    std::cout << "SI scaling: dx=" << units.dx_m << " m, dt=" << units.dt_s
              << " s, velocity scale=" << units.velocity_scale_m_s() << " m/s\n";
    std::cout << "Output directory: " << output_directory.string() << "\n\n";
    std::cout << std::setw(8) << "step" << std::setw(16) << "max_ux" << std::setw(16)
              << "avg_ux" << std::setw(16) << "L2_error" << std::setw(16) << "rel_change"
              << '\n';

    double previous_max_ux = 0.0;
    auto log_stats = [&](const int step) {
      const Stats stats = compute_stats(current, geometry, config.wall_velocity);
      const double relative_change =
          (std::abs(stats.max_ux) > 0.0)
              ? std::abs(stats.max_ux - previous_max_ux) / std::abs(stats.max_ux)
              : 0.0;
      previous_max_ux = stats.max_ux;

      history << step << ',' << static_cast<double>(step) * units.dt_s << ','
              << stats.max_ux << ',' << stats.max_ux * units.velocity_scale_m_s() << ','
              << stats.average_ux << ',' << stats.average_ux * units.velocity_scale_m_s()
              << ',' << stats.relative_l2_error << ',' << stats.fluid_mass << ','
              << stats.fluid_mass * units.mass_scale_kg() << ',' << relative_change << '\n';

      std::cout << std::setw(8) << step << std::setw(16) << std::scientific << stats.max_ux
                << std::setw(16) << stats.average_ux << std::setw(16)
                << stats.relative_l2_error << std::setw(16) << relative_change << '\n';
    };

    log_stats(0);
    lbm::write_legacy_vtk(
        (output_directory / "couette_step_000000.vtk").string(),
        current,
        geometry,
        "Couette flow step 0");

    for (int step = 1; step <= config.steps; ++step) {
      lbm::collide_and_stream_moving_bounce_back(
          current, scratch, geometry, wall_velocity, omega);

      if (step % config.report_interval == 0 || step == 1 || step == config.steps) {
        log_stats(step);
      }

      if (step % config.vtk_interval == 0) {
        lbm::write_legacy_vtk(
            (output_directory / ("couette_step_" + step_tag(step) + ".vtk")).string(),
            current,
            geometry,
            "Couette flow");
      }
    }

    lbm::write_legacy_vtk(
        (output_directory / "couette_final.vtk").string(),
        current,
        geometry,
        "Couette flow final");
    write_profile_csv(
        output_directory / "couette_profile.csv",
        current,
        geometry,
        config.wall_velocity);

    std::cout << "\nWrote:\n";
    std::cout << "  " << (output_directory / "couette_history.csv").string() << '\n';
    std::cout << "  " << (output_directory / "couette_profile.csv").string() << '\n';
    std::cout << "  " << (output_directory / "couette_final.vtk").string() << '\n';
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    print_usage(argv[0]);
    return 1;
  }

  return 0;
}
