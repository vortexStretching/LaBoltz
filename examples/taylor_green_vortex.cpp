#include "lbm/equilibrium.hpp"
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
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct Config {
  std::size_t nx = 64;
  std::size_t ny = 64;
  std::size_t nz = 4;
  int steps = 2000;
  int report_interval = 250;
  int vtk_interval = 500;
  double tau = 0.8;
  double amplitude = 0.04;
  double dx_m = 1.0;
  double dt_s = 1.0;
  double rho0_kg_m3 = 1.0;
};

struct Stats {
  double max_speed{};
  double average_speed{};
  double kinetic_energy{};
  double analytical_kinetic_energy{};
  double relative_energy_error{};
  double relative_velocity_l2_error{};
  double enstrophy{};
  double fluid_mass{};
};

void print_usage(const char* executable) {
  std::cout << "Usage: " << executable
            << " [--steps N] [--report N] [--vtk-interval N] [--nx N] [--ny N] [--nz N]"
               " [--tau VALUE] [--amplitude VALUE] [--dx-m VALUE] [--dt-s VALUE]"
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
    } else if (arg == "--amplitude") {
      config.amplitude = std::stod(require_value(arg));
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

  if (config.nx < 8 || config.ny < 8 || config.nz == 0) {
    throw std::invalid_argument("Taylor-Green vortex requires nx >= 8, ny >= 8, and nz > 0");
  }
  if (config.nx != config.ny) {
    throw std::invalid_argument("The first Taylor-Green benchmark expects nx == ny");
  }
  if (config.steps < 0 || config.report_interval <= 0 || config.vtk_interval <= 0) {
    throw std::invalid_argument("Step and interval values must be positive");
  }
  if (config.tau <= 0.5) {
    throw std::invalid_argument("Tau must be greater than 0.5 for positive viscosity");
  }
  if (config.amplitude <= 0.0 || config.amplitude >= 0.2) {
    throw std::invalid_argument("Amplitude should be in (0, 0.2) for this low-Mach benchmark");
  }
  if (config.dx_m <= 0.0 || config.dt_s <= 0.0 || config.rho0_kg_m3 <= 0.0) {
    throw std::invalid_argument("SI scaling values dx, dt, and rho0 must be positive");
  }

  return config;
}

std::array<double, 3> analytical_velocity(
    const std::size_t x,
    const std::size_t y,
    const lbm::Extent3D& extent,
    const double amplitude,
    const double viscosity,
    const int step) {
  const double wave_number = 2.0 * std::numbers::pi / static_cast<double>(extent.nx);
  const double decay = std::exp(-2.0 * viscosity * wave_number * wave_number * step);
  const double x_position = static_cast<double>(x);
  const double y_position = static_cast<double>(y);

  return {
      amplitude * decay * std::sin(wave_number * x_position) *
          std::cos(wave_number * y_position),
      -amplitude * decay * std::cos(wave_number * x_position) *
          std::sin(wave_number * y_position),
      0.0,
  };
}

template <typename Lattice>
void initialize_taylor_green(
    lbm::PopulationField<Lattice>& field,
    const double amplitude,
    const double viscosity) {
  const lbm::Extent3D& extent = field.extent();
  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t cell = lbm::cell_index(extent, x, y, z);
        const auto velocity = analytical_velocity(x, y, extent, amplitude, viscosity, 0);
        for (std::size_t q = 0; q < Lattice::q; ++q) {
          field(q, cell) = lbm::equilibrium<Lattice>(q, 1.0, velocity);
        }
      }
    }
  }
}

template <typename Lattice>
std::array<double, 2> velocity_xy(
    const lbm::PopulationField<Lattice>& field,
    const std::size_t x,
    const std::size_t y,
    const std::size_t z) {
  const auto macro = lbm::macroscopic_at(field, lbm::cell_index(field.extent(), x, y, z));
  return {macro.velocity[0], macro.velocity[1]};
}

template <typename Lattice>
double vorticity_z(
    const lbm::PopulationField<Lattice>& field,
    const std::size_t x,
    const std::size_t y,
    const std::size_t z) {
  const lbm::Extent3D& extent = field.extent();
  const std::size_t xp = (x + 1) % extent.nx;
  const std::size_t xm = (x + extent.nx - 1) % extent.nx;
  const std::size_t yp = (y + 1) % extent.ny;
  const std::size_t ym = (y + extent.ny - 1) % extent.ny;

  const auto u_yp = velocity_xy(field, x, yp, z)[0];
  const auto u_ym = velocity_xy(field, x, ym, z)[0];
  const auto v_xp = velocity_xy(field, xp, y, z)[1];
  const auto v_xm = velocity_xy(field, xm, y, z)[1];

  return 0.5 * (v_xp - v_xm) - 0.5 * (u_yp - u_ym);
}

template <typename Lattice>
Stats compute_stats(
    const lbm::PopulationField<Lattice>& field,
    const double amplitude,
    const double viscosity,
    const int step) {
  Stats stats{};
  double error_norm = 0.0;
  double reference_norm = 0.0;
  double speed_sum = 0.0;
  const lbm::Extent3D& extent = field.extent();

  for (std::size_t z = 0; z < extent.nz; ++z) {
    for (std::size_t y = 0; y < extent.ny; ++y) {
      for (std::size_t x = 0; x < extent.nx; ++x) {
        const std::size_t cell = lbm::cell_index(extent, x, y, z);
        const auto macro = lbm::macroscopic_at(field, cell);
        const auto reference = analytical_velocity(x, y, extent, amplitude, viscosity, step);

        const double speed2 = macro.velocity[0] * macro.velocity[0] +
                              macro.velocity[1] * macro.velocity[1] +
                              macro.velocity[2] * macro.velocity[2];
        const double speed = std::sqrt(speed2);
        const double error_x = macro.velocity[0] - reference[0];
        const double error_y = macro.velocity[1] - reference[1];
        const double reference2 = reference[0] * reference[0] + reference[1] * reference[1];
        const double omega_z = vorticity_z(field, x, y, z);

        stats.max_speed = std::max(stats.max_speed, speed);
        speed_sum += speed;
        stats.kinetic_energy += 0.5 * macro.rho * speed2;
        stats.fluid_mass += macro.rho;
        stats.enstrophy += 0.5 * omega_z * omega_z;
        error_norm += error_x * error_x + error_y * error_y;
        reference_norm += reference2;
      }
    }
  }

  const double cell_count = static_cast<double>(extent.cell_count());
  stats.average_speed = speed_sum / cell_count;
  stats.kinetic_energy /= cell_count;
  stats.enstrophy /= cell_count;
  stats.analytical_kinetic_energy =
      0.25 * amplitude * amplitude *
      std::exp(-4.0 * viscosity *
               (2.0 * std::numbers::pi / static_cast<double>(extent.nx)) *
               (2.0 * std::numbers::pi / static_cast<double>(extent.nx)) * step);
  if (stats.analytical_kinetic_energy > 0.0) {
    stats.relative_energy_error =
        std::abs(stats.kinetic_energy - stats.analytical_kinetic_energy) /
        stats.analytical_kinetic_energy;
  }
  if (reference_norm > 0.0) {
    stats.relative_velocity_l2_error = std::sqrt(error_norm / reference_norm);
  }

  return stats;
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
    const double reynolds = config.amplitude * static_cast<double>(extent.nx) / viscosity;
    const auto output_directory =
        lbm::create_simulation_output_directory(
            "outputs", "taylor_green_vortex", 2, lbm::lattice_name<Lattice>());

    lbm::write_metadata_json(
        output_directory / "metadata.json",
        {
            {"case", lbm::json_string("taylor_green_vortex")},
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
            {"amplitude_lattice", lbm::json_number(config.amplitude)},
            {"amplitude_m_s", lbm::json_number(config.amplitude * units.velocity_scale_m_s())},
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
    initialize_taylor_green(current, config.amplitude, viscosity);

    const auto history_path = output_directory / "taylor_green_vortex_history.csv";
    std::ofstream history(history_path);
    if (!history) {
      throw std::runtime_error("Could not open convergence history CSV");
    }
    history << "step,time_s,max_speed_lattice,max_speed_m_s,average_speed_lattice,"
               "average_speed_m_s,kinetic_energy_lattice,analytical_kinetic_energy_lattice,"
               "relative_energy_error,relative_velocity_l2_error,enstrophy_lattice,"
               "fluid_mass_lattice,fluid_mass_kg\n";

    std::cout << "Taylor-Green vortex benchmark\n";
    std::cout << "Domain: " << extent.nx << " x " << extent.ny << " x " << extent.nz << '\n';
    std::cout << "tau=" << config.tau << " omega=" << omega << " viscosity=" << viscosity
              << " amplitude=" << config.amplitude << " Re=" << reynolds << '\n';
    std::cout << "SI scaling: dx=" << units.dx_m << " m, dt=" << units.dt_s
              << " s, velocity scale=" << units.velocity_scale_m_s() << " m/s\n";
    std::cout << "Output directory: " << output_directory.string() << "\n\n";
    std::cout << std::setw(8) << "step" << std::setw(16) << "max_speed" << std::setw(16)
              << "kinetic_E" << std::setw(16) << "E_exact" << std::setw(16)
              << "rel_error" << std::setw(16) << "enstrophy" << '\n';

    auto log_stats = [&](const int step) {
      const Stats stats = compute_stats(current, config.amplitude, viscosity, step);
      history << step << ',' << static_cast<double>(step) * units.dt_s << ','
              << stats.max_speed << ',' << stats.max_speed * units.velocity_scale_m_s()
              << ',' << stats.average_speed << ','
              << stats.average_speed * units.velocity_scale_m_s() << ','
              << stats.kinetic_energy << ',' << stats.analytical_kinetic_energy << ','
              << stats.relative_energy_error << ',' << stats.relative_velocity_l2_error << ','
              << stats.enstrophy << ',' << stats.fluid_mass << ','
              << stats.fluid_mass * units.mass_scale_kg() << '\n';

      std::cout << std::setw(8) << step << std::setw(16) << std::scientific
                << stats.max_speed << std::setw(16) << stats.kinetic_energy
                << std::setw(16) << stats.analytical_kinetic_energy << std::setw(16)
                << stats.relative_energy_error << std::setw(16) << stats.enstrophy << '\n';
    };

    log_stats(0);
    lbm::write_legacy_vtk(
        (output_directory / "taylor_green_vortex_step_000000.vtk").string(),
        current,
        "Taylor-Green vortex step 0");

    for (int step = 1; step <= config.steps; ++step) {
      lbm::collide_and_stream_periodic(current, scratch, omega);

      if (step % config.report_interval == 0 || step == 1 || step == config.steps) {
        log_stats(step);
      }

      if (step % config.vtk_interval == 0) {
        lbm::write_legacy_vtk(
            (output_directory / ("taylor_green_vortex_step_" + step_tag(step) + ".vtk")).string(),
            current,
            "Taylor-Green vortex");
      }
    }

    lbm::write_legacy_vtk(
        (output_directory / "taylor_green_vortex_final.vtk").string(),
        current,
        "Taylor-Green vortex final");

    std::cout << "\nWrote:\n";
    std::cout << "  " << (output_directory / "taylor_green_vortex_history.csv").string() << '\n';
    std::cout << "  " << (output_directory / "taylor_green_vortex_final.vtk").string() << '\n';
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    print_usage(argv[0]);
    return 1;
  }

  return 0;
}
