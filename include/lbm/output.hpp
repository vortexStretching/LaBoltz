#pragma once

#include <chrono>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lbm {

struct UnitScale {
  double dx_m = 1.0;
  double dt_s = 1.0;
  double rho0_kg_m3 = 1.0;

  [[nodiscard]] double velocity_scale_m_s() const {
    return dx_m / dt_s;
  }

  [[nodiscard]] double acceleration_scale_m_s2() const {
    return dx_m / (dt_s * dt_s);
  }

  [[nodiscard]] double viscosity_scale_m2_s() const {
    return dx_m * dx_m / dt_s;
  }

  [[nodiscard]] double mass_scale_kg() const {
    return rho0_kg_m3 * dx_m * dx_m * dx_m;
  }
};

inline std::string simulation_timestamp() {
  const auto now = std::chrono::system_clock::now();
  const auto time = std::chrono::system_clock::to_time_t(now);

  std::tm local_time{};
#if defined(_WIN32)
  localtime_s(&local_time, &time);
#else
  localtime_r(&time, &local_time);
#endif

  std::ostringstream stream;
  stream << std::put_time(&local_time, "%Y%m%d_%H%M%S");
  return stream.str();
}

inline std::filesystem::path create_simulation_output_directory(
    const std::filesystem::path& root,
    const std::string& case_name,
    const int physical_dimension,
    const std::string& lattice_name) {
  std::filesystem::create_directories(root);

  const std::string timestamp = simulation_timestamp();
  for (int suffix = 0; suffix < 100; ++suffix) {
    std::ostringstream name;
    name << timestamp << '_' << case_name << '_' << physical_dimension << "d_" << lattice_name;
    if (suffix > 0) {
      name << '_' << std::setw(2) << std::setfill('0') << suffix;
    }

    const auto path = root / name.str();
    if (!std::filesystem::exists(path)) {
      std::filesystem::create_directories(path);
      return path;
    }
  }

  throw std::runtime_error("Could not create a unique simulation output directory");
}

inline std::string json_escape(const std::string& value) {
  std::ostringstream stream;
  for (const char character : value) {
    switch (character) {
      case '"':
        stream << "\\\"";
        break;
      case '\\':
        stream << "\\\\";
        break;
      case '\n':
        stream << "\\n";
        break;
      case '\r':
        stream << "\\r";
        break;
      case '\t':
        stream << "\\t";
        break;
      default:
        stream << character;
        break;
    }
  }
  return stream.str();
}

inline std::string json_string(const std::string& value) {
  return "\"" + json_escape(value) + "\"";
}

inline std::string json_number(const double value) {
  std::ostringstream stream;
  stream << std::setprecision(17) << value;
  return stream.str();
}

inline std::string json_integer(const std::size_t value) {
  return std::to_string(value);
}

inline void write_metadata_json(
    const std::filesystem::path& path,
    const std::vector<std::pair<std::string, std::string>>& entries) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open metadata file: " + path.string());
  }

  out << "{\n";
  for (std::size_t index = 0; index < entries.size(); ++index) {
    out << "  \"" << json_escape(entries[index].first) << "\": " << entries[index].second;
    out << (index + 1 == entries.size() ? "\n" : ",\n");
  }
  out << "}\n";
}

}  // namespace lbm
