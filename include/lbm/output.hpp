#pragma once

#include <chrono>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lbm {

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
    const std::string& case_name) {
  std::filesystem::create_directories(root);

  const std::string timestamp = simulation_timestamp();
  for (int suffix = 0; suffix < 100; ++suffix) {
    std::ostringstream name;
    name << timestamp << '_' << case_name;
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

}  // namespace lbm
