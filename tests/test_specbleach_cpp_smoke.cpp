/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/**
 * @file test_specbleach_cpp_smoke.cpp
 * @brief Compile-and-run smoke test for the header-only C++ wrappers.
 * Deliberately avoids assertions that vanish under NDEBUG.
 */

#include <specbleach_stereo.hpp>

#include <cstdio>
#include <utility>

namespace {

bool expect(const bool condition, const char* message) {
  if (!condition) {
    std::fprintf(stderr, "TEST FAILED: %s\n", message);
    return false;
  }
  return true;
}

} // namespace

int main() {
  auto denoiser = specbleach::make_denoiser(48000U, 46.0F);
  if (!expect(static_cast<bool>(denoiser), "denoiser created")) {
    return 1;
  }

  auto group = specbleach::make_stereo_group(48000U, 46.0F, 2U);
  if (!expect(static_cast<bool>(group), "stereo group created")) {
    return 1;
  }
  if (!expect(specbleach_stereo_get_channel_count(group.get()) == 2U,
              "group forwards calls through .get()")) {
    return 1;
  }

  specbleach::DenoiserPtr moved = std::move(denoiser);
  if (!expect(!denoiser && static_cast<bool>(moved),
              "move transfers ownership; copy is ill-formed")) {
    return 1;
  }

  std::printf("✅ All C++ wrapper smoke tests passed!\n");
  return 0;
}
