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
 * @file specbleach.hpp
 * @brief Header-only C++ ownership wrapper for the core single-channel API.
 *
 * Framework-agnostic (standard library only): wraps OWNERSHIP only — RAII
 * lifetime via std::unique_ptr — and adds nothing else: no threading
 * policy, no buffers, no framework types. All processing calls remain the
 * plain C functions in specbleach_denoiser.h, safe to call on .get().
 *
 * Usage:
 * @code
 *   #include <specbleach.hpp>
 *
 *   auto denoiser = specbleach::make_denoiser(48000U, 46.0F);
 *   if (!denoiser) { // handle allocation failure }
 *   specbleach_denoiser_process(denoiser.get(), nframes, in, out);
 *   // freed automatically, exception-safe, non-copyable by construction
 * @endcode
 *
 * Multi-channel users: #include <specbleach_stereo.hpp> instead (it pulls
 * this header in and adds make_stereo_group).
 */

#ifndef SPECBLEACH_HPP_INCLUDED
#define SPECBLEACH_HPP_INCLUDED

#include <cstdint>
#include <memory>

#include "specbleach_denoiser.h"

namespace specbleach {

struct DenoiserDeleter {
  void operator()(specbleach_denoiser* handle) const noexcept {
    specbleach_denoiser_free(handle);
  }
};

using DenoiserPtr = std::unique_ptr<specbleach_denoiser, DenoiserDeleter>;

/// Creates a single-channel spectral denoiser; null on failure.
inline DenoiserPtr make_denoiser(const uint32_t sample_rate,
                                 const float frame_size_ms) noexcept {
  return DenoiserPtr(
      specbleach_denoiser_initialize(sample_rate, frame_size_ms));
}

} // namespace specbleach

#endif /* SPECBLEACH_HPP_INCLUDED */
