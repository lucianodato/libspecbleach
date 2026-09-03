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
 * @file specbleach_stereo.hpp
 * @brief Header-only C++ ownership wrapper for the extras stereo API.
 *
 * Includes <specbleach.hpp> (core wrapper) and adds the multi-channel
 * group. Handles are movable, never copyable, so a double-free cannot be
 * expressed.
 *
 * Usage:
 * @code
 *   #include <specbleach_stereo.hpp>
 *
 *   auto group = specbleach::make_stereo_group(48000U, 46.0F, 2U);
 *   specbleach_stereo_process(group.get(), nframes, in, out);
 * @endcode
 */

#ifndef SPECBLEACH_STEREO_HPP_INCLUDED
#define SPECBLEACH_STEREO_HPP_INCLUDED

#include <cstdint>
#include <memory>

#include "specbleach.hpp"
#include "specbleach_stereo.h"

namespace specbleach {

struct StereoGroupDeleter {
  void operator()(specbleach_stereo* handle) const noexcept {
    specbleach_stereo_free(handle);
  }
};

using StereoGroupPtr = std::unique_ptr<specbleach_stereo, StereoGroupDeleter>;

/// Creates a multi-channel engine group (extras); null on failure.
inline StereoGroupPtr make_stereo_group(const uint32_t sample_rate,
                                        const float frame_size_ms,
                                        const uint32_t channels) noexcept {
  return StereoGroupPtr(
      specbleach_stereo_initialize(sample_rate, frame_size_ms, channels));
}

} // namespace specbleach

#endif /* SPECBLEACH_STEREO_HPP_INCLUDED */
