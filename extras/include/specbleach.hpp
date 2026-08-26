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
 * @brief Header-only C++ ownership wrappers for the libspecbleach C API.
 *
 * This header is deliberately framework-agnostic (standard library only):
 * it serves JUCE plugin authors and raw-VST3/CLAP/LV2 developers,
 * DAW/editor codebases, and standalone audio applications equally. It
 * wraps OWNERSHIP only — RAII lifetime management via std::unique_ptr —
 * and adds nothing else: no threading policy, no buffers, no framework
 * types. All processing calls remain the plain C functions, which are
 * safe to call directly on the owned handles.
 *
 * Usage:
 * @code
 *   #include <specbleach.hpp>
 *
 *   auto denoiser = specbleach::make_denoiser(48000U, 46.0F);
 *   if (!denoiser) { / * handle allocation failure * / }
 *   specbleach_denoiser_process(denoiser.get(), nframes, in, out);
 *   // freed automatically, exception-safe, non-copyable by construction
 * @endcode
 *
 * The same shape exists for every handle type in the library:
 * make_denoiser, make_2d_denoiser, make_stereo_group, and
 * make_delay_line (extras). Handles are movable, never copyable, so a
 * double-free cannot be expressed.
 */

#ifndef SPECBLEACH_HPP
#define SPECBLEACH_HPP

#include <memory>

#include "specbleach_2d_denoiser.h"
#include "specbleach_delay_line.h"
#include "specbleach_denoiser.h"
#include "specbleach_stereo.h"

namespace specbleach {

struct DenoiserDeleter {
  /**
   * @brief Releases a denoiser handle.
   *
   * @param handle Denoiser handle to release.
   */
  void operator()(specbleach_denoiser* handle) const noexcept {
    specbleach_denoiser_free(handle);
  }
};

struct Denoiser2dDeleter {
  /**
   * @brief Releases a 2D denoiser handle.
   *
   * @param handle The 2D denoiser handle to release.
   */
  void operator()(specbleach_2d_denoiser* handle) const noexcept {
    specbleach_2d_free(handle);
  }
};

struct StereoGroupDeleter {
  /**
   * @brief Releases a stereo processing handle.
   *
   * @param handle Stereo processing handle to release.
   */
  void operator()(specbleach_stereo* handle) const noexcept {
    specbleach_stereo_free(handle);
  }
};

struct DelayLineDeleter {
  void operator()(specbleach_delay_line* handle) const noexcept {
    specbleach_delay_line_free(handle);
  }
};

using DenoiserPtr = std::unique_ptr<specbleach_denoiser, DenoiserDeleter>;
using Denoiser2dPtr =
    std::unique_ptr<specbleach_2d_denoiser, Denoiser2dDeleter>;
using StereoGroupPtr = std::unique_ptr<specbleach_stereo, StereoGroupDeleter>;
using DelayLinePtr = std::unique_ptr<specbleach_delay_line, DelayLineDeleter>;

/**
 * @brief Creates a single-channel spectral denoiser.
 *
 * @return An owning denoiser pointer, or null if initialization fails.
 */
inline DenoiserPtr make_denoiser(const uint32_t sample_rate,
                                 const float frame_size_ms) noexcept {
  return DenoiserPtr(
      specbleach_denoiser_initialize(sample_rate, frame_size_ms));
}

/**
 * @brief Creates a single-channel 2D non-local means denoiser.
 *
 * @param sample_rate Audio sample rate in hertz.
 * @param frame_size_ms Processing frame size in milliseconds.
 * @return An owning pointer to the denoiser, or a null pointer if initialization fails.
 */
inline Denoiser2dPtr make_2d_denoiser(const uint32_t sample_rate,
                                      const float frame_size_ms) noexcept {
  return Denoiser2dPtr(specbleach_2d_initialize(sample_rate, frame_size_ms));
}

/**
 * @brief Creates an owned multi-channel stereo processing group.
 *
 * @param sample_rate Audio sample rate in hertz.
 * @param frame_size_ms Processing frame size in milliseconds.
 * @param channels Number of audio channels.
 * @param engine Stereo processing engine to use.
 * @return Stereo group handle, or a null pointer if initialization fails.
 */
inline StereoGroupPtr make_stereo_group(
    const uint32_t sample_rate, const float frame_size_ms,
    const uint32_t channels, const SpecbleachStereoEngine engine) noexcept {
  return StereoGroupPtr(specbleach_stereo_initialize(sample_rate, frame_size_ms,
                                                     channels, engine));
}

/**
 * @brief Creates an owned delay line for the specified channel count.
 *
 * @param max_delay_samples Maximum delay capacity in samples.
 * @param channels Number of audio channels.
 * @return An owning delay-line pointer, or a null pointer if initialization fails.
 */
inline DelayLinePtr make_delay_line(const uint32_t max_delay_samples,
                                    const uint32_t channels) noexcept {
  return DelayLinePtr(
      specbleach_delay_line_initialize(max_delay_samples, channels));
}

} // namespace specbleach

#endif /* SPECBLEACH_HPP */
