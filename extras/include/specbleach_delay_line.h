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
 * @file specbleach_delay_line.h
 * @brief Multi-channel single-tap delay line for time-aligning processor
 *        streams (the structural-alignment recipe's building block).
 *
 * Why this exists: libspecbleach engine families report different
 * algorithmic latencies. The documented strategy for artifact-free live
 * switching between them (examples/README.md, "Switching between
 * processors", Recipe A) keeps BOTH families running permanently aligned
 * at max(latency_a, latency_b): the shorter-latency family renders through
 * this delay line, the application reports the maximum latency once, and
 * switching becomes a plain equal-power crossfade between two streams that
 * never shift relative to each other.
 *
 * Contract notes:
 * - This is deliberately ONLY a delay line. It owns no policy: what delay
 *   to use, when to change it, and what you report to your host remain
 *   application decisions.
 * - process() is real-time safe: no allocation, no locks, no I/O.
 *   set_delay() is a plain member write and is safe between process()
 *   calls; do not move the tap while a block is being processed from
 *   another thread.
 * - Output semantics: out[n] == in[n - delay]. With delay == 0 the line is
 *   a pass-through. The first `delay` samples of output are silence.
 * - Any combination of block sizes is supported across calls; the internal
 *   ring keeps full stream continuity.
 */

#ifndef SPECBLEACH_DELAY_LINE_H
#define SPECBLEACH_DELAY_LINE_H

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_export.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct specbleach_delay_line specbleach_delay_line;

/**
 * Creates a delay line. Allocation happens here only.
 *
 * @param max_delay_samples Upper bound for set_delay(). May be 0 (the line
 *                          is then permanently a pass-through).
 * @param channels Number of interleaved-in-time channel pointer slots.
 * @return NULL on invalid arguments or allocation failure.
 */
SPECBLEACH_API specbleach_delay_line* specbleach_delay_line_initialize(
    uint32_t max_delay_samples, uint32_t channels);

/**
 * Frees a delay line. NULL is accepted and ignored.
 */
SPECBLEACH_API void specbleach_delay_line_free(specbleach_delay_line* instance);

/**
 * Sets the delay in samples applied by subsequent process() calls.
 * Values above the configured maximum return false and leave the current
 * delay unchanged.
 */
SPECBLEACH_API bool specbleach_delay_line_set_delay(
    specbleach_delay_line* instance, uint32_t delay_samples);

/**
 * Current tap position in samples.
 */
SPECBLEACH_API uint32_t
specbleach_delay_line_get_delay(specbleach_delay_line* instance);

/**
 * Configured upper bound passed at initialization.
 */
SPECBLEACH_API uint32_t
specbleach_delay_line_get_max_delay(specbleach_delay_line* instance);

/**
 * Number of channel pointers expected by process().
 */
SPECBLEACH_API uint32_t
specbleach_delay_line_get_channels(specbleach_delay_line* instance);

/**
 * Pushes a block in and writes the delayed block out. In-place operation
 * (out aliases in per channel) is allowed and well-defined. Real-time
 * safe.
 */
SPECBLEACH_API bool specbleach_delay_line_process(
    specbleach_delay_line* instance, uint32_t num_samples,
    const float* const* in, float* const* out);

#ifdef __cplusplus
}
#endif

#endif /* SPECBLEACH_DELAY_LINE_H */
