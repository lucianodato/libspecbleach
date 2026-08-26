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

#ifndef SPECBLEACH_TRANSITION_H_INCLUDED
#define SPECBLEACH_TRANSITION_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

#include "specbleach_export.h"

/**
 * Equal-power crossfade between two already-configured engine groups.
 *
 * This module belongs to the optional extras layer. It is engine-agnostic:
 * it knows nothing about denoisers or parameters, it only blends two
 * streams of audio. The caller renders BOTH engines' wet outputs itself
 * and feeds them here every block during a transition.
 *
 * Behavior:
 * - Fades are equal-power (sin/cos) over SPECBLEACH_TRANSITION_FADE_TIME_MS,
 *   so loudness stays constant mid-fade (a linear fade would dip -3 dB).
 * - If the two engines report different latencies, whichever side has the
 *   SHORTER latency is read through an internal alignment tap so both blend
 *   sides share one time origin, and a slew phase slides the tap back out
 *   after landing on the shorter engine.
 * - Callers must update host delay compensation when they call begin():
 *   specbleach_transition_get_latency() reports the target latency from
 *   that moment on.
 *
 * KNOWN LIMITATION — latency-mismatched engines in hosted plugins:
 * changing the reported latency mid-stream forces the host to re-anchor its
 * delay compensation, and that host-side event cannot be softened from
 * inside the plugin. Empirically (DAWs such as REAPER/Ardour/Bitwig) the
 * result is audible skips or gaps regardless of how clean the internal
 * blend is. If your engines have different algorithmic latencies and you
 * need seamless switching, do what Noise Repellent does: run the
 * shorter-latency family through a permanent delay ring so both families
 * share one time origin, report max(latency) constantly, and use plain
 * weighted blending between the aligned streams. This module then reduces
 * to that final blend step, which it performs correctly.
 *
 * Best suited to: engines with EQUAL latency, non-hosted processing
 * (offline files, custom players without PDC), or as the blend stage of a
 * structurally aligned pipeline.
 */

/**
 * Duration of the equal-power fade and of the post-transition delay
 * ramp-down, in milliseconds.
 */
#define SPECBLEACH_TRANSITION_FADE_TIME_MS 40.0f

/**
 * Opaque handle to a transition state machine.
 */
typedef struct specbleach_transition // NOLINT(readability-identifier-naming)
    specbleach_transition;

/**
 * Creates a transition processor.
 *
 * @param max_block_size Largest block that will ever be passed to
 * specbleach_transition_process() or specbleach_transition_feed().
 * @param max_alignment_delay Largest latency difference (in samples) the
 * instance must be able to align; begin() fails if asked for more.
 * @return New instance or NULL on failure.
 */
SPECBLEACH_API specbleach_transition* specbleach_transition_initialize(
    uint32_t sample_rate, uint32_t max_block_size, uint32_t channels,
    uint32_t max_alignment_delay);

/**
 * Frees the instance. NULL is a no-op.
 *
 * Not real-time safe.
 */
SPECBLEACH_API void specbleach_transition_free(specbleach_transition* instance);

/**
 * Begins a fade from one engine group to another.
 *
 * Call whenever parameters/engines change behind the scenes. May be called
 * again mid-fade; the fade restarts toward the new target.
 *
 * Not real-time safe (may allocate the alignment delay line once).
 *
 * @param latency_from Algorithmic latency of the source outputs.
 * @param latency_to Algorithmic latency of the target outputs.
 */
SPECBLEACH_API bool specbleach_transition_begin(specbleach_transition* instance,
                                                uint32_t latency_from,
                                                uint32_t latency_to);

/**
 * Blends one block.
 *
 * @param from_out Wet output of the source engine group. Must be non-NULL
 * while a fade is running; may be NULL afterwards (pass-through of to_out,
 * with delay ramp-down applied if required).
 * @param to_out Wet output of the target engine group. Must not be NULL.
 * @param blended Destination buffer; may alias to_out.
 *
 * Real-time safe once specbleach_transition_begin() has completed.
 * @return false on invalid arguments.
 */
SPECBLEACH_API bool specbleach_transition_process(
    specbleach_transition* instance, uint32_t number_of_samples,
    const float** from_out, const float** to_out, float** blended);

/**
 * Returns true while a fade or a delay ramp-down is in progress.
 */
SPECBLEACH_API bool specbleach_transition_active(
    const specbleach_transition* instance);

/**
 * Returns the latency hosts should report from this moment on: the target
 * latency once begin() was called, 0 before. Callers must update host delay
 * compensation when they call begin().
 */
SPECBLEACH_API uint32_t
specbleach_transition_get_latency(const specbleach_transition* instance);

#ifdef __cplusplus
}
#endif
#endif /* SPECBLEACH_TRANSITION_H_INCLUDED */
