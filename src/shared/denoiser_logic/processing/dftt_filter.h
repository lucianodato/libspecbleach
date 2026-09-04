/*
libspecbleach - A spectral processing library

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

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

#ifndef DFTT_FILTER_H
#define DFTT_FILTER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct DfttFilter DfttFilter;

/**
 * DFT-thresholding post-filter for the 2D NLM chain.
 * Lukin & Todd AES123 §4.2: overlapping Hann tiles of the SNR map undergo a
 * 2D DFT; per-coefficient Wiener gain with the threshold set from the
 * speckle residual between the noisy and the NLM-smoothed tile (the NLM
 * tile acts as the structure prior; speckle is white in the tile-DFT
 * domain). Then inverse DFT and overlap-add. Tiles are past-only in time,
 * so no extra latency is introduced. All memory is allocated at init;
 * processing never allocates, locks, or performs I/O.
 */

/**
 * Initialize a DFTT filter.
 * @param spectrum_size Number of frequency bins (real spectrum size)
 * @param time_span_frames Time extent of a tile in frames (clamped to
 *                         DFTT_MAX_TIME_FRAMES); 0 fails
 * @param block_freq_bins Frequency extent of a tile in bins, hop is
 *                        DFTT_FREQ_OVERLAP-divided (clamped to
 *                        DFTT_MAX_DIM); 0 fails
 * @return Pointer to initialized filter, or NULL on failure
 */
DfttFilter* dftt_filter_initialize(uint32_t spectrum_size,
                                   uint32_t time_span_frames,
                                   uint32_t block_freq_bins);

/**
 * Free resources associated with a DFTT filter.
 * @param filter Pointer to filter to free (NULL-safe)
 */
void dftt_filter_free(DfttFilter* filter);

/**
 * Push one frame pair into the history rings.
 * @param filter Pointer to filter
 * @param noisy_snr Current noisy SNR frame (size: spectrum_size)
 * @param smoothed_snr Current NLM-smoothed SNR frame (size: spectrum_size)
 */
void dftt_filter_push(DfttFilter* filter, const float* noisy_snr,
                      const float* smoothed_snr);

/**
 * Process the current frame pair and compute the refined SNR output.
 * @param filter Pointer to filter
 * @param refined_snr Output array for refined SNR (size: spectrum_size)
 * @return true if processing succeeded, false if history not yet full or no
 *         fresh frame was pushed since the last call
 */
bool dftt_filter_process(DfttFilter* filter, float* refined_snr);

/**
 * Check if the filter history is ready for processing.
 * @param filter Pointer to filter
 * @return true if enough frames were pushed
 */
bool dftt_filter_is_ready(DfttFilter* filter);

/**
 * Reset the filter state, clearing all buffered frames.
 * @param filter Pointer to filter (NULL-safe)
 */
void dftt_filter_reset(DfttFilter* filter);

#endif /* DFTT_FILTER_H */
