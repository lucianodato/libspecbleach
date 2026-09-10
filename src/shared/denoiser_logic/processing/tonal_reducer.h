/*
libspecbleach - A spectral processing library
*/

#ifndef TONAL_REDUCER_H
#define TONAL_REDUCER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct TonalReducer TonalReducer;

/**
 * Initialize the tonal reducer.
 *
 * Encapsulates tonal detection and the dual-path profile split (broadband vs
 * tonal residual) plus the parallel tonal gain path.
 *
 * @param real_spectrum_size Number of spectral bins (fft_size/2 + 1)
 * @param sample_rate Audio sample rate in Hz
 * @param fft_size FFT size used for analysis
 */
TonalReducer* tonal_reducer_initialize(uint32_t real_spectrum_size,
                                       uint32_t sample_rate, uint32_t fft_size);

void tonal_reducer_free(TonalReducer* self);
void tonal_reducer_reset(TonalReducer* self);

/**
 * Split the noise profile into a broadband and a tonal part (dual-path
 * decoupling). Also detects/publishes the tonal mask exactly as the old
 * tonal_reducer_run did (used downstream by the noise floor manager's
 * dual-path budget, the profile tonal offset and the peak reporter).
 *
 * Broadband output `noise_bb`: at mask bins the profile is replaced by the
 * morphological-opening envelope of the profile (tonal peaks lifted out),
 * blended by mask strength; everywhere else it is a copy of the input.
 * Tonal output `noise_tonal`: the extracted peak residual (N - N_bb).
 *
 * When the tonal path is inactive (tonal_reduction_gain >= ~1.0 or an all
 * zero mask / negligible noise) both outputs degrade safely: `noise_bb`
 * equals the input profile and `noise_tonal` is all zeros, so downstream
 * min() combining is a no-op and behavior matches the legacy single path.
 *
 * @param self               TonalReducer instance
 * @param noise_spectrum     Current noise estimate (morphed profile)
 * @param cv_mask_profile    Coefficient of Variation mask profile
 * @param cv_mask_available  Whether the CV mask profile is available
 * @param tonal_reduction_gain Linear reduction coefficient (0.0–1.0)
 * @param noise_bb           Output broadband noise profile
 * @param noise_tonal        Output tonal noise residual
 */
void tonal_reducer_compute_split(TonalReducer* self,
                                 const float* noise_spectrum,
                                 const float* cv_mask_profile,
                                 bool cv_mask_available,
                                 float tonal_reduction_gain, float* noise_bb,
                                 float* noise_tonal);

/**
 * Compute the parallel tonal gain path: a second Wiener evaluation of the
 * chain's (already smoothed) signal magnitude against the tonal residual.
 * Depth follows the same mask-weighted alpha mapping the legacy alpha boost
 * used (up to ALPHA_MAX_TONAL), so signal harmonics coinciding with hum bins
 * survive when they exceed the tonal noise. Output is stabilized with a
 * light per-bin one-pole state (TONAL_GAIN_STABILIZATION_HOPS).
 *
 * @param slot                Per-chain state slot (0 = active/chain A,
 *                            1 = transition/chain B)
 * @param smoothed_magnitude  Smoothed signal magnitude of the calling chain
 * @param noise_tonal         Tonal noise residual from compute_split
 * @param tonal_reduction_gain Linear reduction coefficient (0.0–1.0)
 * @param gain_tonal          Output tonal gain spectrum (1.0 where no notch)
 */
void tonal_reducer_compute_tonal_gains(TonalReducer* self, uint32_t slot,
                                       const float* smoothed_magnitude,
                                       const float* noise_tonal,
                                       float tonal_reduction_gain,
                                       float* gain_tonal);

/**
 * Legacy (coupled-path) alpha boost: raise alpha toward the reduction-depth
 * alpha at published tonal mask bins, only ever boosting (never reducing).
 * Used by the build with TONAL_DUAL_PATH=0 to reproduce the pre-decoupling
 * behavior for the A/B measurement; not part of the default path.
 */
void tonal_reducer_apply_alpha_boost(TonalReducer* self, float* alpha,
                                     float tonal_reduction_gain);

/**
 * Get the tonal mask from the last run (for downstream use like
 * noise_floor_manager dual-path floor blending).
 *
 * @return Read-only pointer to the tonal mask array, or NULL if not run yet
 */
const float* tonal_reducer_get_mask(const TonalReducer* self);

/**
 * Extract detected peak frequencies in Hz directly from the tonal reducer.
 */
uint32_t tonal_reducer_get_peaks(const TonalReducer* self, float* peak_freqs_hz,
                                 uint32_t max_peaks);

#endif
