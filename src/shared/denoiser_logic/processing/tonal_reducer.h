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
 * Encapsulates tonal detection and alpha boosting. Must run BEFORE the
 * masking veto so the veto can naturally protect signal harmonics.
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
 * Detect tonal components and boost alpha at tonal bins.
 *
 * Internally calls the tonal detector on the noise spectrum, then boosts
 * alpha proportionally so the Wiener filter suppresses tonal noise by
 * the requested dB amount.
 *
 * @param self             TonalReducer instance
 * @param noise_spectrum   Current noise estimate (morphed profile)
 * @param cv_mask_profile      Coefficient of Variation mask profile
 * @param cv_mask_available    Whether the CV mask profile is available
 * @param alpha                Per-bin oversubtraction array (modified in place)
 * @param tonal_reduction_gain Linear reduction coefficient (0.0–1.0)
 */
void tonal_reducer_run(TonalReducer* self, const float* noise_spectrum,
                       const float* cv_mask_profile, bool cv_mask_available,
                       float* alpha, float tonal_reduction_gain);

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
