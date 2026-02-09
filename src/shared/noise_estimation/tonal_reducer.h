/*
libspecbleach - A spectral processing library
*/

#ifndef TONAL_REDUCER_H
#define TONAL_REDUCER_H

#include <stdint.h>

/**
 * Boost the oversubtraction factor (alpha) at tonal bins so each detected
 * tonal component gets the requested dB of reduction from its own level.
 *
 * This must run BEFORE the masking veto so the veto can naturally reduce
 * alpha when wanted signal harmonics coincide with tonal noise.
 *
 * @param alpha          Per-bin oversubtraction array (modified in place)
 * @param tonal_mask     Per-bin tonal mask (0.0–1.0 from tonal_detector)
 * @param real_spectrum_size  Number of spectral bins (fft_size/2 + 1)
 * @param tonal_reduction_amount  Linear reduction coefficient (0.0–1.0,
 *                                from dB conversion: 0.004 ≈ -48 dB)
 */
void tonal_reducer_apply(float* alpha, const float* tonal_mask,
                         uint32_t real_spectrum_size,
                         float tonal_reduction_amount);

#endif
