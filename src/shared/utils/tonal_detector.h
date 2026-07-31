/*
libspecbleach - A spectral processing library
*/

#ifndef TONAL_DETECTOR_H
#define TONAL_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

/**
 * Identify tonal peaks in a noise profile using frequency-adaptive detection
 * with parabolic peak interpolation.
 *
 * At low frequencies, the detector uses wider background estimation windows
 * and sideband spread to account for the dense harmonic spacing. Parabolic
 * interpolation corrects for FFT scalloping loss to improve threshold accuracy.
 *
 * @param profile Noise profile (magnitude or power)
 * @param max_profile Maximum captured profile (for stationarity check)
 * @param median_profile Median captured profile (for stationarity check)
 * @param size Spectrum size (real_spectrum_size = fft_size/2 + 1)
 * @param sample_rate Audio sample rate in Hz
 * @param fft_size FFT size used for analysis
 * @param tonal_mask Output mask (0.0–1.0, higher = more tonal)
 */
void detect_tonal_components(const float* profile, const float* max_profile,
                             const float* median_profile, uint32_t size,
                             uint32_t sample_rate, uint32_t fft_size,
                             float* tonal_mask);

/**
 * Extract distinct peak frequencies in Hz from a tonal mask.
 *
 * @param tonal_mask Array of tonal mask values (0.0 to 1.0)
 * @param size Spectrum size (real_spectrum_size)
 * @param sample_rate Audio sample rate in Hz
 * @param fft_size STFT analysis size
 * @param peak_freqs_hz Output array to store detected peak frequencies in Hz
 * @param max_peaks Maximum number of peaks to return
 * @return Number of peak frequencies written into peak_freqs_hz
 */
uint32_t tonal_detector_get_peaks(const float* tonal_mask, uint32_t size,
                                  uint32_t sample_rate, uint32_t fft_size,
                                  float* peak_freqs_hz, uint32_t max_peaks);

/**
 * Extract peak frequencies directly from a noise profile array.
 */
uint32_t tonal_detector_get_peaks_from_profile(
    const float* profile, uint32_t size, uint32_t sample_rate,
    uint32_t fft_size, float* peak_freqs_hz, uint32_t max_peaks);

#endif
