/*
libspecbleach - A spectral processing library
*/

#ifndef TONAL_DETECTOR_H
#define TONAL_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

/**
 * Identify tonal peaks in a noise profile.
 * @param profile Noise profile (magnitude or power)
 * @param max_profile Maximum captured profile (for stationarity check)
 * @param median_profile Median captured profile (for stationarity check)
 * @param size Spectrum size
 * @param sensitivity User sensitivity (0.0 to 1.0)
 * @param tonal_mask Output mask (1.0 for tonal, 0.0 for broadband, or
 * intermediate)
 */
void detect_tonal_components(const float* profile, const float* max_profile,
                             const float* median_profile, uint32_t size,
                             float* tonal_mask);

#endif
