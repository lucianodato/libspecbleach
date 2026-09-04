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

#ifndef SPECBLEACH_CONFIGURATIONS_H
#define SPECBLEACH_CONFIGURATIONS_H

#include "utils/spectral_utils.h"
#include <stdbool.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_PIf
#define M_PIf (3.14159265358979323846F)
#endif

/* --------------------------------------------------------------------- */
/* Sections below follow the runtime signal pipeline: STFT front-end ->   */
/* noise estimation -> tonal/transient analysis -> NLM/HPSS processing -> */
/* gain computation -> masking -> gain smoothing -> core plumbing.        */
/* --------------------------------------------------------------------- */

/* --------------------------------------------------------------------- */
/* 1. STFT front-end: frame/overlap/windowing and FFT sizing.             */
/* The 1D/2D smoothing choice lives in SpecbleachSmoothingMode, not here. */
/* --------------------------------------------------------------------- */

// STFT configurations
#define OVERLAP_FACTOR 4
#define INPUT_WINDOW_TYPE HANN_WINDOW
#define OUTPUT_WINDOW_TYPE HANN_WINDOW

// FFT configuration: smallest PFFFT-valid size >= frame, frame centered,
// the rest zero-padded. Padding never affects latency (latency =
// frame_size + NLM look-ahead, constant across smoothing modes).
#define PADDING_CONFIGURATION PAD_TO_VALID_SIZE
#define ZEROPADDING_AMOUNT 0

// Spectral Type
#define SPECTRAL_TYPE POWER_SPECTRUM

/* --------------------------------------------------------------------- */
/* 2. Noise estimation: block counts, adaptive estimators, and profile    */
/* scaling.                                                              */
/* --------------------------------------------------------------------- */

// Noise Estimator Constants
#define MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED 5
#define NUMBER_OF_MEDIAN_SPECTRUM 25
#define NUMBER_OF_MEDIAN_SPECTRUM_MAX 64U
#define MEDIAN_WINDOW_MS (287.5F) // 25 frames * 11.5ms Lukin hop
#define MEDIAN_UPDATE_DECIMATION 8
#define MEDIAN_UPDATE_MS (92.0F) // 8 frames * 11.5ms Lukin hop
#define NOISE_ESTIMATION_INTERPOLATION_THRESHOLD (1e-9F)
#define NOISE_ESTIMATION_SMOOTHING_FACTOR (0.5F)
#define ADAPTIVE_NOISE_FLOOR_SMOOTHING (0.5F)
#define RAYLEIGH_NOISE_CV_REFERENCE (0.523F)
#define RAYLEIGH_CV_REFERENCE RAYLEIGH_NOISE_CV_REFERENCE
#define STANDARD_DEVIATION_MORPH_SCALE (2.0F)

// Adaptive Estimator
#define ESTIMATOR_SILENCE_THRESHOLD (1e-10F) // Roughly -100dB in power

// Martin (2001) Constants
#define MARTIN_SUBWIN_COUNT 8     // Number of sub-windows
#define MARTIN_SUBWIN_LEN 12      // Sub-window length (96/8)
#define MARTIN_SUBWIN_MS (138.0F) // 12 frames * 11.5ms Lukin hop
#define MARTIN_BIAS_CORR 1.5F // Conservative bias correction for min tracking
#define MARTIN_SMOOTH_ALPHA 0.8F // Baseline smoothing for PSD

// SPP-MMSE Estimator Constants
#define SPP_PRIOR_H1 (0.5F)      // P(H1) - Speech present prior
#define SPP_PRIOR_H0 (0.5F)      // P(H0) - Speech absent prior
#define SPP_FIXED_XI_H1 (31.62F) // Fixed a priori SNR (15 dB in linear)
#define SPP_ALPHA_POW (0.8F)     // Power spectrum smoothing factor
#define SPP_SMOOTH_SPP (0.9F)    // SPP smoothing for stagnation control
#define SPP_CURRENT_SPP (0.1F)   // Current SPP weighting for stagnation control
#define SPP_STAGNATION_CAP (0.99F) // Maximum SPP value to prevent locking

// Brandt (Trimmer Mean) Constants
#define BRANDT_DEFAULT_HISTORY_MS 5000.0f
#define BRANDT_DEFAULT_PERCENTILE 0.5f
#define BRANDT_MIN_CONFIDENCE                                                  \
  0.90f // Lowered from 0.98 for better learning speed
#define BRANDT_ESTIMATOR_BIAS_EPSILON                                          \
  (1e-6F) // Precision for bias correction calc
#define BRANDT_ESTIMATOR_MIN_HISTORY_FRAMES                                    \
  5U // Minimum frames for history-based tracking
#define BRANDT_ESTIMATOR_MIN_DURATION_MS                                       \
  (0.1F) // Safety floor for duration calcs
#define BRANDT_ESTIMATOR_STATS_UPDATE_INTERVAL_FRAMES 4U
#define BRANDT_ESTIMATOR_STATS_UPDATE_MS (46.0F) // 4 frames * 11.5ms Lukin hop

// Noise Profile Offset
#define NOISE_PROFILE_OFFSET_DEFAULT_DB 0.0f
#define NOISE_PROFILE_OFFSET_MIN_DB (-12.0F)
#define NOISE_PROFILE_OFFSET_MAX_DB (12.0F)

/* --------------------------------------------------------------------- */
/* 3. Tonal detection and reduction.                                      */
/* --------------------------------------------------------------------- */

// Tonal Detector Constants
#define PEAK_THRESHOLD 2.5f
#define MIN_PEAK_PROMINENCE                                                    \
  1e-7f // Minimum absolute prominence above background
#define TONAL_MEDIAN_FILTER_WINDOW                                             \
  15 // 15-bin frequency-domain median-filter window
#define TONAL_PEAK_MIN_SIGNIFICANCE                                            \
  0.08f // Minimum tonal mask value to report as peak
#define TONAL_PEAK_MIN_LOCAL_PROMINENCE                                        \
  0.01f // Minimum sharpness over immediate neighbors
#define TONAL_PEAK_MIN_WIDER_PROMINENCE                                        \
  0.0f // Minimum sharpness over wider neighborhood (k-2, k+2)
#define TONAL_PEAK_OCTAVE_DYNAMIC_RANGE_DB                                     \
  20.0f // Max dB below strongest peak within same octave band
#define TONAL_PEAK_MIN_OCTAVE_RELATIVE_POWER                                   \
  0.01f // Power ratio within octave (10^(-20/10))
#define TONAL_PEAK_GLOBAL_DYNAMIC_RANGE_DB                                     \
  40.0f // Max dB below global broadband peak energy
#define TONAL_PEAK_MIN_GLOBAL_RELATIVE_POWER                                   \
  1e-4f                              // Power ratio globally (10^(-40/10))
#define MAX_TONAL_PEAKS_REPORTED 32U // Maximum candidate peaks returned
#define TONAL_OCTAVE_LOWER_RATIO 0.7071f
#define TONAL_OCTAVE_UPPER_RATIO 1.4142f
#define TONAL_OCTAVE_MIN_WIDTH_BINS 8
#define TONAL_OCTAVE_FALLBACK_HALF_WIDTH 4
#define TONAL_PEAK_MIN_FREQ_HZ 20.0f
#define TONAL_PEAK_NYQUIST_SAFETY_FACTOR 0.48f
#define TONAL_DETECTOR_DEQUE_CAPACITY 4096U
#define TONAL_REDUCER_NEGLIGIBLE_NOISE_THRESHOLD (1e-12F)
#define ALPHA_MAX_TONAL (10.F)

// Transient Detector Constants
#define UPPER_LIMIT (5.F)
#define DEFAULT_TRANSIENT_THRESHOLD (2.F)
#define MIN_INNOVATION_ENERGY                                                  \
  1e-4F // Minimum clean energy threshold in band (~-40dB magnitude)
#define ONSET_RATIO_SENSITIVITY 0.25F // Innovation required for full weight
#define TRANSIENT_SMOOTH_ALPHA 0.8F   // Reference smoothing alpha
#define TRANSIENT_CLEAN_NOISE_SCALE                                            \
  4.0F // Noise over-subtraction factor for clean transient estimation

/* --------------------------------------------------------------------- */
/* 4. NLM 2D smoothing (Lukin Algorithm B) and HPSS transient protection. */
/* --------------------------------------------------------------------- */

// NLM (Lukin Algorithm B) Parameters
// Time geometry is fixed in ms and recomputed per frame size (see
// frame_rate_norm.h) so temporal smear stays constant across 23-93ms frames.
// Patch size follows Lukin & Todd AES123 (~46ms frame, 75% overlap,
// hop ~11.5ms: PATCH 8 == 92ms). Search context is past-heavy per the paper's
// [-16..+4] mapping (184ms past / 46ms future at the Lukin hop): this cuts
// NLM lookahead latency vs a symmetric context and favors post- over
// pre-echoes of noise. Paste follows the paper's 4x4 block (~80Hz freq
// stride, time resolution is one frame per NLM pass).
#define NLM_PATCH_SIZE 8U
#define NLM_PASTE_BLOCK_SIZE 8U
#define NLM_SEARCH_RANGE_FREQ 8U
#define NLM_SEARCH_RANGE_TIME_PAST 16U
#define NLM_SEARCH_RANGE_TIME_FUTURE 4U
#define NLM_PATCH_TIME_MS (92.0F)
#define NLM_SEARCH_PAST_MS (184.0F)
#define NLM_SEARCH_FUTURE_MS (46.0F)
#define NLM_SEARCH_FREQ_HZ (170.0F)
#define NLM_PASTE_FREQ_HZ (80.0F)
// DFTT-lite post-NLM cleanup (Lukin & Todd AES123 S4.2, lite adaptation:
// the paper uses 32x16 blocks with 8/4-bin hops at a fixed STFT size; we use
// half-width freq blocks at DFTT_BLOCK_FREQ_HZ with hop block/2 and a
// past-only time span of DFTT_TIME_MS so no extra latency is introduced and
// the analysis stays invariant across frame sizes. Dual input per the paper:
// the noisy SNR tile undergoes analysis/modification/synthesis while the
// NLM-smoothed tile sets the per-coefficient suppression threshold under a
// white-quefrency assumption).
#define DFTT_BLOCK_FREQ_HZ (340.0F)
#define DFTT_MIN_BLOCK_FREQ (8U)
#define DFTT_MAX_BLOCK_FREQ (32U)
#define DFTT_TIME_MS (92.0F)
#define DFTT_MIN_TIME_FRAMES (4U)
#define DFTT_MAX_TIME_FRAMES (16U)
#define DFTT_SILENCE_EPS (1e-9F)
#define DFTT_KILL_K                                                            \
  1.0F // Kill threshold as a multiple of the tile's diffuse floor
       // energy: quefrency bins below it are treated as diffuse speckle
#define DFTT_VETO_GATE                                                         \
  2.0F // NLM bins at/above this are confident speech: the veto skips them
       // and keeps NLM's exact values (map floor sits near 1.0)
// Vectorized-distance patch ceiling (matches the frame_rate_norm.h fallback
// clamp 4..16) and pointer-cache halo (half of the ceiling each side).
#define NLM_MAX_PATCH_FRAMES 16U
#define NLM_HALO_FRAMES 8U
// Reference hops: legacy 50ms-frame tuning anchor (12.5ms hop) for tau
// conversion of per-frame IIRs; Lukin hop (~11.5ms) for NLM geometry.
#define LEGACY_REF_HOP_SEC (0.0125F)
#define LUKIN_REF_HOP_SEC (0.0115F)
#define NLM_DEFAULT_H_PARAMETER 1.0F
#define NLM_MAX_H_PARAMETER 5.0F
#define NLM_MIN_WEIGHT 1e-10F
#define NLM_SNR_NOISE_FLOOR_MIN 1e-9F
#define NLM_DISTANCE_THRESHOLD_MULTIPLIER 4.0F
#define NLM_NUM_THREADS_DEFAULT 4U
#define NLM_MAX_THREADS 16U

// Must be >= search_time_past + search_time_future + patch_size for NLM caching
// Using power-of-two (64U) for efficient modulo wrap-around and future headroom
#define DELAY_BUFFER_FRAMES 64U

#define HPSS_SLIDING_ITERATIONS 3U
#define HPSS_SLIDING_SMOOTH_FACTOR 0.50F

/* --------------------------------------------------------------------- */
/* 5. Gain computation: spectral subtraction exponents and                */
/* oversubtraction limits.                                                */
/* --------------------------------------------------------------------- */

// Gain Estimators
#define GSS_EXPONENT                                                           \
  2.0F // 2 Power Subtraction / 1 Magnitude Subtraction / 0.5 Spectral
       // Subtraction

// Oversubtraction criteria
#define ALPHA_MAX (4.F)
#define ALPHA_MIN (1.F)
#define DEFAULT_OVERSUBTRACTION (ALPHA_MIN)
#define DEFAULT_UNDERSUBTRACTION (0.01F)
#define SUPPRESSION_LOWER_SNR_DB (-5.0F)
#define SUPPRESSION_HIGHER_SNR_DB (20.0F)

/* --------------------------------------------------------------------- */
/* 6. Masking and psychoacoustics.                                        */
/* --------------------------------------------------------------------- */

// Noise Scaling strategy
#define CRITICAL_BANDS_TYPE BARK_SCALE
#define GAIN_ESTIMATION_TYPE WIENER

// Masking Veto defaults
#define USE_TEMPORAL_MASKING_DEFAULT true

// Absolute hearing thresholds
#define REFERENCE_SINE_WAVE_FREQ (1000.F)
#define REFERENCE_LEVEL (90.F)
#define SINE_AMPLITUDE (1.F)

// Johnston Psychoacoustic Model (1988)
#define DB_FS_TO_SPL_REF 96.0F    // 0dBFS -> 96dB SPL reference level
#define POWER_LAW_EXPONENT (0.6F) // Johnston power-law integration exponent

// Default Additivity Exponents
#define SPECTRAL_ADDITIVITY_EXPONENT_STANDARD (1.0F) // Pure Johnston
#define SPECTRAL_ADDITIVITY_EXPONENT_PEAQ (0.4F)     // Advanced precision

// Schroeder spreading function constants
#define S_MIN_UPWARD (5.0F)  // Minimum upward spreading slope (dB/Bark)
#define S_MAX_UPWARD (15.0F) // Maximum upward spreading slope (dB/Bark)
#define S_DOWNWARD (25.0F)   // Constant downward spreading slope (dB/Bark)

// Threshold Offsets (NMT = Noise-Masking-Tone, TMN = Tone-Masking-Noise)
#define NMT_OFFSET_DB (6.0F)    // Standard offset for NMT
#define TMN_OFFSET_BASE (14.5F) // Base offset for TMN (Bark-dependent)

// Johnston SFM (Spectral Flatness Measure) constants
#define SFM_MIN_DB (-60.0F) // Minimum expected SFM (highly tonal)
#define SFM_MAX_DB (0.0F)   // Maximum expected SFM (random noise)
// Temporal Masking Constants (seconds; the _SEC suffix is load-bearing:
// hop_sec is seconds, so decays are expf(-hop_sec / tau_sec))
// Frequency-dependent forward masking (Low: 100ms, High: 25ms)
#define FORWARD_MASKING_TAU_LOW_SEC (0.100F)  // 100ms decay for low frequencies
#define FORWARD_MASKING_TAU_HIGH_SEC (0.025F) // 25ms decay for high frequencies
#define BACKWARD_MASKING_TAU_SEC (0.010F)     // 10ms decay for pre-masking

// Schroeder Slope Adaptation Constants
#define S_LEVEL_REF_DB 40.0F // Reference level for slope adaptation (dB SPL)
#define S_SLOPE_FACTOR 0.2F  // Slope broadening factor per dB

// Veto Parameters
#define MASKING_VETO_SMOOTHING 0.5F  // Stabilization alpha for clean signal
#define MASKING_VETO_NMR_RANGE 20.0F // NMR range for protection mapping (dB)
#define MASKING_VETO_NOISE_GATE                                                \
  1.5F // Noise floor multiplier for clean signal estimation (gates noise
       // variance)
#define MASKING_VETO_SNR_THRESHOLD                                             \
  3.0F // ~5dB SNR local bin-level threshold for full protection

/* --------------------------------------------------------------------- */
/* 7. Gain smoothing envelopes and runtime smoothing-mode crossfade.      */
/* --------------------------------------------------------------------- */

#define GAIN_SMOOTHING_MIN_RELEASE_SEC (0.010F)
#define GAIN_SMOOTHING_MAX_RELEASE_SEC (0.350F)
#define GAIN_SMOOTHING_MIN_ATTACK_SEC (0.001F)
#define GAIN_SMOOTHING_MAX_ATTACK_SEC (0.040F)

// Duration of the internal crossfade when switching smoothing modes at runtime
#define SMOOTHING_TRANSITION_SECONDS (0.030F)
#define SMOOTHING_TRANSITION_MIN_FRAMES 4U

/* --------------------------------------------------------------------- */
/* 8. Core plumbing: numeric floors and circular-buffer capacity.         */
/* --------------------------------------------------------------------- */

#define SPECTRAL_EPSILON (1e-12F)
#define MAX_SPECTRAL_CIRCULAR_BUFFER_LAYERS 8

/* --------------------------------------------------------------------- */
/* Compile-time assertions for configuration validity                     */
/* --------------------------------------------------------------------- */
_Static_assert(HANN_WINDOW >= 0 && HANN_WINDOW <= 3,
               "HANN_WINDOW must be between 0 and 3");
_Static_assert(HAMMING_WINDOW >= 0 && HAMMING_WINDOW <= 3,
               "HAMMING_WINDOW must be between 0 and 3");
_Static_assert(BLACKMAN_WINDOW >= 0 && BLACKMAN_WINDOW <= 3,
               "BLACKMAN_WINDOW must be between 0 and 3");
_Static_assert(VORBIS_WINDOW >= 0 && VORBIS_WINDOW <= 3,
               "VORBIS_WINDOW must be between 0 and 3");

// Additional C17 compile-time validations
_Static_assert(sizeof(uint32_t) == 4, "uint32_t must be exactly 32 bits");

#endif // SPECBLEACH_CONFIGURATIONS_H
