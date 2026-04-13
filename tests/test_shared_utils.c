/*
 * Unit tests for Pre-estimation modules
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/utils/absolute_hearing_thresholds.h"
#include "shared/utils/critical_bands.h"
#include "shared/utils/masking_estimator.h"
#include "shared/utils/simd_utils.h"
#include "shared/utils/spectral_smoother.h"
#include "shared/utils/transient_detector.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_absolute_hearing_thresholds(void) {
  printf("Testing Absolute Hearing Thresholds...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  AbsoluteHearingThresholds* aht = absolute_hearing_thresholds_initialize(
      sample_rate, fft_size, POWER_SPECTRUM);

  TEST_ASSERT(aht != NULL, "AHT initialization should succeed");

  float spectrum[513] = {0.0f}; // 513 = fft_size/2 + 1
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 0.001f; // Set low values that should be floored
  }

  TEST_ASSERT(apply_thresholds_as_floor(aht, spectrum),
              "Apply thresholds should succeed");

  // Check that some values were floored to hearing thresholds
  bool some_values_changed = false;
  for (int i = 0; i < 513; i++) {
    if (spectrum[i] > 0.001f) {
      some_values_changed = true;
      break;
    }
  }
  TEST_ASSERT(some_values_changed,
              "Some spectrum values should be floored to hearing thresholds");

  absolute_hearing_thresholds_free(aht);
  printf("✓ Absolute Hearing Thresholds tests passed\n");
}

void test_critical_bands(void) {
  printf("Testing Critical Bands (Shared Utils)...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  CriticalBands* cb =
      critical_bands_initialize(sample_rate, fft_size, BARK_SCALE);
  TEST_ASSERT(cb != NULL, "Critical bands initialization should succeed");

  uint32_t num_bands = get_number_of_critical_bands(cb);
  TEST_ASSERT(num_bands > 0, "Should have some critical bands");

  // Test getting band indexes
  for (uint32_t i = 0; i < num_bands; i++) {
    CriticalBandIndexes indexes = get_band_indexes(cb, i);
    TEST_ASSERT(indexes.start_position < indexes.end_position,
                "Band start should be before end");
    TEST_ASSERT(indexes.end_position <= (fft_size / 2U) + 1U,
                "Band end should be within spectrum bounds");
  }

  // Test computing critical bands spectrum
  float spectrum[513] = {0.0f};
  float* critical_bands_spectrum = (float*)calloc(num_bands, sizeof(float));
  TEST_ASSERT(critical_bands_spectrum != NULL, "Allocation failed");

  for (int i = 0; i < 513; i++) {
    spectrum[i] = (float)i * 0.1f;
  }

  TEST_ASSERT(
      compute_critical_bands_spectrum(cb, spectrum, critical_bands_spectrum),
      "Compute critical bands should succeed");

  free(critical_bands_spectrum);
  critical_bands_free(cb);
  printf("✓ Critical Bands tests passed\n");
}

void test_masking_estimator(void) {
  printf("Testing Masking Estimator...\n");

  uint32_t fft_size = 1024;
  uint32_t sample_rate = 44100;

  MaskingEstimator* me = masking_estimation_initialize(
      fft_size, sample_rate, OPUS_SCALE, POWER_SPECTRUM, true, true);
  TEST_ASSERT(me != NULL, "Masking estimator initialization should succeed");

  float spectrum[513] = {0.0f};
  float masking_thresholds[513] = {0.0f};

  // Create a test spectrum with some energy
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 0.1f + ((float)i * 0.01f);
  }

  TEST_ASSERT(
      compute_masking_thresholds(me, spectrum, spectrum, masking_thresholds),
      "Compute masking thresholds should succeed");

  // Check that masking thresholds are reasonable
  for (int i = 0; i < 513; i++) {
    TEST_ASSERT(masking_thresholds[i] >= 0.0f,
                "Masking thresholds should be non-negative");
  }

  masking_estimation_free(me);
  printf("✓ Masking Estimator tests passed\n");
}

void test_spectral_smoother(void) {
  printf("Testing Spectral Smoother...\n");

  uint32_t fft_size = 1024;

  // Test all smoothing types
  for (int type = NO_SMOOTHING; type <= TRANSIENT_AWARE; type++) {
    SpectralSmoother* ss =
        spectral_smoothing_initialize(fft_size, 44100, (TimeSmoothingType)type);
    TEST_ASSERT(ss != NULL, "Spectral smoother initialization should succeed");

    float spectrum[513] = {0.0f};
    for (int i = 0; i < 513; i++) {
      spectrum[i] = 1.0f + (0.5f * sinf((float)i * 0.1f));
    }

    TimeSmoothingParameters params = {.smoothing = 0.8f};
    TEST_ASSERT(spectral_smoothing_run(ss, params, spectrum),
                "Spectral smoothing should succeed");

    // Run again to test previous spectrum logic
    TEST_ASSERT(spectral_smoothing_run(ss, params, spectrum),
                "Spectral smoothing should succeed on second run");

    // Check that output is reasonable
    for (int i = 0; i < 513; i++) {
      TEST_ASSERT(spectrum[i] >= 0.0f,
                  "Smoothed spectrum should be non-negative");
    }

    spectral_smoothing_free(ss);
  }
  printf("✓ Spectral Smoother tests passed\n");
}

void test_transient_detector(void) {
  printf("Testing Transient Detector...\n");

  uint32_t fft_size = 1024;
  uint32_t real_size = (fft_size / 2) + 1;

  TransientDetector* td = transient_detector_initialize(real_size);
  TEST_ASSERT(td != NULL, "Transient detector initialization should succeed");

  float spectrum[513] = {0.0f};
  float different_spectrum[513] = {0.0f};

  // Create initial spectrum
  for (int i = 0; i < 513; i++) {
    spectrum[i] = 1.0f;
    different_spectrum[i] = 2.0f; // Different spectrum to create change
  }

  // First run - should not detect transient (no previous data, initialized with
  // current)
  bool result1 = transient_detector_process(td, spectrum, NULL);

  // Second run with different spectrum - likely transient due to jump from 1.0
  // to 2.0 Ratio = 2.0 / 1.0 = 2.0. Weight = (2-1)/0.25 = 4.0 -> Clamped
  // to 1.0. Should return true.
  bool result2 = transient_detector_process(td, different_spectrum, NULL);

  // The function should run without error (return value is boolean indicating
  // transient presence) We just test that it doesn't crash and returns a valid
  // boolean
  TEST_ASSERT(result1 == false, "First run should not detect transient (init)");
  TEST_ASSERT(result2 == true, "Second run should detect transient (jump)");

  transient_detector_free(td);
  printf("✓ Transient Detector tests passed\n");
}

void test_simd_utils(void) {
  printf("Testing SIMD Utils...\n");

  sb_simd_state_t state = sb_simd_enable_ftz_daz();

  float buf[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
  float out[8] = {0.0f};

  sb_vec8_t v1 = sb_load8(buf);
  sb_vec8_t v2 = sb_set8(2.0f);

  sb_vec8_t vadd = sb_add8(v1, v2);
  sb_vec8_t vsub = sb_sub8(v1, v2);
  sb_vec8_t vmul = sb_mul8(v1, v2);
  sb_vec8_t vdiv = sb_div8(v1, v2);
  sb_vec8_t vmax = sb_max8(v1, v2);

  sb_store8(out, vadd);
  TEST_FLOAT_CLOSE(out[0], 3.0f, 0.001f);

  sb_vec4_t v4_1 = sb_load4(buf);
  sb_vec4_t v4_2 = sb_set4(2.0f);
  sb_vec4_t v4_add = sb_add4(v4_1, v4_2);
  sb_vec4_t v4_sub = sb_sub4(v4_1, v4_2);
  sb_vec4_t v4_mul = sb_mul4(v4_1, v4_2);

  sb_store4(out, v4_add);
  TEST_FLOAT_CLOSE(out[0], 3.0f, 0.001f);

  float ssd4 = sb_vec4_ssd(v4_1, v4_2);
  (void)ssd4;

  sb_vec8_t mask = sb_gt8(v1, v2);
  sb_vec8_t result = sb_sel8(mask, vsub, vmul);

  sb_store8(out, result);

  sb_acc8_t acc = sb_acc8_zero();
  acc = sb_acc8_add_ssd(acc, v1, v2);
  float hsum = sb_acc8_hsum(acc);
  (void)hsum;

  sb_simd_restore_state(state);
  printf("✓ SIMD Utils tests passed\n");
}

void test_fast_exp(void) {
  printf("Testing fast exponential approximation...\n");

  // x = 0 -> exp(0) = 1
  TEST_FLOAT_CLOSE(sb_fast_expf(0.0f), 1.0f, 0.05f);

  // x = -1.0 -> exp(-1) approx 0.3678
  TEST_FLOAT_CLOSE(sb_fast_expf(-1.0f), expf(-1.0f), 0.05f);

  // x = -5.0 -> exp(-5) approx 0.0067
  TEST_FLOAT_CLOSE(sb_fast_expf(-5.0f), expf(-5.0f), 0.05f);

  // Underflow protection: x < -20 should be 0
  TEST_ASSERT(sb_fast_expf(-21.0f) == 0.0f, "Should underflow to 0.0f");
  TEST_ASSERT(sb_fast_expf(-100.0f) == 0.0f, "Should underflow to 0.0f");

  printf("✓ fast exponential tests passed\n");
}

void test_patch_ssd(void) {
  printf("Testing 8x8 patch SSD...\n");

  float target[64];
  float candidate[64];
  float* cand_rows[8];

  for (int i = 0; i < 64; i++) {
    target[i] = (float)i;
    candidate[i] = (float)i + 1.0f; // Difference of 1.0 for each element
  }

  sb_vec8_t target_vecs[8];
  for (int r = 0; r < 8; r++) {
    target_vecs[r] = sb_load8(target + r * 8);
    cand_rows[r] = candidate + r * 8;
  }

  float ssd = sb_vec8_patch_ssd(target_vecs, cand_rows);

  // Expected SSD: 64 elements, each with diff = 1.0 -> SSD = 64 * (1.0^2) = 64.0
  TEST_FLOAT_CLOSE(ssd, 64.0f, 0.001f);

  printf("✓ 8x8 patch SSD tests passed\n");
}

int main(void) {
  printf("Running pre-estimation tests...\n");

  test_absolute_hearing_thresholds();
  test_critical_bands();
  test_masking_estimator();
  test_spectral_smoother();
  test_transient_detector();
  test_simd_utils();
  test_fast_exp();
  test_patch_ssd();

  printf("✅ All pre-estimation tests passed!\n");
  return 0;
}
