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

#include "shared/utils/transient_detector.h"
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(cond, msg)                                                 \
  do {                                                                         \
    if (!(cond)) {                                                             \
      fprintf(stderr, "FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__);          \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

int main(void) {
  uint32_t num_bands = 24;
  TransientDetector* td = transient_detector_initialize(num_bands);
  TEST_ASSERT(td != NULL, "TransientDetector initialization should succeed");

  float steady_energies[24];
  float onset_weights[24];
  for (int i = 0; i < 24; i++) {
    steady_energies[i] = 1.0f;
  }

  // First frame initializes background energy
  float intensity = 0.0f;
  bool is_transient = transient_detector_process(td, steady_energies,
                                                 onset_weights, &intensity);
  TEST_ASSERT(!is_transient, "First frame should not trigger transient");
  TEST_ASSERT(intensity == 0.0f, "First frame intensity should be 0.0");

  // Sudden high-frequency energy burst (e.g. bands 18-20, sharp click/snap)
  float burst_energies[24];
  for (int i = 0; i < 24; i++) {
    burst_energies[i] = (i >= 18 && i <= 20) ? 10.0f : 1.0f;
  }

  is_transient =
      transient_detector_process(td, burst_energies, onset_weights, &intensity);
  TEST_ASSERT(is_transient, "Energy burst should trigger transient");
  TEST_ASSERT(onset_weights[18] > 0.8f,
              "High band onset weight should be high");
  TEST_ASSERT(intensity >= 0.25f,
              "Sharp high-frequency burst should produce strong intensity");

  // Broadband transient burst across all bands
  float broadband_energies[24];
  for (int i = 0; i < 24; i++) {
    broadband_energies[i] = 10.0f;
  }
  is_transient = transient_detector_process(td, broadband_energies,
                                            onset_weights, &intensity);
  TEST_ASSERT(is_transient, "Broadband burst should trigger transient");
  TEST_ASSERT(intensity >= 0.8f,
              "Broadband burst should produce high intensity");

  transient_detector_free(td);
  printf("✓ Transient Detector tests passed\n");
  return 0;
}
