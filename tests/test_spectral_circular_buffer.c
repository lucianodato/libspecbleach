/*
 * Unit tests for the SbSpectralCircularBuffer module
 */

#include "shared/utils/spectral_circular_buffer.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

#define TEST_FLOAT_CLOSE(a, b, tolerance)                                      \
  TEST_ASSERT(fabsf((a) - (b)) < (tolerance), "Float values not close enough")

void test_circular_buffer_init_and_free(void) {
  printf("Testing circular buffer initialization and cleanup...\n");

  uint32_t num_frames = 10;
  SbSpectralCircularBuffer* cb = spectral_circular_buffer_create(num_frames);

  TEST_ASSERT(cb != NULL, "Initialization should succeed");

  uint32_t layer0 = spectral_circular_buffer_add_layer(cb, 128);
  uint32_t layer1 = spectral_circular_buffer_add_layer(cb, 256);
  uint32_t layer2 = spectral_circular_buffer_add_layer(cb, 512);

  TEST_ASSERT(layer0 == 0, "Layer 0 ID check");
  TEST_ASSERT(layer1 == 1, "Layer 1 ID check");
  TEST_ASSERT(layer2 == 2, "Layer 2 ID check");

  spectral_circular_buffer_free(cb);
  printf("✓ Initialization and cleanup passed\n");
}

void test_circular_buffer_push_retrieve(void) {
  printf("Testing circular buffer push and retrieve...\n");

  uint32_t num_frames = 5;
  SbSpectralCircularBuffer* cb = spectral_circular_buffer_create(num_frames);

  uint32_t layer0 = spectral_circular_buffer_add_layer(cb, 4);
  uint32_t layer1 = spectral_circular_buffer_add_layer(cb, 4);

  float data1[] = {1.0f, 2.0f, 3.0f, 4.0f};
  float data2[] = {5.0f, 6.0f, 7.0f, 8.0f};

  // Push to layer 0 and retrieve with 0 delay (current frame)
  spectral_circular_buffer_push(cb, layer0, data1);
  float* retrieved = spectral_circular_buffer_retrieve(cb, layer0, 0);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(retrieved[i], data1[i], 1e-6f);
  }

  // Push to layer 1 and retrieve
  spectral_circular_buffer_push(cb, layer1, data2);
  retrieved = spectral_circular_buffer_retrieve(cb, layer1, 0);
  for (int i = 0; i < 4; i++) {
    TEST_FLOAT_CLOSE(retrieved[i], data2[i], 1e-6f);
  }

  spectral_circular_buffer_free(cb);
  printf("✓ Push and current frame retrieval passed\n");
}

void test_circular_buffer_delay_alignment(void) {
  printf("Testing circular buffer delay alignment (temporal wrap-around)...\n");

  uint32_t num_frames = 4;
  SbSpectralCircularBuffer* cb = spectral_circular_buffer_create(num_frames);
  uint32_t layer0 = spectral_circular_buffer_add_layer(cb, 1);

  // Buffer state progression:
  // Frame 0: push 10.0, advance
  // Frame 1: push 20.0, advance
  // Frame 2: push 30.0, advance
  // Frame 3: push 40.0, now write_index is 3

  float val = 10.0f;
  spectral_circular_buffer_push(cb, layer0, &val);
  spectral_circular_buffer_advance(cb);

  val = 20.0f;
  spectral_circular_buffer_push(cb, layer0, &val);
  spectral_circular_buffer_advance(cb);

  val = 30.0f;
  spectral_circular_buffer_push(cb, layer0, &val);
  spectral_circular_buffer_advance(cb);

  val = 40.0f;
  spectral_circular_buffer_push(cb, layer0, &val);
  // Current write_index is 3.

  // Retrieve current (delay 0): should be 40.0
  TEST_FLOAT_CLOSE(*spectral_circular_buffer_retrieve(cb, layer0, 0), 40.0f,
                   1e-6f);

  // Retrieve delay 1: should be 30.0
  TEST_FLOAT_CLOSE(*spectral_circular_buffer_retrieve(cb, layer0, 1), 30.0f,
                   1e-6f);

  // Retrieve delay 2: should be 20.0
  TEST_FLOAT_CLOSE(*spectral_circular_buffer_retrieve(cb, layer0, 2), 20.0f,
                   1e-6f);

  // Retrieve delay 3: should be 10.0
  TEST_FLOAT_CLOSE(*spectral_circular_buffer_retrieve(cb, layer0, 3), 10.0f,
                   1e-6f);

  // Push another to test circular wrap-around
  spectral_circular_buffer_advance(cb); // write_index -> 0
  val = 50.0f;
  spectral_circular_buffer_push(cb, layer0, &val);

  // Current (delay 0) is 50.0
  TEST_FLOAT_CLOSE(*spectral_circular_buffer_retrieve(cb, layer0, 0), 50.0f,
                   1e-6f);
  // Delay 1 is now 40.0 (previous write_index 3)
  TEST_FLOAT_CLOSE(*spectral_circular_buffer_retrieve(cb, layer0, 1), 40.0f,
                   1e-6f);

  spectral_circular_buffer_free(cb);
  printf("✓ Temporal delay alignment passed\n");
}

void test_circular_buffer_edge_cases(void) {
  printf(
      "Testing circular buffer edge cases (NULL checks, invalid indices)...\n");

  SbSpectralCircularBuffer* cb = spectral_circular_buffer_create(10);
  uint32_t layer0 = spectral_circular_buffer_add_layer(cb, 128);

  // Check NULL handling
  spectral_circular_buffer_push(NULL, layer0, NULL);
  TEST_ASSERT(spectral_circular_buffer_retrieve(NULL, layer0, 0) == NULL,
              "NULL handle retrieve");
  TEST_ASSERT(spectral_circular_buffer_retrieve(cb, 99, 0) == NULL,
              "Invalid layer id retrieve");
  spectral_circular_buffer_advance(NULL);

  spectral_circular_buffer_free(cb);
  spectral_circular_buffer_free(NULL);
  printf("✓ Edge cases passed\n");
}

int main(void) {
  printf("Running SbSpectralCircularBuffer tests...\n\n");

  test_circular_buffer_init_and_free();
  test_circular_buffer_push_retrieve();
  test_circular_buffer_delay_alignment();
  test_circular_buffer_edge_cases();

  printf("\n✅ All spectral circular buffer tests passed!\n");
  return 0;
}
