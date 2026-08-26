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

#include "fft_transform.h"
#include "../utils/general_utils.h"

#include "pffft.h"
#include <stdlib.h>
#include <string.h>

static uint32_t calculate_fft_size(FftTransform* self);
static bool allocate_pffft(FftTransform* self);
static bool is_valid_fft_size(uint32_t n);
static uint32_t get_next_valid_fft_size(uint32_t n);

struct FftTransform {
  PFFFT_Setup* setup;

  uint32_t fft_size;
  uint32_t frame_size;
  uint32_t zeropadding_amount;
  uint32_t copy_position;
  ZeroPaddingType padding_type;
  uint32_t padding_amount;
  float* input_fft_buffer;
  float* output_fft_buffer;
  float* pffft_work_buffer;
  float* pffft_canonical_buffer;
};

static bool is_valid_fft_size(uint32_t n) {
  if (n < 32U || (n % 32U != 0U)) {
    return false;
  }
  uint32_t temp = n;
  while (temp % 2U == 0U) {
    temp /= 2U;
  }
  while (temp % 3U == 0U) {
    temp /= 3U;
  }
  while (temp % 5U == 0U) {
    temp /= 5U;
  }
  return (temp == 1U);
}

static uint32_t get_next_valid_fft_size(uint32_t n) {
  if (n < 32U) {
    return 32U;
  }
  uint32_t candidate = ((n + 31U) / 32U) * 32U;
  while (!is_valid_fft_size(candidate)) {
    candidate += 32U;
  }
  return candidate;
}

FftTransform* fft_transform_initialize(const uint32_t frame_size,
                                       const ZeroPaddingType padding_type,
                                       const uint32_t zeropadding_amount) {
  FftTransform* self = (FftTransform*)calloc(1U, sizeof(FftTransform));
  if (!self) {
    return NULL;
  }

  self->padding_type = padding_type;
  self->zeropadding_amount = zeropadding_amount;
  self->frame_size = frame_size;

  self->fft_size = calculate_fft_size(self);

  self->copy_position = (self->fft_size / 2U) - (self->frame_size / 2U);

  if (!allocate_pffft(self)) {
    fft_transform_free(self);
    return NULL;
  }

  return self;
}

FftTransform* fft_transform_initialize_bins(const uint32_t fft_size) {
  if (fft_size == 0U) {
    return NULL;
  }

  FftTransform* self = (FftTransform*)calloc(1U, sizeof(FftTransform));
  if (!self) {
    return NULL;
  }

  self->fft_size = get_next_valid_fft_size(fft_size);
  self->frame_size = self->fft_size;

  if (!allocate_pffft(self)) {
    fft_transform_free(self);
    return NULL;
  }

  return self;
}

static bool allocate_pffft(FftTransform* self) {
  self->setup = pffft_new_setup((int)self->fft_size, PFFFT_REAL);
  if (!self->setup) {
    return false;
  }

  self->input_fft_buffer =
      (float*)pffft_aligned_malloc(self->fft_size * sizeof(float));
  self->output_fft_buffer =
      (float*)pffft_aligned_malloc(self->fft_size * sizeof(float));
  self->pffft_work_buffer =
      (float*)pffft_aligned_malloc(self->fft_size * sizeof(float));
  self->pffft_canonical_buffer =
      (float*)pffft_aligned_malloc(self->fft_size * sizeof(float));

  if (!self->input_fft_buffer || !self->output_fft_buffer ||
      !self->pffft_work_buffer || !self->pffft_canonical_buffer) {
    return false;
  }

  memset(self->input_fft_buffer, 0, self->fft_size * sizeof(float));
  memset(self->output_fft_buffer, 0, self->fft_size * sizeof(float));
  memset(self->pffft_work_buffer, 0, self->fft_size * sizeof(float));
  memset(self->pffft_canonical_buffer, 0, self->fft_size * sizeof(float));

  return true;
}

static uint32_t calculate_fft_size(FftTransform* self) {
  switch (self->padding_type) {
    case NO_PADDING: {
      self->padding_amount = 0;
      return get_next_valid_fft_size(self->frame_size);
    }
    case NEXT_POWER_OF_TWO: {
      uint32_t next_power_of_two = get_next_power_two((int)self->frame_size);
      if (next_power_of_two < 32U) {
        next_power_of_two = 32U;
      }
      self->padding_amount = (next_power_of_two > self->frame_size)
                                 ? (next_power_of_two - self->frame_size)
                                 : 0U;
      return next_power_of_two;
    }
    case FIXED_AMOUNT: {
      self->padding_amount = self->zeropadding_amount;
      return get_next_valid_fft_size(self->frame_size + self->padding_amount);
    }
    default:
      return get_next_valid_fft_size(self->frame_size);
  }
}

void fft_transform_free(FftTransform* self) {
  if (!self) {
    return;
  }

  if (self->input_fft_buffer) {
    pffft_aligned_free(self->input_fft_buffer);
  }
  if (self->output_fft_buffer) {
    pffft_aligned_free(self->output_fft_buffer);
  }
  if (self->pffft_work_buffer) {
    pffft_aligned_free(self->pffft_work_buffer);
  }
  if (self->pffft_canonical_buffer) {
    pffft_aligned_free(self->pffft_canonical_buffer);
  }

  if (self->setup) {
    pffft_destroy_setup(self->setup);
  }

  free(self);
}

uint32_t get_fft_size(FftTransform* self) {
  if (!self) {
    return 0;
  }
  return self->fft_size;
}

uint32_t get_fft_real_spectrum_size(FftTransform* self) {
  if (!self) {
    return 0;
  }
  return (self->fft_size / 2U) + 1U;
}

bool fft_load_input_samples(FftTransform* self, const float* input) {
  if (!self || !input) {
    return false;
  }

  // Ensure buffer bounds are safe
  if (self->frame_size + self->copy_position > self->fft_size) {
    return false;
  }

  // Clear residual IFFT data from zero-padding regions
  memset(self->input_fft_buffer, 0, self->fft_size * sizeof(float));

  // Copy centered values only
  for (uint32_t i = self->copy_position;
       i < (self->frame_size + self->copy_position); i++) {
    self->input_fft_buffer[i] = input[i - self->copy_position];
  }

  return true;
}

bool fft_get_output_samples(FftTransform* self, float* output) {
  if (!self || !output) {
    return false;
  }

  // Ensure buffer bounds are safe
  if (self->frame_size + self->copy_position > self->fft_size) {
    return false;
  }

  // Copy centered values only
  for (uint32_t i = self->copy_position;
       i < (self->frame_size + self->copy_position); i++) {
    output[i - self->copy_position] = self->input_fft_buffer[i];
  }

  return true;
}

bool fft_accumulate_output_samples(FftTransform* self, float* accumulator) {
  if (!self || !accumulator) {
    return false;
  }

  if (self->frame_size + self->copy_position > self->fft_size) {
    return false;
  }

  for (uint32_t i = self->copy_position;
       i < (self->frame_size + self->copy_position); i++) {
    accumulator[i - self->copy_position] += self->input_fft_buffer[i];
  }

  return true;
}

bool compute_forward_fft(FftTransform* self) {
  if (!self || !self->setup) {
    return false;
  }

  // Perform ordered real forward transform
  pffft_transform_ordered(self->setup, self->input_fft_buffer,
                          self->pffft_canonical_buffer, self->pffft_work_buffer,
                          PFFFT_FORWARD);

  // Convert from PFFFT canonical layout to halfcomplex layout
  const uint32_t n = self->fft_size;
  const uint32_t n2 = n / 2U;

  self->output_fft_buffer[0] = self->pffft_canonical_buffer[0];
  self->output_fft_buffer[n2] = self->pffft_canonical_buffer[1];

  for (uint32_t k = 1U; k < n2; k++) {
    self->output_fft_buffer[k] = self->pffft_canonical_buffer[2U * k];
    self->output_fft_buffer[n - k] =
        self->pffft_canonical_buffer[(2U * k) + 1U];
  }

  return true;
}

bool compute_backward_fft(FftTransform* self) {
  if (!self || !self->setup) {
    return false;
  }

  // Convert from halfcomplex layout to PFFFT canonical layout
  const uint32_t n = self->fft_size;
  const uint32_t n2 = n / 2U;

  self->pffft_canonical_buffer[0] = self->output_fft_buffer[0];
  self->pffft_canonical_buffer[1] = self->output_fft_buffer[n2];

  for (uint32_t k = 1U; k < n2; k++) {
    self->pffft_canonical_buffer[2U * k] = self->output_fft_buffer[k];
    self->pffft_canonical_buffer[(2U * k) + 1U] =
        self->output_fft_buffer[n - k];
  }

  // Perform ordered real backward transform
  pffft_transform_ordered(self->setup, self->pffft_canonical_buffer,
                          self->input_fft_buffer, self->pffft_work_buffer,
                          PFFFT_BACKWARD);

  return true;
}

float* get_fft_input_buffer(FftTransform* self) {
  return self->input_fft_buffer;
}

float* get_fft_output_buffer(FftTransform* self) {
  return self->output_fft_buffer;
}
