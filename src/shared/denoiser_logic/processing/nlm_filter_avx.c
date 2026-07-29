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

// Compile this file only when AVX is enabled globally, or when dynamic dispatch
// compilation flags enable AVX specifically for this translation unit.

#include "shared/denoiser_logic/processing/nlm_filter_internal.h"

bool nlm_filter_process_avx(NlmFilter* filter, float* smoothed_snr) {
  return nlm_filter_process_core(filter, smoothed_snr);
}
