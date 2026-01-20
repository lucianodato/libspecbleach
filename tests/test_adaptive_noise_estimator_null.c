#include "../src/shared/noise_estimation/adaptive_noise_estimator.h"
#include <assert.h>
#include <stdlib.h>

void test_louizou_free_null(void) {
  louizou_estimator_free(NULL);
}

void test_spp_mmse_free_null(void) {
  spp_mmse_estimator_free(NULL);
}

int main(void) {
  test_louizou_free_null();
  test_spp_mmse_free_null();
  return 0;
}
