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

#include "specbleach_version.h"
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

static void test_version_string(void) {
  printf("Testing version string...\n");
  const char* version = specbleach_get_version_string();
  TEST_ASSERT(version != NULL, "Version string should not be NULL");
  TEST_ASSERT(strcmp(version, SPECBLEACH_VERSION_STRING) == 0,
              "Runtime string should match header macro");
}

static void test_version_components_compose_string(void) {
  printf("Testing version component composition...\n");
  char expected[32];
  snprintf(expected, sizeof(expected), "%u.%u.%u",
           (unsigned)SPECBLEACH_VERSION_MAJOR,
           (unsigned)SPECBLEACH_VERSION_MINOR,
           (unsigned)SPECBLEACH_VERSION_PATCH);
  TEST_ASSERT(strcmp(SPECBLEACH_VERSION_STRING, expected) == 0,
              "Version string should be composed of its components");
}

int main(void) {
  test_version_string();
  test_version_components_compose_string();

  printf("✅ All specbleach version tests passed!\n");
  return 0;
}
