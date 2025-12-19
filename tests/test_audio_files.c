/*
 * Audio file validation tests - ensure test data is accessible
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Test that essential test data files exist and are readable
void test_data_files_exist() {
  printf("Testing essential test data file accessibility...\n");

  // Check if Speech.wav exists (required)
  FILE* file = fopen("../tests/test_data/Speech.wav", "rb");
  if (!file) {
    fprintf(stderr, "ERROR: ../tests/test_data/Speech.wav not found\n");
    fprintf(stderr,
            "This file is required for audio testing. See tests/README.md for "
            "details.\n");
    exit(1);
  }
  fclose(file);

  // Check if reference denoised files exist (optional - can be regenerated)
  file = fopen("../tests/test_data/Speech_denoised.wav", "rb");
  if (!file) {
    printf("  Note: Speech_denoised.wav not found (can be regenerated)\n");
  } else {
    fclose(file);
  }

  file = fopen("../tests/test_data/Speech_adaptive_denoised.wav", "rb");
  if (!file) {
    printf(
        "  Note: Speech_adaptive_denoised.wav not found (can be "
        "regenerated)\n");
  } else {
    fclose(file);
  }

  printf("✓ Essential test data files are accessible\n");
}

// Test basic file I/O operations
void test_file_operations() {
  printf("Testing basic file operations...\n");

  const char* test_filename = "../tests/test_data/test_output.tmp";

  // Test writing to test directory
  FILE* file = fopen(test_filename, "wb");
  if (!file) {
    fprintf(stderr, "ERROR: Cannot write to test_data directory\n");
    exit(1);
  }

  // Write some test data
  const char* test_data = "test";
  fwrite(test_data, 1, 4, file);
  fclose(file);

  // Test reading back
  file = fopen(test_filename, "rb");
  if (!file) {
    fprintf(stderr, "ERROR: Cannot read test file\n");
    exit(1);
  }

  char buffer[5] = {0};
  (void)fread(buffer, 1, 4, file);
  fclose(file);

  if (strcmp(buffer, test_data) != 0) {
    fprintf(stderr, "ERROR: File read/write test failed\n");
    exit(1);
  }

  // Clean up
  remove(test_filename);

  printf("✓ File operations work correctly\n");
}

int main() {
  printf("Running audio file validation tests...\n\n");

  test_data_files_exist();
  test_file_operations();

  printf("\n✅ All audio file validation tests passed!\n");
  printf("Test data is ready for comprehensive audio testing.\n");

  return 0;
}
