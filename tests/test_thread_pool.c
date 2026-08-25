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

/*
 * Unit tests for SbThreadPool
 */

#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/utils/thread_pool.h"

#define TEST_ASSERT(condition, message)                                        \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "TEST FAILED: %s\n", message);                           \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

typedef struct {
  atomic_uint* hits; // Per-item execution counter (detects misses/duplicates)
  atomic_uint processed;
  uint32_t num_items;
} FillContext;

static void fill_task(void* arg, uint32_t start, uint32_t count) {
  FillContext* ctx = (FillContext*)arg;
  for (uint32_t i = start; i < start + count; i++) {
    atomic_fetch_add_explicit(&ctx->hits[i], 1U, memory_order_relaxed);
  }
  atomic_fetch_add_explicit(&ctx->processed, count, memory_order_relaxed);
}

static void fill_context_init(FillContext* ctx, uint32_t num_items) {
  ctx->num_items = num_items;
  ctx->hits = (atomic_uint*)calloc(num_items, sizeof(atomic_uint));
  TEST_ASSERT(ctx->hits != NULL, "Failed to allocate hits array");
  for (uint32_t i = 0; i < num_items; i++) {
    atomic_init(&ctx->hits[i], 0U);
  }
  atomic_init(&ctx->processed, 0U);
}

static void fill_context_verify_and_free(FillContext* ctx,
                                         const char* case_name) {
  TEST_ASSERT(atomic_load(&ctx->processed) == ctx->num_items, case_name);
  for (uint32_t i = 0; i < ctx->num_items; i++) {
    if (atomic_load(&ctx->hits[i]) != 1U) {
      fprintf(stderr, "  item %u hit %u times (%s)\n", i,
              atomic_load(&ctx->hits[i]), case_name);
      TEST_ASSERT(atomic_load(&ctx->hits[i]) == 1U, case_name);
    }
  }
  free(ctx->hits);
}

static void test_pool_null_fallback(void) {
  printf("Testing NULL pool sequential fallback...\n");

  FillContext ctx;
  fill_context_init(&ctx, 16U);
  sb_thread_pool_parallel_for(NULL, 16U, fill_task, &ctx);
  fill_context_verify_and_free(&ctx, "NULL pool should run sequentially");

  printf("✓ NULL pool sequential fallback tests passed\n");
}

static void test_pool_edge_cases(void) {
  printf("Testing pool edge cases...\n");

  TEST_ASSERT(sb_thread_pool_create(0U) == NULL,
              "Zero workers should be rejected");
  sb_thread_pool_free(NULL); // Must be a no-op

  FillContext ctx;

  fill_context_init(&ctx, 0U);
  SbThreadPool* pool = sb_thread_pool_create(4U);
  TEST_ASSERT(pool != NULL, "Pool creation failed");
  sb_thread_pool_parallel_for(pool, 0U, fill_task, &ctx);
  TEST_ASSERT(atomic_load(&ctx.processed) == 0U,
              "Zero items should process nothing");
  free(ctx.hits);

  fill_context_init(&ctx, 1U);
  sb_thread_pool_parallel_for(pool, 1U, fill_task, &ctx);
  fill_context_verify_and_free(&ctx, "Single item should process once");

  // Fewer items than participants
  fill_context_init(&ctx, 2U);
  sb_thread_pool_parallel_for(pool, 2U, fill_task, &ctx);
  fill_context_verify_and_free(&ctx, "Two items on five threads");

  sb_thread_pool_free(pool);

  printf("✓ Pool edge case tests passed\n");
}

static void test_pool_parallel_correctness(void) {
  printf("Testing pool parallel correctness...\n");

  SbThreadPool* pool = sb_thread_pool_create(3U);
  TEST_ASSERT(pool != NULL, "Pool creation failed");
  TEST_ASSERT(sb_thread_pool_num_workers(pool) == 3U, "Worker count mismatch");

  FillContext ctx;
  fill_context_init(&ctx, 1000U);
  sb_thread_pool_parallel_for(pool, 1000U, fill_task, &ctx);
  fill_context_verify_and_free(&ctx, "Parallel fill should cover all items");

  // Repeated dispatches reuse the same workers and work item storage
  for (uint32_t round = 0; round < 200U; round++) {
    fill_context_init(&ctx, 64U);
    sb_thread_pool_parallel_for(pool, 64U, fill_task, &ctx);
    fill_context_verify_and_free(&ctx, "Repeated dispatch should stay correct");
  }

  sb_thread_pool_free(pool);

  printf("✓ Pool parallel correctness tests passed\n");
}

int main(void) {
  test_pool_null_fallback();
  test_pool_edge_cases();
  test_pool_parallel_correctness();

  printf("\n✅ All thread pool tests passed!\n");
  return 0;
}
